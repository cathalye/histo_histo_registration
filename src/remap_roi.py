import sys

import numpy as np
import SimpleITK as sitk

from src.histology_data import HistologyData

# https://github.com/pyushkevich/histoannot.git
# git clone the GitHub repository and add the path to the sys.path
# This makes sure you are using the latest version of the code
sys.path.append('/Users/cathalye/Packages/histoannot/')
import phas.client.api as phas
from phas.dltrain import spatial_transform_roi

class RemapROI:

    def __init__(self, registration_dir, moving_slide: HistologyData):
        self.registration_dir = registration_dir
        self.chunk_mask = sitk.ReadImage(f"{registration_dir}/reference_chunk_mask.nii.gz")
        self.moving_slide = moving_slide
        self._get_nearest_chunk_map(self.chunk_mask)
        self.moving_slide_single_channel = self.moving_slide.get_single_channel_image(channel=1)


    def _get_nearest_chunk_map(self, chunk_mask):
        # Explanation of this function in nearest_chunk_map.ipynb
        chunk_mask_arr = sitk.GetArrayFromImage(chunk_mask)[0, :, :]

        chunk_labels = np.unique(chunk_mask_arr)
        chunk_labels = chunk_labels[chunk_labels != 0] # Remove the background label

        def _get_dist_from_chunk(k):
            # Get the distance of every pixel from the boundary of the chunk with label k
            mask = sitk.BinaryThreshold(chunk_mask, int(k), int(k), 1, 0) # Extract only the chunk with label k
            return sitk.SignedDanielssonDistanceMap(mask, insideIsPositive=False, squaredDistance=True)

        dist_maps_all_chunks = { k: _get_dist_from_chunk(k) for k in chunk_labels }
        dist_maps_all_chunks = np.array([ dist_maps_all_chunks[k] for k in chunk_labels ])

        # Find the chunk with the minimum distance to a given point
        nearest_chunk = chunk_labels[np.argmin(dist_maps_all_chunks, axis = 0)]

        self.nearest_chunk_map = np.reshape(nearest_chunk, chunk_mask_arr.shape)


    def get_chunk_transforms(self, x, y):
        chunk = self.nearest_chunk_map[(y, x)]
        chunk_str = f"{chunk:02d}"
        chunk_rigid = np.loadtxt(f"{self.registration_dir}/transforms/piecewise_rigid_{chunk_str}.mat")
        chunk_warp = sitk.ReadImage(f"{self.registration_dir}/transforms/piecewise_deformable_{chunk_str}.nii.gz")

        return chunk_warp, chunk_rigid


    def registration_transform(self, xy):
        # Apply the transformations to the sampling ROI in the opposite order
        # i.e first the piecewise deformable, then the piecewise rigid
        # All transformations are applied in the physical space]

        xy = np.array(xy)
        chunk_warp, chunk_rigid = self.get_chunk_transforms(int(xy[0]), int(xy[1]))

        xy_phys = np.array(chunk_warp.TransformContinuousIndexToPhysicalPoint(xy))

        # Step 1: Evaluate the deformable transform
        displacement = chunk_warp.EvaluateAtContinuousIndex(xy_phys)
        xy_warp = xy_phys + displacement

        # Step 2: Evaluate the piecewise rigid transform
        # The piecewise rigid transform also considers the global rigid transform
        # so no need to apply that separately
        A = chunk_rigid[:2, :2]
        b = chunk_rigid[:2, 2]

        # NOTE:
        # sitk uses LPS coordinate system, and c3d/greedy uses RAS
        # xy_warp is in LPS coordinate system.
        # To go from LPS to RAS the transformation matrix
        # is [[-1, 0], [0, -1]] for a 2D image.
        #
        # X_warp_RAS = -X_warp_LPS
        # X_rigid_RAS = A @ X_warp_RAS + b
        # X_rigid_LPS = -X_rigid_RAS = -A @ X_warp_RAS - b = A @ X_warp_LPS - b
        #
        # So to skip the back and forth conversion between LPS and RAS, we can
        # directly apply A @ X_warp - b to get the coordinates in LPS
        xy_chunk_rigid = A @ xy_warp - b

        # Get coordinates from physical space to index space in the moving image
        xy_remap = np.array(self.moving_slide_single_channel.TransformPhysicalPointToContinuousIndex(xy_chunk_rigid))

        return xy_remap[0], xy_remap[1]


def process_roi_data(roi_data, type, scale=1):
    # Convert the ROI coordinates to from full resolution to thumbnail space
    # Return it in the format that PHAS slie_sampling_roi expects
    if type == 'polygon':
        x_roi, y_roi = zip(*roi_data)
        assert len(x_roi) == len(y_roi), "x_roi and y_roi must have the same length"
        x_roi = np.array(x_roi)
        y_roi = np.array(y_roi)

        x_scaled_roi = (x_roi + 0.5) * scale
        y_scaled_roi = (y_roi + 0.5) * scale

        return [list(pair) for pair in zip(x_scaled_roi, y_scaled_roi)]

    elif type == 'trapezoid':
        x_roi, y_roi, w_roi = zip(*roi_data)
        assert len(x_roi) == len(y_roi) == len(w_roi), "x_roi, y_roi, and w_roi must have the same length"
        x_roi = np.array(x_roi)
        y_roi = np.array(y_roi)
        w_roi = np.array(w_roi)

        x_scaled_roi = (x_roi + 0.5) * scale
        y_scaled_roi = (y_roi + 0.5) * scale
        w_scaled_roi = w_roi * scale

        return [list(pair) for pair in zip(x_scaled_roi, y_scaled_roi, w_scaled_roi)]
    else:
        raise ValueError(f"Unknown ROI type: {type}")
