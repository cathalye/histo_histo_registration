import argparse
import subprocess
import sys
import warnings
warnings.filterwarnings("ignore")

import json
import numpy as np
import SimpleITK as sitk


# https://github.com/pyushkevich/histoannot.git
# git clone the GitHub repository and add the path to the sys.path
# This makes sure you are using the latest version of the code
sys.path.append('/Users/cathalye/Packages/histoannot/')
import phas.client.api as phas
from phas.dltrain import spatial_transform_roi



def connect_to_server(task):
    # PHAS connection parameters
    PHAS_URL='https://histo.itksnap.org'
    PRIVATE_KEY = '/Users/cathalye/.private/histoitk_api_key.json'

    conn = phas.Client(PHAS_URL, PRIVATE_KEY)
    # Create a sampling ROI task object to pass to Slide class for downloading sampling ROI json
    task = phas.SamplingROITask(conn, task)

    return task

def read_json_property(jsonfile, property):
    with open(jsonfile, 'r') as file:
        data = json.load(file)
        slide_id = data.get(f'{property}', None)
    return slide_id

def find_scaling_factor(task, slide_id):
    # TODO: Read the spacing of slide and thumbnail and calculate the scaling factor using that
    # Find the scaling factor given a PHAS server slide ID
    # When downsampling, the largest dimension is scaled to 1000 pixels
    x, y = phas.Slide(task, slide_id).dimensions
    if x > y:
        scale = 1000 / x
    else:
        scale = 1000 / y
    return scale

def process_roi_data(roi, scale):
    # Load the json data from the ROI
    # TODO: How to handle this case elegantly?
    try:
        roi = json.loads(roi['json'])
    except:
        pass
    # It should contain the following fields:
    #   'type': polygon or trapezium dtype = string
    #   'data': [[x1, y1], [x2, y2], ...] dtype = list of lists
    x_roi, y_roi = zip(*roi['data'])
    assert len(x_roi) == len(y_roi), "x_roi and y_roi must have the same length"

    # Convert to numpy array to do element-wise operations
    x_roi = np.array(x_roi)
    y_roi = np.array(y_roi)

    x_scaled_roi = (x_roi + 0.5) * scale
    y_scaled_roi = (y_roi + 0.5) * scale

    # Return in the same format as input - list of lists
    xy_scaled = [list(pair) for pair in zip(x_scaled_roi, y_scaled_roi)]

    roi_scaled = roi.copy()
    roi_scaled['data'] = xy_scaled

    return roi_scaled


def get_chunk_id(registration_dir, x, y):
    # Read the multi chunk mask from the registration directory
    chunk_mask = sitk.ReadImage(f"{registration_dir}/reference_multi_chunk.nii.gz")
    chunk_mask = chunk_mask[:, :,0]
    # sitk.GetArrayFromImage returns a numpy array with flipped coordinates
    chunk_mask_arr = sitk.GetArrayFromImage(chunk_mask)

    chunk_ids = [ x for x in np.unique(chunk_mask) if x != 0 ]

    def chunk_dist_map(k):
        # Extract only the chunk with label k
        mask = sitk.BinaryThreshold(chunk_mask, k-0.5, k+0.5, 1, 0)
        return sitk.SignedDanielssonDistanceMap(mask, insideIsPositive=False, squaredDistance=True)

    # If a point lies outside the mask, then use the nearest chunk
    if chunk_mask_arr[(y, x)] == 0:
        # XXX: understand the logic here
        dist_cmask = { k: chunk_dist_map(k) for k in chunk_ids }
        # changed ((y,x)) see if it works
        dist = np.array([ dist_cmask[k].GetPixel((y,x)) for k in chunk_ids ])
        chunk = chunk_ids[np.argmin(dist)]
    else: # Check value of the label at given location
        chunk = chunk_mask_arr[(y, x)]
    chunk_str = f"{chunk:02d}"
    chunk_rigid = np.loadtxt(f"{registration_dir}/output_piecewise_rigid_{chunk_str}.mat")
    chunk_warp = sitk.ReadImage(f"{registration_dir}/output_piecewise_deformable_{chunk_str}.nii.gz")

    return chunk_warp, chunk_rigid

def map_sampling_roi(task_id, root_dir, block_json):
    task = connect_to_server(task_id)
    # Read from json
    specimen = read_json_property(block_json, 'specimen')
    block = read_json_property(block_json, 'block')
    nissl_slide_id = read_json_property(block_json, 'nissl_slide_id')
    tau_slide_id = read_json_property(block_json, 'tau_slide_id')

    tau_slide_thumbnail = sitk.ReadImage(f"{root_dir}/{specimen}/{block}/tau_slide_thumbnail.nii.gz")[:,:,0]
    print(f"Read thumbnail for slide {tau_slide_id} {tau_slide_thumbnail.GetSize()}")

    # Registration files
    registration_dir = f"{root_dir}/{specimen}/{block}/registration"
    global_rigid = np.loadtxt(f"{registration_dir}/output_global_rigid.mat")

    nissl_scaling = find_scaling_factor(task, nissl_slide_id)
    tau_scaling = find_scaling_factor(task, tau_slide_id)

    def my_transform(roi_xy):
        # Apply the transformations to the sampling ROI in the opposite order
        # i.e first the piecewise deformable, then the piecewise rigid
        # All transformations are applied in the physical space
        xy = np.array(roi_xy)
        # Determing which chunk containes the given point
        chunk_warp, chunk_rigid = get_chunk_id(registration_dir, int(xy[0]), int(xy[1]))
        # Deformable transform is saved as dispalcement field in the physical space
        # Get the coordinates in physical space
        xy_phys = np.array(chunk_warp.TransformContinuousIndexToPhysicalPoint(xy))
        # Step 1 - Piecewise deformable transform
        displacement = chunk_warp.EvaluateAtContinuousIndex(xy)
        xy_warp = xy_phys + displacement
        # Before the next step, there are some transformations to be done
        # XXX: Why are we flipping the coordinates? ITK-SNAP says the orientation
        # is RAI and shows all physical space coordinates as negative
        # xy_warp = -xy_warp # transformation to RAS??
        # Homogeneous coordinates for matrix multiplication

        # TODO: Instead of homogeneous do Ax + b
        A = chunk_rigid[:2, :2]
        b = chunk_rigid[:2, 2]
        # Step 2 - Piecewise rigid transform
        # xy_chunk_rigid = np.dot(chunk_rigid, xy_warp_homogeneous) # [x',y',1]

        xy_chunk_rigid = np.dot(A, xy_warp) - b # -b to ensure all coordinates are in RAS

        # Get coordinates from physical space to index space in the moving image
        xy_remap = np.array(tau_slide_thumbnail.TransformPhysicalPointToContinuousIndex(xy_chunk_rigid))
        # xy_remap = np.array(chunk_warp.TransformPhysicalPointToIndex(xy_chunk_rigid[0], xy_chunk_rigid[1]))

        return xy_remap[0], xy_remap[1]


    # Get the sampling ROI coordinates
    for roi in task.slide_sampling_rois(nissl_slide_id):
        roi_scaled = process_roi_data(roi, nissl_scaling)
        print("Done processing ROI data")
        roi_warped = spatial_transform_roi(roi_scaled, my_transform)
        print("Done warping ROI")
        roi_fullres_warped = process_roi_data(roi_warped, 1/tau_scaling) # we are scaling up hence 1/scale

        roi_id = task.create_sampling_roi(tau_slide_id, roi['label'], roi_fullres_warped)
        print(f"Created ROI {roi['label']} for slide {tau_slide_id}")


if __name__ == '__main__':
    parse = argparse.ArgumentParser(description="Map sampling ROI from NISSL to another stain")
    parse.add_argument('--task_id', type=int, default=39, help='Task ID on the PHAS server')
    parse.add_argument('--root_dir', type=str, default='./data/', help='Root directory of the data')
    parse.add_argument('--block_json', type=str, help='Path to block json file')
    args = parse.parse_args()

    map_sampling_roi(args.task_id, args.root_dir, args.block_json)
