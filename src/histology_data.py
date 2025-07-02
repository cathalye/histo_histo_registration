"""
This module contains the HistologyData class, which is used to load and process
histology data.

Inputs
(1) Slide ID

Processes
(1) Get the scaling factor between full resolution slide and thumbnail
(2) Get the binary mask of the slide
(3) Get the chunk mask of the slide
(4) Get the nearest chunk map of the slide

Outputs


"""

import subprocess
import sys

import numpy as np
import SimpleITK as sitk

# https://github.com/pyushkevich/histoannot.git
# git clone the GitHub repository and add the path to the sys.path
# This makes sure you are using the latest version of the code
sys.path.append('/Users/cathalye/Packages/histoannot/')
import phas.client.api as phas
from phas.dltrain import spatial_transform_roi


class HistologyData:

    def __init__(self, task: phas.SamplingROITask, slide_id, thumbnail_path=None):
        self.slide = phas.Slide(task, slide_id)
        self.full_size = self.slide.dimensions
        if thumbnail_path is None:
            self.thumbnail = None
        else:
            self.thumbnail = sitk.ReadImage(thumbnail_path)[:, :, 0]
            self.thumbnail_size = self.thumbnail.GetSize()

        self.get_scaling_factor()


    def get_scaling_factor(self):
        if self.thumbnail is None:
            self.scaling_factor = 1
        # When downsampling, the largest dimension is scaled to 1000 pixels
        elif self.full_size[0] > self.full_size[1]:
            self.scaling_factor = 1000 / self.full_size[0]
        else:
            self.scaling_factor = 1000 / self.full_size[1]


    def get_thumbnail_coord_from_full(self, x, y):
        # Integer coordinates in the thumbnail are in the center of the voxel
        # Hence the +0.5 for adjustment
        x_thumbnail = (x + 0.5) * self.scaling_factor
        y_thumbnail = (y + 0.5) * self.scaling_factor

        return x_thumbnail, y_thumbnail


    def get_full_coord_from_thumbnail(self, x, y):
        x_full = (x + 0.5) / self.scaling_factor
        y_full = (y + 0.5) / self.scaling_factor

        return x_full, y_full


    def get_single_channel_image(self, channel=0):
        single_channel_image = sitk.VectorIndexSelectionCast(self.thumbnail, channel)
        return single_channel_image


    def _remove_mask_border(self, mask, border=25):
        # A lot of slides have shadow/dark artefacts on the border that affect
        # the mask thresholding. We remove the border of the mask to avoid this.
        mask_arr = sitk.GetArrayFromImage(mask)[0, :, :]
        mask_arr[:border, :] = 0
        mask_arr[-border:, :] = 0
        mask_arr[:, :border] = 0
        mask_arr[:, -border:] = 0

        return sitk.GetImageFromArray(mask_arr)


    def get_binary_mask(self, channel=1):
        # get the single channel image
        single_channel = self.get_single_channel_image(channel=channel)

        otsu_filter = sitk.OtsuThresholdImageFilter()
        otsu_filter.SetInsideValue(1)
        otsu_filter.SetOutsideValue(0)
        binary_mask = otsu_filter.Execute(single_channel)

        # Otsu thresholding has holes in the mask due to the irregular staining
        # Dilation to fill the holes
        structuring_element = sitk.sitkBall
        # XXX: hard coded radius
        radius = [2, 2, 2]
        closed_mask = sitk.BinaryMorphologicalClosing(binary_mask, radius, structuring_element)

        closed_mask = self._remove_mask_border(closed_mask)

        return closed_mask


    def get_chunk_mask(self, binary_mask_path, chunk_mask_path, n_parts=10):
        # BUG: There are no registered IO factories
        # RuntimeError: /Users/runner/work/image-graph-cut/image-graph-cut/be/install/include/ITK-5.4/itkImageFileReader.hxx:135:
        #  Could not create IO object for reading file /Users/cathalye/Projects/proj_histo_mri_greedy_registration/scratch/binary.nii.gz
        #  There are no registered IO factories.
        #  Please visit https://www.itk.org/Wiki/ITK/FAQ#NoFactoryException to diagnose the problem.

        # XXX: hard coded parameters
        # image_graph_cut(
        #     fn_input=binary_mask_path,
        #     fn_output=chunk_mask_path,
        #     n_parts=n_parts,
        #     tolerance=1.2,
        #     n_metis_iter=100,
        #     max_comp=4,
        #     min_comp_frac=0.1
        # )

        subprocess.run([
            "image_graph_cut",
            "-u", "1.2",
            "-n", "100",
            "-c", "4", "0.1",
            binary_mask_path,
            chunk_mask_path,
            str(n_parts),
        ], stdout=subprocess.DEVNULL)
