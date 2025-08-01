"""
This script remaps PHAS sampling ROI from one stain to another. The requirement
is that the two slides are registered and the registration transform files are
available.

Inputs:
(1) PHAS server URL, private key, and task ID
(2) Fixed slide ID (stain with the sampling ROI)
(3) Moving slide ID (stain to remap to)
(4) Transform files directory path (where the registration transform files are stored)
Assumes that the registration transform files are named as follows:
- WORK_DIR/transforms/piecewise_rigid_00.mat
- WORK_DIR/transforms/piecewise_deformable_00.nii.gz
- WORK_DIR/reference_chunk_mask.nii.gz

Process:
(1) Get the sampling ROI coordinates in thumbnail resolution
(2) Find which chunk the ROI coordinates lies in
(3) Apply the respective chunk transforms to the ROI coordinates
(4) Get the ROI coordinates in the moving slide thumbnail
(5) Upscale the ROI coordinates to full resolution
(6) Create a new sampling ROI on the moving slide on PHAS

Outputs:
(1) Sampling ROI on the moving slide on PHAS
"""

import argparse
import json
import os
import sys

# https://github.com/pyushkevich/histoannot.git
# git clone the GitHub repository and add the path to the sys.path
# This makes sure you are using the latest version of the code
sys.path.append('/Users/cathalye/Packages/histoannot/')
import phas.client.api as phas
from phas.dltrain import spatial_transform_roi

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from src.histology_data import HistologyData
from src.remap_roi import RemapROI, process_roi_data


def connect_to_server(task_id, phas_url, private_key):
    conn = phas.Client(phas_url, private_key, verify=False)
    # Create a sampling ROI task object to pass to Slide class for downloading sampling ROI json
    task = phas.SamplingROITask(conn, task_id)

    return task


if __name__ == '__main__':
    parse = argparse.ArgumentParser(description="Map PHAS sampling ROI from one stain to another")
    parse.add_argument('--phas_url', type=str, help='PHAS server URL')
    parse.add_argument('--private_key', type=str, help='Path to PHAS private key')
    parse.add_argument('--task_id', type=int, help='Task ID on the PHAS server')
    parse.add_argument('--fixed_slide_id', type=str, help='Slide ID of the fixed slide on PHAS')
    parse.add_argument('--moving_slide_id', type=str, help='Slide ID of the moving slide on PHAS')
    parse.add_argument('--moving_slide_thumbnail_path', type=str, help='Path to the thumbnail of the moving slide')
    parse.add_argument('--registration_dir', type=str, help='Path to the directory containing the registration files')
    args = parse.parse_args()

    # Parse the arguments to get the task ID, specimen, block, and slide IDs
    task = connect_to_server(args.task_id, args.phas_url, args.private_key)
    fixed_slide_id = args.fixed_slide_id
    moving_slide_id = args.moving_slide_id
    registration_dir = args.registration_dir

    fixed_slide = HistologyData(task, fixed_slide_id, thumbnail_path=None)
    moving_slide = HistologyData(task, moving_slide_id, thumbnail_path=args.moving_slide_thumbnail_path)

    # Create the Remap object
    remap = RemapROI(registration_dir, moving_slide)

    # Reset the sampling ROI on the moving slide
    task.delete_sampling_rois_on_slide(moving_slide_id)

    # Get the sampling ROI coordinates
    for roi in task.slide_sampling_rois(fixed_slide_id):
        print(f"Processing ROI {roi['id']}")

        roi_json = json.loads(roi['json'])
        roi_type = roi_json['type'] # polygon or trapezoid
        # If type is polygon, the data [[x1, y1], [x2, y2], ...]
        # If type is trapezoid, the data [[x1, y1, w1], [x2, y2, w2], ...]
        roi_data = roi_json['data'] # list of listsroi_data = roi_json['data'] # list of lists

        # Get the ROI coordinates in the thumbnail space
        roi_thumbnail_json = roi_json.copy()
        roi_thumbnail_json['data'] = process_roi_data(roi_data, roi_type, fixed_slide.scaling_factor)

        # Apply the remapping to the ROI
        roi_warped = spatial_transform_roi(roi_thumbnail_json, remap.registration_transform)

        # Get the ROI coordinates in the full resolution space for the tau slide
        roi_fullres_warped = process_roi_data(roi_warped['data'], roi_type, 1/moving_slide.scaling_factor)

        # Create the json for the tau slide
        roi_warped = roi_json.copy()
        roi_warped['data'] = roi_fullres_warped
        roi_moving = task.create_sampling_roi(moving_slide_id, roi['label'], roi_warped)

        print(f"Created ROI on the moving slide")
