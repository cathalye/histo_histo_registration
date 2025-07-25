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

sys.path.append('../src')
from histology_data import HistologyData
from remap_roi import RemapROI
import utils


if __name__ == '__main__':
    parse = argparse.ArgumentParser(description="Map sampling ROI from NISSL to another stain")
    parse.add_argument('--task_id', type=int, help='Task ID on the PHAS server')
    parse.add_argument('--root_dir', type=str, default='./data/', help='Root directory of the data')
    parse.add_argument('--block_json', type=str, help='Path to block json file')
    args = parse.parse_args()

    # Parse the arguments to get the task ID, specimen, block, and slide IDs
    task = utils.connect_to_server(args.task_id)
    specimen = utils.read_json_property(args.block_json, 'specimen')
    block = utils.read_json_property(args.block_json, 'block')
    nissl_slide_id = utils.read_json_property(args.block_json, 'nissl_slide_id')
    tau_slide_id = utils.read_json_property(args.block_json, 'tau_slide_id')

    print(f"Processing {specimen} {block}")

    if not os.path.exists(f"{args.root_dir}/{specimen}/{block}/result_workspace.itksnap"):
        print(f"Skipping {specimen} {block}. Registration workspace not found.")
        exit()
    else:
        registration_dir = f"{args.root_dir}/{specimen}/{block}/registration"

    # Get the thumbnail paths
    nissl_thumbnail_path = f"{args.root_dir}/{specimen}/{block}/nissl_slide_thumbnail.nii.gz"
    tau_thumbnail_path = f"{args.root_dir}/{specimen}/{block}/tau_slide_thumbnail.nii.gz"

    # Create the HistologyData objects for the Nissl and Tau slides
    nissl_slide = HistologyData(task, nissl_slide_id, nissl_thumbnail_path)
    tau_slide = HistologyData(task, tau_slide_id, tau_thumbnail_path)

    # Create the Remap object
    remap = RemapROI(registration_dir, tau_slide)

    # Reset the sampling ROI on the Tau slide
    utils.reset_slide(task, tau_slide_id)

    # Get the sampling ROI coordinates
    for roi in task.slide_sampling_rois(nissl_slide_id):
        print(f"Processing ROI {roi['id']}")

        roi_json = json.loads(roi['json'])
        roi_type = roi_json['type'] # polygon or trapezoid
        # If type is polygon, the data [[x1, y1], [x2, y2], ...]
        # If type is trapezoid, the data [[x1, y1, w1], [x2, y2, w2], ...]
        roi_data = roi_json['data'] # list of lists

        # Get the ROI coordinates in the thumbnail space
        roi_thumbnail_json = roi_json.copy()
        roi_thumbnail_json['data'] = utils.process_roi_data(roi_data, roi_type, nissl_slide.scaling_factor)

        # Apply the remapping to the ROI
        roi_warped = spatial_transform_roi(roi_thumbnail_json, remap.registration_transform)

        # Get the ROI coordinates in the full resolution space for the tau slide
        roi_fullres_warped = utils.process_roi_data(roi_warped['data'], roi_type, 1/tau_slide.scaling_factor)

        # Create the json for the tau slide
        roi_warped = roi_json.copy()
        roi_warped['data'] = roi_fullres_warped
        roi_tau = task.create_sampling_roi(tau_slide_id, roi['label'], roi_warped)

        print(f"Created ROI on the tau slide")
