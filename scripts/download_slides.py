import argparse
import json
import os
import sys

import pandas as pd

# https://github.com/pyushkevich/histoannot.git
# git clone the GitHub repository and add the path to the sys.path
# This makes sure you are using the latest version of the code
sys.path.append('/Users/cathalye/Packages/histoannot/')
import phas.client.api as phas

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--server", type=str, required=True, help="The server to use")
    parser.add_argument("--private_key", type=str, required=True, help="The path to server private key")
    parser.add_argument("--task_id", type=int, required=True, help="The task id on the server")
    parser.add_argument("--ref_stain", type=str, required=True, help="Reference stain with the annotation")
    parser.add_argument("--mov_stain", type=str, required=True, help="Moving stain to be registered")
    parser.add_argument("--root_dir", type=str, required=True, help="The root directory to download the data")

    args = parser.parse_args()

    task_id = args.task_id
    ref_stain = args.ref_stain
    mov_stain = args.mov_stain
    root_dir = args.root_dir

    conn = phas.Client(args.server, args.private_key)
    task = phas.Task(conn, task_id)

    manifest = pd.DataFrame(task.manifest())
    # Get all the reference slides with ROIs
    manifest_ref_stain = manifest[(manifest['n_sampling_rois'] != 0) & (manifest['stain']==ref_stain)]

    for index, row in manifest_ref_stain.iterrows():
        specimen = row["specimen_private"]
        os.makedirs(f"{root_dir}/{specimen}", exist_ok=True)
        block = row["block_name"]
        os.makedirs(f"{root_dir}/{specimen}/{block}", exist_ok=True)
        section = row["section"]

        mov_slide_id = None
        ref_slide_id = None
        try:
            # Check if there is a mov_stain slide for the same specimen, block and section as the ref_stain slide
            mov_slide_id = manifest[
                (manifest["specimen_private"] == specimen)
                & (manifest["block_name"] == block)
                & (manifest["section"] == section)
                & (manifest["stain"] == mov_stain)
            ]["id"].values[0]
            ref_slide_id = row["id"]
        except IndexError:
            print(f"Moving slide not found for Specimen: {specimen}, Block: {block}, Section: {section}")

        # Make the stain name a valid filename without spaces or dashes
        ref_stain_str = ref_stain.lower().replace(" ", "").replace("-", "")
        mov_stain_str = mov_stain.lower().replace(" ", "").replace("-", "")

        ref_fname = f"{root_dir}/{specimen}/{block}/{ref_stain_str}_slide_thumbnail.nii.gz"
        mov_fname = f"{root_dir}/{specimen}/{block}/{mov_stain_str}_slide_thumbnail.nii.gz"
        # Download the slides if they don't exist
        if not os.path.exists(ref_fname):
            ref_slide = phas.Slide(task, ref_slide_id)
            ref_slide.thumbnail_nifti_image(ref_fname)
        if not os.path.exists(mov_fname):
            mov_slide = phas.Slide(task, mov_slide_id)
            mov_slide.thumbnail_nifti_image(mov_fname)

        # Download the ROI for the reference slide
        ROI_task = phas.SamplingROITask(conn, task_id)
        ROI_fname = f"{root_dir}/{specimen}/{block}/{ref_stain_str}_sampling_roi.nii.gz"
        if not os.path.exists(ROI_fname):
            try:
                ROI_task.slide_sampling_roi_nifti_image(ref_slide_id, ROI_fname)
            except Exception as e:
                print(f"Error downloading ROI for Reference slide {ref_slide_id}")

        # Save the block metadata as json file for future reference
        block_dict = {
            "specimen": specimen,
            "block": block,
            "section": section,
            "ref_slide_id": int(ref_slide_id),
            "mov_slide_id": int(mov_slide_id),
        }
        with open(f"{root_dir}/{specimen}/{block}/block_info.json", "w", encoding="utf-8") as f:
            json.dump(block_dict, f, indent=4)
