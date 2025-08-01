"""
This script registers two histology slides using greedy.

Inputs:
(1) Reference slide thumbnail path
(2) Moving slide thumbnail path
(3) Working directory path (where the output files will be stored)

Process:
(1) Get scalar images of the reference and moving slides
(2) Get reference binary and chunk masks
(3) Global rigid registration
(4) Piecewise rigid registration
(5) Piecewise deformable registration
(6) Apply the transforms to the moving slide

Outputs:
(1) Transformed moving slide
(2) ITK-SNAP workspace to check the registration results
(3) Intermediate transform files, scalar images, and masks
"""

import argparse
import os
import subprocess
import sys

from picsl_greedy import Greedy2D, MultiChunkGreedy2D

greedy = Greedy2D()
multi_chunk_greedy = MultiChunkGreedy2D()

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from src.histology_data import HistologyData

if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="Register two histology slides")
    parse.add_argument('--reference_slide', type=str, help='Path to reference slide thumbnail')
    parse.add_argument('--moving_slide', type=str, help='Path to moving slide thumbnail')
    parse.add_argument('--working_dir', type=str, default='../work/', help='Working directory to store output files')
    args = parse.parse_args()

    reference_slide_path = args.reference_slide
    moving_slide_path = args.moving_slide
    working_dir = args.working_dir
    os.makedirs(working_dir, exist_ok=True)
    transforms_dir = os.path.join(working_dir, "transforms")
    os.makedirs(transforms_dir, exist_ok=True)

    reference_slide = HistologyData(task=None, slide_id=None, thumbnail_path=reference_slide_path)
    moving_slide = HistologyData(task=None, slide_id=None, thumbnail_path=moving_slide_path)

    # Step 1 - Get scalar images of the reference and moving slides
    reference_scalar_path = os.path.join(working_dir, "reference_scalar.nii.gz")
    reference_slide.get_single_channel_image(channel=1, save_path=reference_scalar_path, return_img=False)

    moving_scalar_path = os.path.join(working_dir, "moving_scalar.nii.gz")
    moving_slide.get_single_channel_image(channel=1, save_path=moving_scalar_path, return_img=False)

    # Step 2 - Get reference binary and chunk masks
    reference_binary_mask_path = os.path.join(working_dir, "reference_binary_mask.nii.gz")
    reference_slide.get_binary_mask(save_path=reference_binary_mask_path, return_img=False)

    reference_chunk_mask_path = os.path.join(working_dir, "reference_chunk_mask.nii.gz")
    reference_slide.get_chunk_mask(reference_binary_mask_path, reference_chunk_mask_path)


    # Step 3 - Global rigid registration
    global_rigid_path = os.path.join(transforms_dir, "global_rigid.mat")
    if not os.path.exists(global_rigid_path):
        greedy.execute('-d 2 -a -dof 6 '
                       '-i {} {} '
                       '-gm {} '
                       '-ia-image-centers '
                       '-m WNCC 4x4 '
                       '-wncc-mask-dilate '
                       '-n 200x200x40x0x0 '
                       '-search 20000 any 10 '
                       '-o {} '.format(reference_scalar_path, moving_scalar_path,
                                       reference_binary_mask_path,
                                       global_rigid_path))


    # Step 4 - Piecewise rigid registration
    piecewise_rigid_path = os.path.join(transforms_dir, "piecewise_rigid_%02d.mat")
    # XXX: multi_chunk_greedy uses run not execute
    if not os.path.exists(piecewise_rigid_path):
        multi_chunk_greedy.run('-d 2 '
                               '-a '
                               '-dof 6 '
                               '-i {} {} '
                               '-cm {} '
                               '-ia {} '
                               '-m WNCC 4x4 '
                               '-wncc-mask-dilate '
                               '-n 600x600x200x0 '
                               '-search 10000 10 5 '
                               '-wreg 0.05 '
                               '-o {} '.format(reference_scalar_path, moving_scalar_path,
                                               reference_chunk_mask_path, global_rigid_path,
                                               piecewise_rigid_path))


    # Step 5 - Piecewise deformable registration
    piecewise_deformable_path = os.path.join(transforms_dir, "piecewise_deformable_%02d.nii.gz")
    if not os.path.exists(piecewise_deformable_path):
        multi_chunk_greedy.run('-d 2 '
                               '-i {} {} '
                               '-cm {} '
                               '-it {} '
                               '-m WNCC 4x4 '
                               '-wncc-mask-dilate '
                               '-n 400x200x100x20 '
                               '-sv '
                               '-s 0.6mm 0.1mm '
                               '-e 0.25 '
                               '-o {} '.format(reference_scalar_path, moving_scalar_path,
                                               reference_chunk_mask_path,
                                               piecewise_rigid_path,
                                               piecewise_deformable_path))


    # Step 6 - Apply the deformation to the moving slide
    registration_result_path = os.path.join(working_dir, "registered_moving_slide.nii.gz")
    if not os.path.exists(registration_result_path):
        multi_chunk_greedy.run('-d 2 '
                               '-rf {} '
                               '-cm {} '
                               '-r {} {} '
                               '-rb 255 '
                               '-rm {} {} '.format(reference_slide_path,
                                                   reference_chunk_mask_path,
                                                   piecewise_deformable_path, piecewise_rigid_path,
                                                   moving_slide_path, registration_result_path))


    # Save the results to an ITK-SNAP workspace for easy visualization
    workspace_path = os.path.join(working_dir, "registration_result.itksnap")
    if not os.path.exists(workspace_path):
        subprocess.run('itksnap-wt '
                       '-layers-add-anat {} '
                       '-psn "Registered Tau" '
                       '-layers-add-anat {} '
                       '-psn "Nissl" '
                       '-o {} '.format(registration_result_path,
                                       reference_slide_path,
                                       workspace_path),
                       shell=True,
                       stdout=subprocess.DEVNULL,
                    )
