import os
import subprocess
import sys

import numpy as np
import SimpleITK as sitk

from picsl_greedy import Greedy2D
from picsl_greedy import MultiChunkGreedy2D

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from src.histology_data import HistologyData


greedy = Greedy2D()
multi_chunk_greedy = MultiChunkGreedy2D()

DATA_DIR = "../data/HNL-33-18/HR2a"
FIXED_SLIDE = "../data/HNL-33-18/HR2a/nissl_slide_thumbnail.nii.gz"
MOVING_SLIDE = "../data/HNL-33-18/HR2a/tau_slide_thumbnail.nii.gz"

WORKING_DIR = os.path.join(DATA_DIR, "work")
os.makedirs(WORKING_DIR, exist_ok=True)

reference_slide = HistologyData(task=None, slide_id=None, thumbnail_path=FIXED_SLIDE)
print("read reference slide")
moving_slide = HistologyData(task=None, slide_id=None, thumbnail_path=MOVING_SLIDE)
print("read moving slide")


# Step 1 - Get scalar images of the reference and moving slides
reference_scalar_img = reference_slide.get_single_channel_image(channel=1)
print("got single channel image")
reference_scalar_path = os.path.join(WORKING_DIR, "reference_scalar.nii.gz")
sitk.WriteImage(reference_scalar_img, reference_scalar_path)

moving_scalar_img = moving_slide.get_single_channel_image(channel=1)
print("got single channel image")
moving_scalar_path = os.path.join(WORKING_DIR, "moving_scalar.nii.gz")
sitk.WriteImage(moving_scalar_img, moving_scalar_path)


# Step 2 - Get reference binary and chunk masks
reference_binary_mask = reference_slide.get_binary_mask()
print("got binary mask")
reference_binary_mask_path = os.path.join(WORKING_DIR, "reference_binary_mask.nii.gz")
sitk.WriteImage(reference_binary_mask, reference_binary_mask_path)
# NOTE:
# image_graph_cut is run as a subprocess, so it needs the file paths
# change this once python bindings are available to input and output the masks
# direcrtly instead of passing the file paths
# reference_chunk_mask = reference_slide.get_chunk_mask(reference_binary_mask)
reference_chunk_mask_path = os.path.join(WORKING_DIR, "reference_chunk_mask.nii.gz")
reference_slide.get_chunk_mask(reference_binary_mask_path, reference_chunk_mask_path)
print("got chunk mask")
reference_chunk_mask = sitk.ReadImage(reference_chunk_mask_path)
print("read chunk mask")

print("sizes")
print(reference_scalar_img.GetSize())
print(moving_scalar_img.GetSize())
print(reference_binary_mask.GetSize())
print(reference_chunk_mask.GetSize())

# reference_binary_mask_path = "/Users/cathalye/Downloads/histo_histo_registration-main/data/HNL-33-18/HR2a/registration/reference_binary_mask.nii.gz"

# Step 3 - Global rigid registration
output_global_rigid_path = os.path.join(WORKING_DIR, "output_global_rigid.mat")
# if not os.path.exists(output_global_rigid_path):
#     greedy.execute('-d 2 -a -dof 6 '
#                    '-i reference moving '
#                    '-gm binary_mask '
#                    '-ia-image-centers '
#                    '-m WNCC 4x4 '
#                    '-wncc-mask-dilate '
#                    '-n 200x200x40x0x0 '
#                    '-search 20000 any 10 '
#                    '-o global_rigid ',
#                    reference=reference_scalar_img, moving=moving_scalar_img,
#                    binary_mask=reference_binary_mask,
#                    global_rigid=None)
#     global_rigid = greedy['global_rigid']
#     print(global_rigid)
#     print("got global rigid")
#     np.savetxt(output_global_rigid_path, global_rigid)


if not os.path.exists(output_global_rigid_path):
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
                                   output_global_rigid_path))


# Step 4 - Piecewise rigid registration
output_piecewise_rigid_path = os.path.join(WORKING_DIR, "output_piecewise_rigid_%02d.mat")
# XXX: multi_chunk_greedy uses run not execute
if not os.path.exists(output_piecewise_rigid_path):
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
                                           reference_chunk_mask_path, output_global_rigid_path,
                                           output_piecewise_rigid_path))


# Step 5 - Piecewise deformable registration
output_piecewise_deformable_path = os.path.join(WORKING_DIR, "output_piecewise_deformable_%02d.nii.gz")
if not os.path.exists(output_piecewise_deformable_path):
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
                                           output_piecewise_rigid_path,
                                           output_piecewise_deformable_path))


# Step 6 - Apply the deformation to the moving slide
registration_result_path = os.path.join(WORKING_DIR, "registration_result.nii.gz")
if not os.path.exists(registration_result_path):
    multi_chunk_greedy.run('-d 2 '
                           '-rf {} '
                           '-cm {} '
                           '-r {} {} '
                           '-rb 255 '
                           '-rm {} {} '.format(FIXED_SLIDE,
                                               reference_chunk_mask_path,
                                               output_piecewise_deformable_path, output_piecewise_rigid_path,
                                               MOVING_SLIDE, registration_result_path))


# Save the results to an ITK-SNAP workspace for easy visualization
workspace_path = os.path.join(WORKING_DIR, "registration_results.itksnap")
if not os.path.exists(workspace_path):
    subprocess.run(
        'itksnap-wt -layers-add-anat {} -psn "Registered Tau" -layers-add-anat {} -psn "Nissl" -o {}'.format(
            registration_result_path, FIXED_SLIDE, workspace_path
        ),
        shell=True,
        stdout=subprocess.DEVNULL,
    )
