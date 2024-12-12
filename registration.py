import os
import subprocess
import sys

import SimpleITK as sitk

# from picsl_c3d import Convert3D
# from picsl_c3d import Convert2D
# from picsl_greedy import Greedy3D
# from picsl_greedy import Greedy2D

ROOT_DIR = "./data"

for specimen in os.listdir(ROOT_DIR):
    specimen_dir = os.path.join(ROOT_DIR, specimen)
    for block in os.listdir(specimen_dir):
        block_dir = os.path.join(specimen_dir, block)

        print(f"Processing {specimen} {block}")

        reference_thumbnail = os.path.join(block_dir, "nissl_slide_thumbnail.nii.gz")
        moving_thumbnail = os.path.join(block_dir, "tau_slide_thumbnail.nii.gz")

        # Only proceed if both the Nissl and Tau thumbnails are available
        if not os.path.exists(reference_thumbnail) or not os.path.exists(
            moving_thumbnail
        ):
            print(f"Skipping {specimen} {block}")
            continue

        # read the nissl and tau thumbnails
        reference_img = sitk.ReadImage(reference_thumbnail)
        moving_img = sitk.ReadImage(moving_thumbnail)

        # Step 1 - Extract green channel to get the scalar images
        reference_scalar_img = sitk.VectorIndexSelectionCast(reference_img, 1)
        moving_scalar_img = sitk.VectorIndexSelectionCast(moving_img, 1)

        registration_dir = os.path.join(block_dir, "registration")
        os.makedirs(registration_dir, exist_ok=True)

        reference_scalar = f"{registration_dir}/reference_scalar.nii.gz"
        moving_scalar = f"{registration_dir}/moving_scalar.nii.gz"

        sitk.WriteImage(reference_scalar_img, reference_scalar)
        sitk.WriteImage(moving_scalar_img, moving_scalar)

        # Step 2 - Get binanry mask of the reference image
        # Otsu thresholding to get a binary mask
        reference_binary_mask = f"{registration_dir}/reference_binary_mask.nii.gz"
        # If manual mask exists, don't overwrite with automatic mask
        if not os.path.exists(reference_binary_mask):
            otsu_filter = sitk.OtsuThresholdImageFilter()
            otsu_filter.SetInsideValue(1)
            otsu_filter.SetOutsideValue(0)
            binary_image = otsu_filter.Execute(reference_scalar_img)

            # Otsu thresholding segmentation has holes in the mask due to the irregular
            # staining. Closing to fill in holes.
            structuring_element = sitk.sitkBall
            radius = [5, 5, 5]
            closed_image = sitk.BinaryMorphologicalClosing(
                binary_image, radius, structuring_element
            )
            sitk.WriteImage(closed_image, reference_binary_mask)

        # Step 3 - Divide the mask into chunks using image graph cut
        # For simplicity, we will define a fixed number of chunks
        N_CHUNKS = 10
        reference_multi_chunk = f"{registration_dir}/reference_multi_chunk.nii.gz"
        if not os.path.exists(reference_multi_chunk):
            subprocess.run(
                "image_graph_cut -u 1.2 -n 100 -c 4 0.1 {}/reference_binary_mask.nii.gz {}/reference_multi_chunk.nii.gz {}".format(
                    registration_dir, registration_dir, N_CHUNKS
                ),
                shell=True,
                stdout=subprocess.DEVNULL,
            )

        # Step 4 - Global rigid registration
        output_global_rigid = f"{registration_dir}/output_global_rigid.mat"
        if not os.path.exists(output_global_rigid):
            subprocess.run(
                "greedy -d 2 -a -dof 6 -i {} {} -gm {} -ia-image-centers -m WNCC 4x4 -wncc-mask-dilate -n 200x200x40x0x0 -search 20000 any 10 -o {}".format(
                    reference_scalar,
                    moving_scalar,
                    reference_binary_mask,
                    output_global_rigid,
                ),
                shell=True,
                stdout=subprocess.DEVNULL,
            )

        # Step 5 - Piecewise rigid registration
        # We will use the chunks from the image graph cut as the regions
        output_piecewise_rigid = f"{registration_dir}/output_piecewise_rigid_%02d.mat"
        if not os.path.exists(output_piecewise_rigid):
            subprocess.run(
                "multi_chunk_greedy -d 2 -a -dof 6 -i {} {} -cm {} -ia {} -m WNCC 4x4 -wncc-mask-dilate -n 600x600x200x0 -search 10000 10 5 -wreg 0.05 -o {}".format(
                    reference_scalar,
                    moving_scalar,
                    reference_multi_chunk,
                    output_global_rigid,
                    output_piecewise_rigid,
                ),
                shell=True,
                stdout=subprocess.DEVNULL,
            )

        # # Step 6 - Piecewise deformable registration
        output_piecewise_deformable = (
            f"{registration_dir}/output_piecewise_deformable_%02d.nii.gz"
        )
        if not os.path.exists(output_piecewise_deformable):
            subprocess.run(
                "multi_chunk_greedy -d 2 -i {} {} -cm {} -it {} -m WNCC 4x4 -wncc-mask-dilate -n 400x200x100x20 -sv -s 0.6mm 0.1mm -e 0.25 -o {}".format(
                    reference_scalar,
                    moving_scalar,
                    reference_multi_chunk,
                    output_piecewise_rigid,
                    output_piecewise_deformable,
                ),
                shell=True,
                stdout=subprocess.DEVNULL,
            )

        # Step 7 - Apply the deformation to the tau image
        registration_output = f"{block_dir}/tau_to_nissl.nii.gz"
        if not os.path.exists(registration_output):
            subprocess.run(
                "multi_chunk_greedy -d 2 -rf {} -cm {} -r {} {} -rb 255 -rm {} {}".format(
                    reference_thumbnail,
                    reference_multi_chunk,
                    output_piecewise_deformable,
                    output_piecewise_rigid,
                    moving_thumbnail,
                    registration_output,
                ),
                shell=True,
                stdout=subprocess.DEVNULL,
            )

        result_workspace = f"{block_dir}/result_workspace.itksnap"
        if not os.path.exists(result_workspace):
            subprocess.run(
                'itksnap-wt -layers-add-anat {} -psn "Registered Tau" -layers-add-anat {} -psn "Nissl" -o {}'.format(
                    registration_output, reference_thumbnail, result_workspace
                ),
                shell=True,
                stdout=subprocess.DEVNULL,
            )
