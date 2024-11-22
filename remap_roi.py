import subprocess
import sys

import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
import SimpleITK as sitk


# https://github.com/pyushkevich/histoannot.git
# git clone the GitHub repository and add the path to the sys.path
# This makes sure you are using the latest version of the code
sys.path.append('/Users/cathalye/Packages/histoannot/')
import phas.client.api as phas
from phas.dltrain import spatial_transform_roi, draw_sampling_roi, compute_sampling_roi_bounding_box


# To get the API key you need to go to https://histo.itksnap.org/auth/api/generate_key
PHAS_URL='https://histo.itksnap.org'
PRIVATE_KEY = '/Users/cathalye/.private/histoitk_api_key.json'
conn = phas.Client(PHAS_URL, PRIVATE_KEY)
TASK_ID = 39
# Create a sampling ROI task object to pass to Slide class for downloading sampling ROI json
task = phas.SamplingROITask(conn, TASK_ID)

def read_json_property(jsonfile, property):
    with open(jsonfile, 'r') as file:
        data = json.load(file)
        slide_id = data.get(f'{property}', None)
    return slide_id

def find_scaling_factor(task, slide_id):
    x, y = phas.Slide(task, slide_id).dimensions
    # When downsampling, the largest dimension is scaled to 1000
    if x > y:
        scale = 1000 / x
    else:
        scale = 1000 / y
    return scale

# Single test subject
ROOT_DIR = "/Users/cathalye/Library/CloudStorage/Box-Box/chinmayee/202411_pmam_late/data/"
# Define this for dummy data
block_json = "/Users/cathalye/Library/CloudStorage/Box-Box/chinmayee/202411_pmam_late/data/INDD102373/HL3a/block_info.json"
# Read from json
specimen = read_json_property(block_json, 'specimen')
block = read_json_property(block_json, 'block')
nissl_slide_id = read_json_property(block_json, 'nissl_slide_id')
tau_slide_id = read_json_property(block_json, 'tau_slide_id')
# Registration files
registration_dir = f"{ROOT_DIR}/{specimen}/{block}/registration"
global_rigid = np.loadtxt(f"{registration_dir}/output_global_rigid.mat")
multi_chunk_mask = sitk.ReadImage(f"{registration_dir}/reference_multi_chunk.nii.gz")
multi_chunk_mask = np.squeeze(sitk.GetArrayFromImage(multi_chunk_mask))

# Read the downsampled nissl slide and get the spacing
# nissl_slide_downsampled = sitk.ReadImage(f"{ROOT_DIR}/{specimen}/{block}/nissl_slide_thumbnail.nii.gz")
nissl_slide = phas.Slide(task, nissl_slide_id)
nissle_slide_spacing = nissl_slide.spacing
nissl_thumbnail = sitk.ReadImage(f"{ROOT_DIR}/{specimen}/{block}/nissl_slide_thumbnail.nii.gz")
nissl_thumbnail_spacing = nissl_thumbnail.GetSpacing()

chunk = '01'
chunk_rigid = np.loadtxt(f"{registration_dir}/output_piecewise_rigid_{chunk}.mat")
chunk_warp = sitk.ReadImage(f"{registration_dir}/output_piecewise_deformable_{chunk}.nii.gz")

# Get the sampling ROI coordinates
# XXX: roi_task.slide_sampling_rois() returns a list of ROIs, but we are only using the first one
roi = task.slide_sampling_rois(nissl_slide_id)[1]
roi_data = json.loads(roi['json'])

# If we are scaling the ROI
x_roi, y_roi = zip(*roi_data['data'])
assert len(x_roi) == len(y_roi), "x_roi and y_roi must have the same length"
# Convert to numpy array to do element-wise operations
x_roi = np.array(x_roi)
y_roi = np.array(y_roi)

scale = find_scaling_factor(task, nissl_slide_id)
x_scaled_roi = ((x_roi+0.5) * scale).astype(int)
y_scaled_roi = ((y_roi+0.5) * scale).astype(int)

print(x_scaled_roi)

chunk_label = []
for i in range(len(x_scaled_roi)):
    x = x_scaled_roi[i]
    y = y_scaled_roi[i]
    chunk_label.append(multi_chunk_mask[y, x])

chunk_label = np.array(chunk_label)
print(f"The ROI lies in chunk # {np.unique(chunk_label)}")

# xy_scaled coordinates are correct - double checked with the ROI
# segmentation masks in ITK-SNAP
xy_scaled = [list(pair) for pair in zip(x_scaled_roi, y_scaled_roi)]
roi_data_scaled = roi_data.copy()
roi_data_scaled['data'] = xy_scaled
print(roi_data_scaled)

# displacement = chunk_warp.EvaluateAtPhysicalPoint((-20, -40))
# print("Displacement at (-20, -40):", displacement)

def my_transform(roi_xy):
    xy = np.array(roi_xy)

    # Get the coordinates in physical space
    # xy_phys = (xy + 0.5) * nissl_slide.spacing *(1) # if using full resolution xy
    # xy_phys = (xy + 0.5) * nissl_thumbnail_spacing[0] * (-1) # if using downsampled xy
    xy_phys = np.array(chunk_warp.TransformIndexToPhysicalPoint([int(xy[0]), int(xy[1])]))

    displacement = chunk_warp.EvaluateAtPhysicalPoint(xy_phys)
    xy_warp = xy_phys + displacement

    xy_warp = -xy_warp # transformation to RAS??

    xy_warp_homogeneous = np.append(xy_warp, 1) # [x,y,1]
    xy_warp_homogeneous = xy_warp_homogeneous.reshape(3,1)

    # Apply the chunk rigid transform
    xy_chunk_rigid = np.dot(chunk_rigid, xy_warp_homogeneous) # [x',y',1]

    # Apply the global rigid transform
    xy_global = np.dot(global_rigid, xy_chunk_rigid) # [x'',y'',1]

    xy_global = -xy_global # transformation to RAS??

    # Convert from physical space to image space
    # xy_warp_img = xy_global / nissl_slide.spacing * (-1) - 0.5 # if using full resolution xy
    # xy_warp_img = xy_warp / nissl_thumbnail_spacing[0] * (-1) - 0.5 # if using downsampled xy

    xy_global = xy_global[:2] # [x'',y'']

    xy_remap = np.array(chunk_warp.TransformPhysicalPointToIndex([float(xy_global[0]), float(xy_global[1])]))

    return xy_remap[0], xy_remap[1]


# roi_warped = spatial_transform_roi(roi_data, my_transform) # if using full resolution xy
roi_warped = spatial_transform_roi(roi_data_scaled, my_transform) # if using downsampled xy
print(roi_warped)


# Scale the transformed points back to full resolution
x_warped, y_warped = zip(*roi_warped['data'])
assert len(x_warped) == len(y_warped), "x_roi and y_roi must have the same length"
# Convert to numpy array to do element-wise operations
x_warped = np.array(x_warped)
y_warped = np.array(y_warped)


scale_ihc = find_scaling_factor(task, tau_slide_id)
x_full_res = x_warped / scale_ihc
y_full_res = y_warped / scale_ihc

xy_full_res = [list(pair) for pair in zip(x_full_res, y_full_res)]
roi_warped['data'] = xy_full_res

print("origina", roi_data)
print("warped", roi_warped)

roi_id = task.create_sampling_roi(tau_slide_id, roi['label'], roi_warped)
# print(roi_id)
