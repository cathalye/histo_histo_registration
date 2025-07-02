import os
import sys
import json

import numpy as np

sys.path.append('/Users/cathalye/Packages/histoannot/')
import phas.client.api as phas
from phas.dltrain import spatial_transform_roi


def connect_to_server(task_id):
    # PHAS connection parameters
    PHAS_URL='https://histo.itksnap.org'
    PRIVATE_KEY = '/Users/cathalye/.private/histoitk_api_key.json'

    conn = phas.Client(PHAS_URL, PRIVATE_KEY)
    # Create a sampling ROI task object to pass to Slide class for downloading sampling ROI json
    task = phas.SamplingROITask(conn, task_id)

    return task


def reset_slide(task: phas.SamplingROITask, slide_id):
    task.delete_sampling_rois_on_slide(slide_id)


def read_json_property(jsonfile, property):
    with open(jsonfile, 'r') as file:
        data = json.load(file)
        slide_id = data.get(f'{property}', None)
    return slide_id


def process_roi_data(roi_data, type, scale=1):
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
