import argparse
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import json
import matplotlib.pyplot as plt
import numpy as np
import SimpleITK as sitk

from histology_data import HistologyData


# https://github.com/pyushkevich/histoannot.git
# git clone the GitHub repository and add the path to the sys.path
# This makes sure you are using the latest version of the code
sys.path.append('/Users/cathalye/Packages/histoannot/')
import phas.client.api as phas
from phas.dltrain import spatial_transform_roi


class Remap:

    def __init__(self, fixed: HistologyData, moving: HistologyData, task: phas.SamplingROITask):
        self.fixed = fixed
        self.moving = moving
        self.task = task

    
