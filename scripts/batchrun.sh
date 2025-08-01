#!/bin/bash

python remap.py --phas_url https://chead.uphs.upenn.edu \
                --private_key /Users/cathalye/.private/chead_api_key.json \
                --task_id 26 \
                --fixed_slide_id 33568 \
                --moving_slide_id 33632 \
                --moving_slide_thumbnail_path ../data/at8.nii.gz \
                --registration_dir ../data/work
