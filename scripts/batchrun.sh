#!/bin/bash

ROOT_DIR=/Users/cathalye/Projects/202506_amygdala_sampling/data

find "$ROOT_DIR" -type f -name 'block_info.json' | while read -r file; do
    python runme.py --task_id 41 --root_dir "$ROOT_DIR" --block_json "$file"
done
