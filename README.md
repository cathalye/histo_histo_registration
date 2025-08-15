# Histology-to-Histology Registration

A comprehensive toolkit for registering and remapping regions of interest (ROIs) between different histology slide stains using advanced image registration techniques.

## Overview

This project provides tools for:
- **Downloading histology slides** from [PHAS (PICSL Histology Annotation Server)](https://github.com/pyushkevich/histoannot)
- **Multi-stage image registration** between reference and moving histology slides
- **ROI remapping** from one stain to another using registration transforms
- **Batch processing** capabilities for large-scale histology analysis

The registration pipeline uses a hierarchical approach combining global rigid, piecewise rigid, and piecewise deformable registration to achieve high-quality alignment between different histology stains.

## Architecture

### Core Components

- **`HistologyData`**: Core class for handling histology slide data, including thumbnail management and coordinate transformations
- **`RemapROI`**: Handles ROI coordinate transformation and remapping between registered slides
- **Registration Pipeline**: Multi-stage registration using Greedy algorithm with chunk-based processing

### Registration Pipeline

1. **Global Rigid Registration**: Initial alignment using 6-degree-of-freedom transformation
2. **Piecewise Rigid Registration**: Local rigid transformations applied to image chunks
3. **Piecewise Deformable Registration**: Local non-linear deformations for fine alignment
4. **Transform Application**: Final registration result generation

## Installation

### Dependencies

```bash
# Create and activate conda environment
conda env create -f environment.yml
conda activate greedy_registration
```

### External Tools

1. **git clone**: All the PICSL packages are open-source. You can clone the github repository and it to path using `sys.path` to use the latest development version:
   - `phas` - https://github.com/pyushkevich/histoannot
   - `picsl_greedy` - https://github.com/pyushkevich/greedy_python
   - `picsl_image_graph_cut` - https://github.com/pyushkevich/image-graph-cut


2. **pip install**: All PICSL packages are also available via pip and automatically installed when you create conda environment above:
   - `phas`: API for PICSL histology annotation server
   - `picsl_greedy`: For image registration
   - `picsl_image_graph_cut`: For chunk mask generation

3. **ITK-SNAP**: For visualization of registration results
   - Download from [ITK-SNAP website](http://www.itksnap.org/)

## Usage

### 1. Download Histology Slides

For chead or histo-itk servers, generate api_key here:
   - `https://chead.uphs.upenn.edu/auth/api/generate_key`
   - `https://histo.itksnap.org/auth/api/generate_key`

Download slides and ROIs from a PHAS server:

```bash
python scripts/download_slides.py \
    --server https://your-server.com \
    --private_key /path/to/private_key.json \
    --task_id 123 \
    --ref_stain "H&E" \
    --mov_stain "IHC" \
    --root_dir ./data
```

**Parameters:**
- `--server`: PHAS server URL
- `--private_key`: Path to server authentication key
- `--task_id`: Task ID on the server
- `--ref_stain`: Reference stain with annotations
- `--mov_stain`: Moving stain to be registered
- `--root_dir`: Output directory for downloaded data

**Output Structure:**
```
root_dir/
├── specimen_1/
│   └── block_1/
│       ├── he_slide_thumbnail.nii.gz
│       ├── ihc_slide_thumbnail.nii.gz
│       ├── he_sampling_roi.nii.gz
│       └── block_info.json
└── ...
```

### 2. Register Histology Slides

Perform multi-stage registration between slides:

```bash
python scripts/registration.py \
    --reference_slide ./data/specimen_1/block_1/he_slide_thumbnail.nii.gz \
    --moving_slide ./data/specimen_1/block_1/ihc_slide_thumbnail.nii.gz \
    --working_dir ./work
```

**Registration Stages:**
1. **Global Rigid**: 6-DOF transformation with WNCC metric
2. **Piecewise Rigid**: Local rigid transformations per chunk
3. **Piecewise Deformable**: Local non-linear deformations
4. **Transform Application**: Final registered image generation

**Output Files:**
```
working_dir/
├── transforms/
│   ├── global_rigid.mat
│   ├── piecewise_rigid_00.mat
│   ├── piecewise_deformable_00.nii.gz
│   └── ...
├── reference_scalar.nii.gz
├── moving_scalar.nii.gz
├── reference_binary_mask.nii.gz
├── reference_chunk_mask.nii.gz
├── registered_moving_slide.nii.gz
└── registration_result.itksnap
```

### 3. Remap ROIs

Transfer ROIs from reference to moving stain using registration transforms:

```bash
python scripts/remap.py \
    --phas_url https://your-server.com \
    --private_key /path/to/private_key.json \
    --task_id 123 \
    --fixed_slide_id 33568 \
    --moving_slide_id 33632 \
    --moving_slide_thumbnail_path ./data/specimen_1/block_1/ihc_slide_thumbnail.nii.gz \
    --registration_dir ./work
```

**Process:**
1. Extract ROI coordinates from reference slide
2. Apply registration transforms to remap coordinates
3. Create new ROIs on moving slide in PHAS

### 4. Batch Processing

Use the provided batch script for automated processing:

```bash
bash scripts/batchrun.sh
```

**Customize the script with your parameters:**
- PHAS server URL and credentials
- Slide IDs for reference and moving stains
- File paths for thumbnails and registration results

## API Reference

### HistologyData Class

Core class for managing histology slide data and operations.

```python
from src.histology_data import HistologyData

# Initialize with PHAS task and slide ID
slide = HistologyData(task, slide_id)

# Or initialize with local thumbnail file
slide = HistologyData(None, None, thumbnail_path="./slide.nii.gz")
```

**Key Methods:**
- `get_single_channel_image(channel, save_path, return_img)`: Extract single channel
- `get_binary_mask(channel, save_path, return_img)`: Generate binary mask using Otsu thresholding
- `get_chunk_mask(binary_mask_path, chunk_mask_path, n_parts)`: Create chunk-based segmentation
- `get_thumbnail_coord_from_full(x, y)`: Convert full-resolution to thumbnail coordinates
- `get_full_coord_from_thumbnail(x, y)`: Convert thumbnail to full-resolution coordinates

### RemapROI Class

Handles ROI coordinate transformation and remapping.

```python
from src.remap_roi import RemapROI

remap = RemapROI(registration_dir, moving_slide)
transformed_coords = remap.registration_transform([x, y])
```

**Key Methods:**
- `registration_transform(xy)`: Apply registration transforms to coordinates
- `get_chunk_transforms(x, y)`: Get transforms for specific image chunk

## File Formats

### Input Formats
- **NIfTI (.nii.gz)**: Histology slide thumbnails and masks
- **JSON**: ROI coordinate data and block metadata
- **MAT**: Rigid transformation matrices

### Output Formats
- **NIfTI**: Registered images, masks, and deformation fields
- **ITK-SNAP Workspace**: Visualization workspace for registration results
- **JSON**: Updated ROI coordinates and metadata

## Coordinate Systems

The system handles coordinate transformations between:
- **Full Resolution**: Original slide dimensions
- **Thumbnail**: Downsampled (max dimension = 1000px)
- **Physical Space**: ITK physical coordinates
- **Index Space**: Pixel/voxel coordinates

**Important Notes:**
- SimpleITK uses LPS coordinate system
- Greedy uses RAS coordinate system
- How the conversion is handled is documents in the `RemapROI.registration_transform` function

## Chunk-Based Processing

The registration uses a chunk-based approach:
1. **Binary Mask Generation**: Otsu thresholding with morphological closing
2. **Chunk Segmentation**: Graph-cut based partitioning (default: 10 chunks)
3. **Local Registration**: Independent transforms per chunk
4. **Nearest Chunk Mapping**: Distance-based chunk assignment for ROI transformation

## Troubleshooting

### Common Issues

1. **PHAS Connection Errors**
   - Verify server URL and private key
   - Check network connectivity
   - Ensure proper PHAS installation

2. **Registration Failures**
   - Verify input image formats
   - Check working directory permissions
   - Ensure Greedy installation

3. **ROI Remapping Errors**
   - Verify registration completion
   - Check transform file existence
   - Validate coordinate system consistency

### Debug Mode

Enable verbose output for debugging:
```bash
# Add --verbose flag to scripts where supported
python scripts/registration.py --verbose --reference_slide ... --moving_slide ...
```

## Citation

If you use this software in your research, please cite:

```
@article{athalye2025,
  title = {Operationalizing Postmortem Pathology-{{MRI}} Association Studies in {{Alzheimer}}’s Disease and Related Disorders with {{MRI-guided}} Histology Sampling},
  author = {Athalye, Chinmayee and Bahena, Alejandra and Khandelwal, Pulkit and Emrani, Sheina and Trotman, Winifred and Levorse, Lisa M. and Khodakarami, Zahra and Ohm, Daniel T. and Teunissen-Bermeo, Eric and Capp, Noah and Sadaghiani, Shokufeh and Arezoumandan, Sanaz and Lim, Sydney A. and Prabhakaran, Karthik and Ittyerah, Ranjit and Robinson, John L. and Schuck, Theresa and Lee, Edward B. and Tisdall, M. Dylan and Das, Sandhitsu R. and Wolk, David A. and Irwin, David J. and Yushkevich, Paul A.},
  date = {2025-05-28},
  journaltitle = {Acta Neuropathologica Communications},
  volume = {13},
  number = {1},
  pages = {120},
  issn = {2051-5960},
  doi = {10.1186/s40478-025-02030-y},
  url = {https://doi.org/10.1186/s40478-025-02030-y},
}

```
