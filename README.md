# histo_histo_registration
Code for IHC to reference histology slide registration

TODOs
- move registration code to python_c3d, modules, and executable file
- clean up imports
- document the nearest_chunk_map notebook
- redundancies in scaling, dingle channel images, etc.
- script for downloading files instead of notebook
  - handle cases where there are multiple annotated nissl slides but only one of them has a corresponding tau slide. Right now, the code is overwriting nissl tau slide numbers
