# genmesh

Generate meshes from packing files (xyzd). Depends on OpenCASCADE and GMSH (v 4.2.2). 

## Installation

1. Install OpenCASCADE (v 7.30 tested)
2. Install GMSH with OpenCASCADE link (v 4.2.2)
3. Compile genmesh with GMSH link.

Note: Ensure that \$LD_LIBRARY_PATH points to the OpenCASCADE libs.

## Workflow

A rough concept of execution is as follows:

- Read input file (default.in) for parameters.
- Read packing.xyzd file within limits specified and store bead data.
- Generate geometries (Cylinder, Beads, and Bridges).
- Perform specified boolean operations: fuse/cut + fragment.
- Generate Named Physical Groups for different features.
- Generate mesh.
- Write full mesh to specified output file. 
- Write individual domains (interstitial and beads) to vtk files. 

Here, the term 'bridges' is used to denote cylindrical geometries between individual neighbouring beads. Bridges may be used to 'cap' or 'bridge' the beads they connect. 

## Usage

``` 
./genmesh <optional output filename with extension> 
```

This should create the mesh in the required format in the `outpath` directory. Additionally, two vtk files of the two domains are generated to allow easier examination of the mesh. 

## Todo

1. Use GMSH::logger.
2. Improved error handling.
3. Easy switch between cap/bridge/reduce/enlarge
4. Easy handling of db_dp and offset calculation
5. Easy handling of nbeads
6. Switch off outputs of mesh fragments
7. Mesh.MaxNumThreads3D, Mesh.CharacteristicLengthFactor, Mesh.OptimizeThreshold, 
8. Implement geometry creation in Built in kernel 
    1. Create routines for spheres and cylinders
    2. Modularize and refactor existing code

## Known Issues
1. After mesh size constraints were applied to only surfaces (and points were used inside beads), the Netgen optimizer crashes randomly.
2. Geometry.ScalingFactor doesn't work. Use dilateFactor instead.
3. Large beds either are stuck or crash with errors.

