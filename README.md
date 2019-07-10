# genmesh

Generate meshes from packing files (xyzd). Depends on OpenCASCADE and GMSH (v4+). 

## Installation

1. Install OpenCASCADE (v 7.30 tested)
2. Install GMSH with OpenCASCADE link (v 4.2.2)
3. Compile genmesh with GMSH link.

Note: Ensure that \$LD_LIBRARY_PATH points to the OpenCASCADE libs.

## Workflow

A rough concept of execution is as follows:

- Read input file (default.in) for parameters.
- Read packing.xyzd file and store bead data.
- Sort beads by z co-ordinate and save only nBeads beads.
- Generate geometries (Cylinder, Beads, and Bridges) with applied preScalingFactor. 
- Perform specified boolean operations: fuse/cut + fragment.
- Generate Named Physical Groups for different features.
- Generate mesh.
- Write full mesh to specified output file. 
- Write individual domains (interstitial and beads) to vtk files. 

Here, the term 'bridges' is used to denote cylindrical objects between individual neighbouring beads. Bridges may be used to 'cap' or 'bridge' the beads they connect. 

## Usage

``` 
./genmesh <input file> <optional output filename with extension> 
```

This should create the mesh in the required format in the `outpath` directory. Additionally, two vtk files of the two domains are generated to allow easier examination of the mesh. 

## Todo

- [x] Number of bridges, minimum bead size
- [x] Scale polydisperse bead meshes to smaller beads. -> set lc for individual beads.
- [x] bounding box for beads
- [x] Generate polydisperse beads with bridges.
- [x] Implement Translation
- [ ] Implement PreTranslation
- [ ] Improved error handling.
- [ ] Switch off outputs of mesh fragments
- [ ] Bridges at the cylinder-bead interface
- [ ] Check for memleaks. No errors. But some blocks are reachable.
- [ ] Implement enlarged beads
- [ ] Write test inputs to run and verify code.
- [ ] Implement memory diagnostics. How much memory/cpu-time does the code use?
- [ ] Investigate capped meshes. Why did 400 beads take 5-10 hours? 
- [ ] Cleaner outputs/logging.
- [ ] Check timings for fields vs boundary mesh sizes
- [ ] Better makefile
- [ ] Better argument handling?
- [ ] Switch for fields vs brep meshing
- [ ] Manual node placement at contact points 
- [ ] expert mode? 

Known Issues
- After mesh size constraints were applied to only surfaces (and points were used inside beads), the Netgen optimizer crashes randomly.
    - Does this still happen with fields? 
- Geometry.ScalingFactor doesn't work. Use preScalingFactor instead.
- In 7k-pre and 6k-pre cases, only a handful of beads were actually captured and meshed. [prescaling issue?]
