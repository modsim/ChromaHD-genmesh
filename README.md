genmesh

Generate meshes from packing files (xyzd). Depends on OpenCASCADE (7.3) and GMSH (v4+). 

# Installation

1. Install OpenCASCADE (v 7.30 tested)
3. Install GMSH with OpenCASCADE linked (v 4.4.1)
4. Compile genmesh with GMSH linked.

Note: Ensure that \$LD_LIBRARY_PATH points to the OpenCASCADE libs.

## Snippets
`freetype2, mesa/opengl dev pacakges, libXmu, libXi`
`./configure --enable-gcc --enable-shared --enable-threads --prefix=/home/IBT/rao/tools/tcl --enable-64bit`
`./configure --enable-gcc --enable-shared --enable-threads --prefix=/home/IBT/rao/tools/tk --enable-64bit --with-tcl=/home/IBT/rao/tools/tcl/lib`                                                                                
`cmake -DCMAKE_INSTALL_PREFIX=../../occt -D3RDPARTY_TCL_LIBRARY_DIR=../../tcl/lib/ -D3RDPARTY_TK_LIBRARY_DIR=../../tk/lib/ -D3RDPARTY_TK_INCLUDE_DIR=../../tk/include/ ..`                                                       

# Workflow

A rough concept of execution is as follows:

- Read input file (default.in) for parameters.
- Read packing.xyzd file and store bead data.
- Sort beads by z co-ordinate and save only nBeads beads.
- transform beads
- Generate geometries (Cylinder, Beads, and Bridges) 
- Perform specified boolean operations: fuse/cut + fragment.
- Generate Named Physical Groups for different features.
- save/load geometry
- Generate mesh.
- Write full mesh to specified output file. 
- Write individual domains (interstitial and beads) to vtk files. 

Here, the term 'bridges' is used to denote conical/cylindrical objects between individual neighbouring beads. Bridges may be used to 'cap' or 'bridge' the beads they connect. 

# Usage

``` 
./genmesh <input file> <optional output filename with extension> 
```

This should create the mesh in the required format in the `outpath` directory. Additionally, two vtk files of the two domains are generated to allow easier examination of the mesh. 

I use preScalingFactor to convert meshes to a size such that bead size = 1, consistent with the mono-full

# Todo

- [x] Number of bridges, minimum bead size
- [x] Scale polydisperse bead meshes to smaller beads. -> set lc for individual beads.
- [x] bounding box for beads
- [x] Generate polydisperse beads with bridges.
- [x] Implement Translation
- [x] Check timings for fields vs boundary mesh sizes
- [x] Implement PreTranslation
- [x] Better preScalingFactor implementation
- [x] Implement auto cylinder boundaries
- [x] Cleaner outputs/logging.
- [x] Improve default.in: folds should allow quick selection of kugelpackung
- [x] Export geometry before mesh start: [ Possible ].
- [x] Time spent on booleans
- [x] Switch for fields vs brep meshing
- [x] Fix poly capping on larger beads: use cones and frustums
- [x] output git commit + state into stdout
- [x] Switch off outputs of mesh fragments 
- [x] Better makefile
- [x] Check if it's possible to import meshes and modify them.
- [x] Figure out how to calculate mesh volume for different components
- [x] git state isn't perfect since it is dirty if default.in is changed before compile
- [x] Allow cylinder constrained creation again(zbot, ztop)
- [x] RCYL is calclated after bead radius modification: Change it to depend on initial packing
- [x] Implement enlarged beads
- [-] Save options (so that it can be reloaded and used with exported geometry?)
- [-] Use OCC-fix for modified beads!!: Caused crashes
- [-] Try volume generating: Unnecessary
- [ ] Improved error handling.
- [ ] Output surfs and volumes to a separate folder?
- [ ] Bridges at the cylinder-bead interface
- [ ] Check for memleaks. No errors. But some blocks are reachable (beads vector not deleted).
- [ ] Write test inputs to run and verify code.
- [ ] {!}Investigate capped meshes. Why did 400 beads take 5-10 hours? 
- [ ] Manual node placement at contact points 
- [ ] Try embedded bead CP and mesh size control
- [ ] Scrap the need for ./create.sh. create <file>.log automatically
- [ ] Enlarged beads don't work at 0.001: Fix the rcyl assertion.
- [ ] meshSizeMethod=0 doesn't work
- [x] write meshes to individual folders. Including logs.
- [ ] allow setting rCyl directly?
- [x] set rCylDelta after transform-scaling
- [ ] Extract mesh volume data into variables to output the mesh-scaled volume in stdout.
- [ ] Mesh Sensitivity
    - [ ] Mesh.RandomFactor(3D)
    - [ ] Mesh.ToleranceInitialDelaunay

Known Issues
- Netgen optimizer crashes sometimes. (After mesh size constraints were applied to surfaces)
- Geometry.ScalingFactor doesn't work. Use preScalingFactor instead.
- In 7k-pre and 6k-pre cases, only a handful of beads were actually captured and meshed. [prescaling issue?]
- poly5 fails (beads peek out of container). Just use poly-full.
- meshSizeMethod = 0 doesn't work
