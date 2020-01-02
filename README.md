genmesh

Generate meshes from packing files (xyzd). Depends on OpenCASCADE (7.3) and GMSH (v4+). 

# Installation

1. Install OpenCASCADE (v 7.30 tested)
2. Install GMSH with OpenCASCADE linked (v 4.4.1)
3. Compile genmesh with GMSH linked.

Check out `install.sh` for a template on what to do.

Note: Ensure that \$LD_LIBRARY_PATH points to the OpenCASCADE libs.
Note: Ensure that you use the same compiler version for all dependencies.

## Snippets
```
freetype2, mesa/opengl dev pacakges, libXmu, libXi
./configure --enable-gcc --enable-shared --enable-threads --prefix=/home/IBT/rao/tools/tcl --enable-64bit
./configure --enable-gcc --enable-shared --enable-threads --prefix=/home/IBT/rao/tools/tk --enable-64bit --with-tcl=/home/IBT/rao/tools/tcl/lib
cmake -DCMAKE_INSTALL_PREFIX=../../occt -D3RDPARTY_TCL_LIBRARY_DIR=../../tcl/lib/ -D3RDPARTY_TK_LIBRARY_DIR=../../tk/lib/ -D3RDPARTY_TK_INCLUDE_DIR=../../tk/include/ ..
cmake -DCMAKE_PREFIX_PATH=<path to occt install>
```


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

Using the `create.sh` script also redirects output to a log file in the output directory.

I use preScalingFactor to convert meshes to a size such that bead size = 1, consistent with the mono-full packing.

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
- [x] Save options (so that it can be reloaded and used with exported geometry?): Using logs now.
- [-] Use OCC-fix for modified beads!!: Caused crashes
- [-] Try volume generating: Unnecessary
- [x] write meshes to individual folders. Including logs.
- [x] set rCylDelta after transform-scaling
- [x] FIX: minimum bead radius is zero
- [x] Implement Higher order mesh generation
- [ ] Output surfs and volumes to a separate folder?
- [ ] Bridges at the cylinder-bead interface
- [ ] Check for memleaks. No errors. But some blocks are reachable (beads vector not deleted).
- [ ] Write test inputs to run and verify code.
- [ ] {!}Investigate capped meshes. Why did 400 beads take 5-10 hours? 
- [ ] Manual node placement at contact points 
- [x] Try embedded bead CP and mesh size control: not supported for HXT
- [ ] Scrap the need for ./create.sh. create <file>.log automatically
- [ ] Enlarged beads don't work at 0.001: Fix the rcyl assertion.
- [ ] meshSizeMethod=0 doesn't work
- [x] Extract mesh volume data into variables to output the mesh-scaled volume in stdout.
- [x] Allow changing bead mesh scaling reference value: max/avg or number
- [ ] Porosity should use bed length from bead centers not ends.
- [ ] Mesh Sensitivity
    - [ ] Mesh.RandomFactor(3D)
    - [ ] Mesh.ToleranceInitialDelaunay
- [ ] clean duplicate local variables in packedbed.transform()
- [x] Allow changing fragment output formats
- [x] Makefile depends on relative folder structure, fix it.
- [ ] boost endian conversion is not available easily on JURECA: remove dependency

Known Issues

- Netgen optimizer crashes sometimes. (After mesh size constraints were applied to surfaces)
- Geometry.ScalingFactor doesn't work. Use preScalingFactor instead.
- In 7k-pre and 6k-pre cases, only a handful of beads were actually captured and meshed. [prescaling issue?]
- poly5 fails (beads peek out of container). Just use poly-full.
- A few features do not work in conjunction with HXT mesh algorithm. This is an upstream issue with GMSH.
- Netgen optimizer might not work well with mesh field gradients: It might make the meshes too uniform somehow
