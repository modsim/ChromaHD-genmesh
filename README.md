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

Note: 

- preScalingFactor is used because packings might be generated at different scales.
- the preScalingFactor is used on zTop and zBot *before* selecting the beads. This means that the zTop and zBot given by the user will roughly correspond to the size of the packed bed *AFTER* prescaling the bed. 
- Here, the term 'bridges' is used to denote conical/cylindrical objects between individual neighbouring beads. Bridges may be used to 'cap' or 'bridge' the beads they connect. 

# Usage

``` 
./genmesh <input file> <optional output filename with extension> 
```

This should create the mesh in the required format in the `outpath` directory. Additionally, two vtk files of the two domains are generated to allow easier examination of the mesh. 

Using the `create.sh` script also redirects output to a log file in the output directory.

I use preScalingFactor to convert meshes to a size such that bead size = 1, consistent with the mono-full packing.

# Todo

- [DROP] Use OCC-fix for modified beads!!: Caused crashes
- [DROP] Try volume generating: Unnecessary
- [DONE] Number of bridges, minimum bead size
- [DONE] Scale polydisperse bead meshes to smaller beads. -> set lc for individual beads.
- [DONE] bounding box for beads
- [DONE] Generate polydisperse beads with bridges.
- [DONE] Implement Translation
- [DONE] Check timings for fields vs boundary mesh sizes
- [DONE] Implement PreTranslation
- [DONE] Better preScalingFactor implementation
- [DONE] Implement auto cylinder boundaries
- [DONE] Cleaner outputs/logging.
- [DONE] Improve default.in: folds should allow quick selection of kugelpackung
- [DONE] Export geometry before mesh start: [ Possible ].
- [DONE] Time spent on booleans
- [DONE] Switch for fields vs brep meshing
- [DONE] Fix poly capping on larger beads: use cones and frustums
- [DONE] output git commit + state into stdout
- [DONE] Switch off outputs of mesh fragments 
- [DONE] Better makefile
- [DONE] Check if it's possible to import meshes and modify them.
- [DONE] Figure out how to calculate mesh volume for different components
- [DONE] git state isn't perfect since it is dirty if default.in is changed before compile
- [DONE] Allow cylinder constrained creation again(zbot, ztop)
- [DONE] RCYL is calclated after bead radius modification: Change it to depend on initial packing
- [DONE] Implement enlarged beads
- [DONE] Save options (so that it can be reloaded and used with exported geometry?): Using logs now.
- [DONE] write meshes to individual folders. Including logs.
- [DONE] set rCylDelta after transform-scaling
- [DONE] FIX: minimum bead radius is zero
- [DONE] Implement Higher order mesh generation
- [DONE] Try embedded bead CP and mesh size control: not supported for HXT yet
- [DONE] Extract mesh volume data into variables to output the mesh-scaled volume in stdout.
- [DONE] Allow changing bead mesh scaling reference value: max/avg or number
- [DONE] Allow changing fragment output formats
- [DONE] Makefile depends on relative folder structure, fix it.
- [DONE] print exact version of gmsh used
- [DONE] Checkpoint saves for 1D & 2D Meshes
- [DONE] Implement porosity control: manipulate porosity by adding/removing beads
    - [DONE] Removing beads
    - [TASK] Adding beads
    - [DONE] Improved algorithm: Deletion zones (% based? diameter based?)
- [TASK] Bridges at the cylinder-bead interface
- [TASK] Check for memleaks. No errors. But some blocks are reachable (beads vector not deleted).
- [TASK] Write test inputs to run and verify code.
- [TASK] Scrap the need for ./create.sh. create <file>.log automatically
- [TASK] Enlarged beads don't work at 0.001: Fix the rcyl assertion.
- [TASK] meshSizeMethod=0 doesn't work
- [TASK] Mesh Sensitivity
    - [TASK] Mesh.RandomFactor(3D)
    - [TASK] Mesh.ToleranceInitialDelaunay
- [TASK] clean duplicate local variables in packedbed.transform()
- [TASK] boost endian conversion is not available easily on JURECA: remove dependency
- [TASK] OCC parallel boolean uses up all cores
- [TASK] Plugin(AnalyseCurvedMesh)
- [TASK] Plugin(DiscretizationError) 
- [TASK] Implement debug/release version handling?
- [TASK] it might be neater to scale bed, updatebounds, then translate bed
- [TASK] Make changes necessary for XNS Generic Implementation of Chromatography
- [TASK] Consider creating a mesh.info output with all the mesh data (lengths, volumes, quality etc) in json format
- [TASK] Use rCyl, xyz Cyl if provided
- [TASK] Check that the output format is viable before starting the mesh
- [TASK] Improved centering algorithm for packed bed translation to origin
- [TASK] makefile git-check/version impedes parallel build.
- [TASK] Write better documentation: Manual etc. 
- [TASK] Clean GeomInFile uses
- [TASK] Check if rectangular mesh sims run or not
- [TASK] Intuitive control & code for mesh sizing
- [TASK] Generate periodic meshes
- [TASK] More intuitive input parameters for dimensions etc. 
    - [TASK] zBot/zTop
    - [TASK] inlet/outlet
    - [TASK] nBeads
    - [TASK] Porosity Control
- [TASK] Fix GeomInFile and GeomOutFile in log output when using as input
- [TASK] Better architecture
- [TASK] let genmesh run in a directory with input file and generate necessary files in the same directory. No output subdir nonsense.

[PROJ: Periodic]
    - [TASK] Create flag
    - [TASK] copy geometry 4x
    - [TASK] Cut with plane
    - [TASK] Generate periodic mesh

Known Issues
- Netgen optimizer crashes sometimes. (After mesh size constraints were applied to surfaces)
- Geometry.ScalingFactor doesn't work. Use preScalingFactor instead.
- poly5 fails (beads peek out of container). Just use poly-full.
    - This problem persists with poly-full for certain sections. Can be fixed by a better centering method instead of only x/y based centering
    - It can be sidestepped by running genmesh with coarse element size, and checking if any beads are meshed AFTER the cylinder/planes & have a higher entity number.
- A few features do not work in conjunction with HXT mesh algorithm. This is an upstream issue with GMSH.
- Netgen optimizer might not work well with mesh field gradients: It might make the meshes too uniform somehow

# Keywords in config.in
- packing: path to packing file xyzd in single precision binary little endian format
- por_target: target value of real porosity to attain. By default, beads are removed from the top of the bed.
- por_eps: tolerance for achieving por_target.
- nBeads: number of beads to slice starting from the bottom of the packed bed. If negative, the zBot/zTop values are used to slice the packing data.
- zBot/zTop: slice limits for the packing data on bead centers.
- refBeadSize: <min|avg|max> reference bead size used while scaling mesh elements in beads
- refBeadRadius: <double> same as above, more precise control.
- rCylDelta: Extra radius added to the cylinder to avoid overlap with beads.
- inlet/outlet: additional space at inlet and outlet. Considers 1 unit = bead diameter.
- geomInfile/geomOutfile: geometry file to input/output while meshing to save time on OCC operations.
- meshSizeMethod: <0|1> global points vs scaled fields
- lc_beads/lc_out/lc_bridge: mesh sizes
- fieldThresholdMinFactor/fieldThresholdMaxFactor: lc_beads vs lc_out threshold along the bead radius. [0 to 1]
- outputFragments: <0|1> output fragments or not
- fragment: use the fragment operation: <0|1> (necessary)
- dryRun: <0|1>
- reduced: multiplicative shrink factor for beads in place (rFactor)
- bridged/capped: <double> relativeBridgeRadius, while bridgeOffsetRatio and bridgeTol are set automatically
- and other GMSH settings such as Mesh.NumThreads that are sent to GMSH.


# Application Design
- These are incomplete notes on how the program works.
- Parameters class holds information about the model. Taken from the input config file and then modified based on those inputs. Defaults are specified in the header file. Some defaults are set during runtime.
- Bead class holds information about a single individual bead: x,y,z,r etc. The class can then apply scale and translate operations on the data. And check for neighbouring beads.
- PackedBed class 
    - reads input file
    - slices packing based on input
    - shrinks beads
    - transforms beads (scale and offset)
    - computes new bounds
    - computes porosities
- Model class
    - create/save/load geometry (using OCCT)
        - create beads/bridges/cylinder
        - create mesh fields for entities
    - generate/save mesh 
        - appropriate boolean operations
        - set named groups
        - meshing

## Notes
- xCyl, yCyl, rCyl are modified in PackedBed::transformBeads()
- if `nBeads > 0`, column length is not restricted
- if `nBeads < 0`, column length AND bead packing is calculated by using zTop and zBot. This is to ensure similar columns between mono/poly packings.
- it is not possible to control the column size if nBeads is positive
- inlet/outlet volumes are added w.r.t. zTop/zBot, not zMax/zMin.
- Porosity control is only done by REMOVING beads. Adding beads is not yet supported. There are three methods of porosity control. 
    - Remove beads by closest radius. This method finds the closest value of the radius to be removed, and removes it. This will be more accurate, but will also remove beads from the bulk of the packed bed. Not recommended. 
    - Remove beads from the top end (default). This method works best as it doesn't modify the packing in the bulk of the packed bed.
    - Remove closest beads from an end zone with length == radius_max (Best)
