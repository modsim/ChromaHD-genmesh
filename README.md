genmesh

Generate meshes from packing files (xyzd). Depends on OpenCASCADE (7.3) and GMSH (v4+). 

# Installation

1. Install OpenCASCADE (v 7.30 tested)
2. Install GMSH with OpenCASCADE linked (v 4.4.1)
3. Compile genmesh with GMSH linked.

**NOTE: Newer versions of GMSH and OpenCASCADE are required for later versions of genmesh. Especially periodic cases**

Check out `install.sh` for a template on what to do.

Note: Ensure that \$LD_LIBRARY_PATH points to the OpenCASCADE libs.
Note: Ensure that \$GMSH_ROOT and \$OCCT_ROOT point to the respective install directories.
Note: Ensure that you use the same compiler version for all dependencies. 

## Snippets
```
freetype2, mesa/opengl dev packages, libXmu, libXi, libfontconfig1-dev
./configure --enable-gcc --enable-shared --enable-threads --prefix=/home/IBT/rao/tools/tcl --enable-64bit
./configure --enable-gcc --enable-shared --enable-threads --prefix=/home/IBT/rao/tools/tk --enable-64bit --with-tcl=/home/IBT/rao/tools/tcl/lib
cmake -D3RDPARTY_TCL_DIR=$HOME/tools/tcl -D3RDPARTY_TK_DIR=$HOME/tools/tk -DCMAKE_INSTALL_PREFIX=/home/IBT/rao/tools/occt-7.5.0 ..
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


# Known Issues
- Netgen optimizer crashes sometimes. (After mesh size constraints were applied to surfaces)
- Geometry.ScalingFactor doesn't work. Use preScalingFactor instead.
- poly5 fails (beads peek out of container). Just use poly-full.
    - This problem persists with poly-full for certain sections. Can be fixed by a better centering method instead of only x/y based centering
    - It can be sidestepped by running genmesh with coarse element size, and checking if any beads are meshed AFTER the cylinder/planes & have a higher entity number.
- A few features do not work in conjunction with HXT mesh algorithm. This is an upstream issue with GMSH.
- Netgen optimizer might not work well with mesh field gradients: It might make the meshes too uniform somehow
- periodic meshing fails for cases where there are no beads intersecting the cut planes (simple test cases)
- periodic meshes with z-dir cuts in beads can contain floating nodes which have 0 nmat info due to not being connected to an element. 
    - the floating nodes are found on the x/y surfaces intersecting with the z surfaces
    - Simulation could be theoretically continued if we modify the nmat file to have non-zero values at these points
- 2021-04-10 : gradient mesh size within beads only works if lc_beads < lc_out

# Keywords in config.in
Look at parameters.h for a more up-to-date set of keywords and default values

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
- autoContainment: <0|1> automatically calculate bounding box for packed bed and create container. Switch this off for periodic meshes. Off => use cyl or box keywords
- containerShape: <0|1> cylinder or box
- cyl: <xCyl yCyl rCyl> 
- box: <xMin yMin zMin xDelta yDelta zDelta> 
- translateOffsets: <auto|x y z> translate packed bed automatically or give x y z offsets
- periodicOffsets: <auto|x y> offsets for stacking the geometry in x-y directions for periodicity
- periodicInlet and periodicOutlet: <double> length of the inlet and outlet sections of the column to be "linked" periodically to main column but generated separately. 
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
- Geometry Class
    - takes PackedBed 
    - Generates geometries of beads, bridges and containers
    - Performs necessary boolean operations
    - cleans up unnecessary entities.
- Column Class
    - Takes a fragmented column dimtags
    - Finds and names surfaces appropriately
    - Sets periodic surfaces appropriately
- Model class
    - Driver class that builds columns, meshes, and writes them out

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
- For periodic meshes:
    - `autoContainment 0`
    - `containerShape 1`
    - `box <coord/dimensions>`
    - `periodicOffsets <auto | pOffX pOffY pOffZ>`
    - `periodicInlet <double>`
    - `periodicOutlet <double>`
    - mind the `translateOffsets` and actual coordinates of beads from xyzd file.
    - Ensure that cut planes are either symmetric or far away from bead surface. Essentially, in case of "xy" periodicity, ensure that there is inlet/outlet void space at the ends of the column. 
- periodic meshing only in the z-direction isn't available yet. Options are "xy" and "xyz" only
- periodicInlet and periodicOutlet keywords create a periodically linked, but separate mesh for the inlet and outlet sections of the column, that are columns in their own right. 
- periodicInlet and periodicOutlet sections derive periodicity from the main column, but discard z-periodicity.
- Memleaks seem to be in GMSH or OCCT.

## Packing Generation
This [Packing Generation](https://github.com/VasiliBaranov/packing-generation) tool can be used to generate periodic packings. As of this writing, this [issue](https://github.com/VasiliBaranov/packing-generation/issues/18) is still unresolved, and confined cylindrical packings don't work. 

More detailed information can be found in the README of the above project, but here I'll elucidate a clear workflow to generate periodic packings. 

Essentially, generating the packing involves running the program multiple times with different algorithms: FBA -> LS -> LSGD 

[FBA]
- Ensure that generation.conf and diameters.txt (optional) are in the current directory.
- Ensure that generation.conf has Generation start: 1
- Run FBA: `PackingGeneration.exe -fba |& tee fba.log`

[LS]
- Remove packing.nfo, optionally move/save all files to a different directory
- Ensure that generation.conf, diameters.txt and packing.xyzd are in the current directory.
- Ensure that generation.conf has Generation start: 0
- Run LS: `PackingGeneration.exe -ls |& tee ls.log`
- If every other line is a warning about mismatch in inner and outer diameter ratios, restart the whole process with modified nbeads and dimensions such that porosity is close to 0.4. (0.4-0.6)

[LSGD]
- Remove packing.nfo, optionally move/save all files to a different directory
- Ensure that generation.conf, diameters.txt and packing.xyzd are in the current directory.
- Ensure that generation.conf has Generation start: 0
- Run LS: `PackingGeneration.exe -lsgd |& tee lsgd.log`

[dumpy]
- Run dumpy to scale the diameters: `dumpy --dpacking packing.xyzd --nfo packing.nfo -w packing_fixed.xyzd`. (or just scale bead diameters manually in genmesh)


