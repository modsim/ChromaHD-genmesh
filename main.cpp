/**
Desc    : Tool to generate meshes from a given packing.xyzd file
Author  : Rao, Jayghosh Subodh
Created : Thu 04 Apr 2019 03:53:10 PM CEST
  **/

#include "Parameters.h"
#include "PackedBed.h"
#include<gmsh.h>

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

int main(int argc, char** argv) {

    gmsh::initialize();

    std::string outfile;

    if (argc != 2)
    {
        outfile = "output.msh2";
        std::cout << "Will store mesh in output.msh!" << std::endl;
    }
    else
        outfile = std::string(argv[1]);

    try{

        Parameters * prm = new Parameters("default.in");
        prm->print();

        gmsh::option::setNumber("General.Terminal", 1);
        gmsh::option::setNumber("General.NumThreads", 8);
        gmsh::option::setNumber("Geometry.OCCParallel", 1);

        gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1); //Default = 1 (longest)
        gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 0  ); //Default = 0
        gmsh::option::setNumber("Mesh.CharacteristicLengthFromCurvature", 1); //Default = 0
        gmsh::option::setNumber("Mesh.CharacteristicLengthFromPoints", 1); //Default = 1

        /* gmsh::option::setNumber("Mesh.MinimumCirclePoints", 7); //Default = 7 */
        /* gmsh::option::setNumber("Mesh.BoundaryLayerFanPoints", 5); //Default = 5 */

        gmsh::option::setNumber("Mesh.Optimize", prm->MeshOptimize); //Default = 1
        gmsh::option::setNumber("Mesh.OptimizeNetgen", prm->MeshOptimizeNetgen); //Default = 0
        gmsh::option::setNumber("Mesh.RefineSteps", prm->MeshRefineSteps); //Default = 10

        //1: MeshAdapt | 2: Auto | 5: Delaunay | 6: Frontal | 7: BAMG | 8: DelQuad
        gmsh::option::setNumber("Mesh.Algorithm", prm->MeshAlgorithm); //Default = 2

        //1: Delaunay | 4: Frontal | 5: Frontal Delaunay | 6: Frontal Hex | 7: MMG3D | 9: RTree | 10: HXT
        gmsh::option::setNumber("Mesh.Algorithm3D", prm->MeshAlgorithm3D); //Default = 1

        PackedBed * packedBed = new PackedBed(prm);
        packedBed->createGeometry();
        packedBed->createBridge(prm->db_dp, prm->bridgeTol);
        packedBed->mesh(outfile);
        delete prm;
        delete packedBed;

    }catch(mixd::MixdException e)
    { std::cout << e.msg() << std::endl; return 1; }

    gmsh::finalize();


}
