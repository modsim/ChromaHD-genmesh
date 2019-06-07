/**
Desc    : Tool to generate meshes from a given packing.xyzd file
Author  : Rao, Jayghosh Subodh
Created : Thu 04 Apr 2019 03:53:10 PM CEST
  **/

#include "Parameters.h"
#include "PackedBed.h"
#include<gmsh.h>

/* #include "TestBed.h" */

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

    /* double lc = 0.06; */

    try{

        Parameters * prm = new Parameters("default.in");
        prm->print();

        gmsh::option::setNumber("General.Terminal", 5);
        gmsh::option::setNumber("General.NumThreads", 8);
        gmsh::option::setNumber("Geometry.OCCParallel", 1);

        gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1); //Default = 1 (longest)
        gmsh::option::setNumber("Mesh.CharacteristicLengthMin", 0  ); //Default = 0
        gmsh::option::setNumber("Mesh.CharacteristicLengthFromCurvature", 1); //Default = 0
        gmsh::option::setNumber("Mesh.CharacteristicLengthFromPoints", 1); //Default = 1
        /* gmsh::option::setNumber("Mesh.MinimumCirclePoints", 7); //Default = 7 */
        /* gmsh::option::setNumber("Mesh.BoundaryLayerFanPoints", 5); //Default = 5 */

        gmsh::option::setNumber("Mesh.OptimizeNetgen", prm->MeshOptimizeNetgen); //Default = 0
        gmsh::option::setNumber("Mesh.RefineSteps", prm->MeshRefineSteps); //Default = 10

        //1: MeshAdapt | 2: Auto | 5: Delaunay | 6: Frontal | 7: BAMG | 8: DelQuad
        gmsh::option::setNumber("Mesh.Algorithm", prm->MeshAlgorithm); //Default = 2
        //1: Delaunay | 4: Frontal | 5: Frontal Delaunay | 6: Frontal Hex | 7: MMG3D | 9: RTree | 10: HXT
        gmsh::option::setNumber("Mesh.Algorithm3D", prm->MeshAlgorithm3D); //Default = 1

        PackedBed * packedBed = new PackedBed(prm);

        packedBed->createGeometry();
        packedBed->createBridge(prm->db_dp, prm->bridgeTol);
        /* packedBed->printPacking(); */
        packedBed->mesh(outfile);
        delete prm;
        delete packedBed;

        /* TestBed * testBed = new TestBed(); */
        /* testBed->makeCyl(); */
        /* testBed->makeBead(); */
        /* testBed->fuse(); */
        /* testBed->mesh(); */

        /* int b1 = factory::addDisk(-1,0,0, 1.000, 1.00); */
        /* int b2 = factory::addDisk( 1,0,0, 1.000, 1.00); */
        /* int r1 = factory::addRectangle(-3,-3,0, 6,6); */
        /* std::vector<std::pair<int, int> > ov; */
        /* std::vector<std::vector<std::pair<int, int> > > ovv; */
        /* factory::fragment({{2,r1}}, {{2,b1},{2,b2}}, ov, ovv); */
        /* factory::synchronize(); */
        /* model::getEntities(ov, 0); */
        /* model::mesh::setSize(ov, 0.1); */
        /* model::mesh::generate(2); */
        /* gmsh::write("twospheres.msh"); */



    }catch(mixd::MixdException e)
    { std::cout << e.msg() << std::endl; return 1; }


    gmsh::finalize();


}
