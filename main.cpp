/**
Desc    : Tool to generate meshes from a given packing.xyzd file
Author  : Rao, Jayghosh Subodh
Created : Thu 04 Apr 2019 03:53:10 PM CEST
  **/

#include "Parameters.h"
#include "PackedBed.h"
#include "version.h"
#include<gmsh.h>

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

int main(int argc, char** argv) {

    std::cout << "# GIT STATE: " << GITCOMMIT << " " << GITSTATE << std::endl;

    gmsh::initialize();

    std::string outfile, infile;

    if (argc == 1)
    {
        infile = "default.in";
        outfile = "output.msh2";
        std::cout << "# Will store mesh in output.msh!" << std::endl;
    }
    else if (argc == 2)
    {
        infile = std::string(argv[1]);
        outfile = "output.msh2";
    }
    else
    {
        infile = std::string(argv[1]);
        outfile = std::string(argv[2]);
    }
    try{


        Parameters * prm = new Parameters(infile);
        prm->setGMSHOptions();

        PackedBed * packedBed = new PackedBed(prm);

        long start = gmsh::logger::time();

        packedBed->createGeometry();
        packedBed->mesh(outfile);

        long duration = gmsh::logger::time() - start;

        gmsh::logger::write("Wall time: " + std::to_string(duration) + " s", "info");
        gmsh::logger::write("CPU  time: " + std::to_string(gmsh::logger::cputime()) + " s", "info");

        delete prm;
        delete packedBed;

    }catch(mixd::MixdException e)
    { std::cout << e.msg() << std::endl; return 1; }

    gmsh::finalize();


}
