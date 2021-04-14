/**
Desc    : Tool to generate meshes from a given packing.xyzd file
Author  : Rao, Jayghosh Subodh
Created : Thu 04 Apr 2019 03:53:10 PM CEST
  **/

#include "Parameters.h"
#include "PackedBed.h"
#include "Model.h"
#include "version.h"
#include "Files.h"
#include <gmsh.h>


namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

int main(int argc, char** argv) {

    std::cout << "# GIT STATE: " << GITCOMMIT << " " << GITSTATE << std::endl;
    std::cout << "# GMSH VERSION: " << GMSHVERSION << std::endl;

    gmsh::initialize();

    std::string outfile, infile, outpath;

    if (argc == 1)
    {
        // If program is called without arguments, use default inputs and outputs

        infile = "default.in";
        outfile = "output.msh2";
        std::cout << "# Will store mesh in output.msh!" << std::endl;
    }
    else if (argc == 2)
    {
        // If one argument is provided,
        infile = std::string(argv[1]);
        outfile = "output.msh2";
    }
    else
    {
        //If both input and output files are provided as commandline arguments
        infile = std::string(argv[1]);
        outfile = std::string(argv[2]);
    }

    try{

        std::cout << std::scientific;
        Parameters * prm = new Parameters(infile);


        prm->outpath += "/" + remove_extension(outfile);
        create_directory(prm->outpath);

        PackedBed * packedBed = new PackedBed(prm);

        Geometry * geom = new Geometry(prm, packedBed);
        Model * defaultModel = new Model(prm, geom);

        long start = gmsh::logger::getWallTime();

        defaultModel->mesh(outfile, prm);
        defaultModel->write(outfile, prm);

        long duration = gmsh::logger::getWallTime() - start;
        gmsh::logger::write("Wall time: " + std::to_string(duration) + " s (" + std::to_string((double)duration/3600) + " h)", "info");
        gmsh::logger::write("CPU  time: " + std::to_string(gmsh::logger::getCpuTime()) + " s", "info");

        delete prm;
        delete packedBed;
        delete geom;
        delete defaultModel;

    }catch(mixd::MixdException e)
    { std::cout << e.msg() << std::endl; return 1; }

    gmsh::finalize();


}
