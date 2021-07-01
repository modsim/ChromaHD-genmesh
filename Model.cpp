/*
 * File: Model.cpp
 * Created by: Rao, Jayghosh Subodh
 * Created on: Sat 06 Apr 2019 01:03:26 PM CEST
 *
 * This class is modeled after and copied from
 * BeadsInside from the old mesher.
 *
 */

#include "Model.h"
#include "Files.h"

#include <float.h>
#include <gmsh.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>


// uncomment to disable assert()
// #define NDEBUG
#include <cassert>

#define PI 3.1415926535897932

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

Model::Model(Parameters * prm, Geometry * geom)
{
    column = Column(geom->dt_fragmented, prm, prm->periodic);

    // remove z from periodic string
    // So that inlet/outlet sections don't have internal z periodicity.
    std::string periodicInOut = prm->periodic;
    periodicInOut.erase(std::remove(periodicInOut.begin(), periodicInOut.end(), 'z'), periodicInOut.end());

    if(prm->periodicInlet > 0)
    {
        columnInlet = Column(geom->dt_fragmentedPeriodicInlet, prm, periodicInOut);
        std::cout << "Linking inlet and column... " << std::endl;
        columnInlet.linkPeriodicZ(column);
    }

    if(prm->periodicOutlet > 0)
    {
        columnOutlet = Column(geom->dt_fragmentedPeriodicOutlet, prm, periodicInOut);
        std::cout << "Linking column and outlet... " << std::endl;
        column.linkPeriodicZ(columnOutlet);
    }

    // TODO: re-implement
    /* if geometry is imported instead: */
    /*     std::cout << "Importing geometry: " << prm->geomInfile << "... " << std::flush; */
    /*     factory::importShapes("geometries/" + prm->geomInfile, dt_fragmented); */
    /*     std::cout << "done!" << std::endl; */

    /*     // Synchronize gmsh model with geometry kernel. */
    /*     std::cout << "Synchronizing... " << std::flush; */
    /*     factory::synchronize(); */
    /*     std::cout << "done!" << std::endl; */

    /*     createNamedGroups(dt_fragmented, prm->containerShape); */

    //if write geometry
    if(!prm->geomOutfile.empty())
        gmsh::write(prm->outpath + "/"+ prm->geomOutfile );

    // Synchronize gmsh model with geometry kernel.
    std::cout << "Synchronizing... " << std::flush;
    factory::synchronize();
    std::cout << "done!" << std::endl;
}


Model::~Model()
{

}


void Model::mesh(std::string outfile, Parameters * prm)
{


    std::vector<std::pair<int,int>> dummy;
    model::getEntities(dummy, 3);
    std::cout << "Number of 3D objects: " << dummy.size() << std::endl;
    model::getEntities(dummy, 2);
    std::cout << "Number of 2D objects: " << dummy.size() << std::endl;
    model::getEntities(dummy, 1);
    std::cout << "Number of 1D objects: " << dummy.size() << std::endl;
    model::getEntities(dummy, 0);
    std::cout << "Number of 0D objects: " << dummy.size() << std::endl;

    if(prm->dryRun)
        return;

    std:: string basefilename = remove_extension(outfile);
    std:: string extension = get_extension(outfile);

    long time_0 = gmsh::logger::getWallTime();

    for(int i=1; i<prm->MeshGenerate; i++)
    {

        std::cout << "Meshing " << i << "D objects..." << std::endl;
        long start = gmsh::logger::getWallTime();

        model::mesh::generate(i);

        long duration = gmsh::logger::getWallTime() - start;

        gmsh::logger::write("Wall time for " + std::to_string(i) + "D mesh: " + std::to_string(duration) + " s (" + std::to_string((double)duration/3600) + " h)", "info");
        gmsh::logger::write("CPU  time for " + std::to_string(i) + "D mesh: " + std::to_string(gmsh::logger::getCpuTime()) + " s", "info");

        gmsh::write(prm->outpath + "/"+ basefilename + "_" + std::to_string(i) + "D" + extension );


    }


    long start = gmsh::logger::getWallTime();

    std::cout << "Meshing " << prm->MeshGenerate << "D objects..." << std::endl;
    model::mesh::generate(prm->MeshGenerate);

    long duration = gmsh::logger::getWallTime() - start;

    long durationAll = gmsh::logger::getWallTime() - time_0;

    gmsh::logger::write("Wall time: " + std::to_string(duration) + " s (" + std::to_string((double)duration/3600) + " h)", "info");
    gmsh::logger::write("CPU  time: " + std::to_string(gmsh::logger::getCpuTime()) + " s", "info");

    std::cout << "Total meshing wall time: " << durationAll << "s (" << durationAll/3600 << " h)" << std::endl;

}

void Model::write(std::string outfile, Parameters * prm)
{
    if(prm->dryRun)
        return;

    std:: string basefilename = remove_extension(outfile);
    std:: string extension = get_extension(outfile);

    // Output Full mesh
    gmsh::write(prm->outpath + "/" + basefilename + "_FULL" + extension);

    std::cout << "Writing main column unit..." << std::endl;
    column.write(prm->outpath, basefilename + "_column",  extension );
    column.meshVolumes(prm);
    column.writeFragments(prm->outpath, basefilename + "_column",  extension );

    if(prm->periodicInlet > 0)
    {
        std::cout << "Writing periodic inlet unit..." << std::endl;
        columnInlet.write(prm->outpath, basefilename + "_inlet",  extension );
        columnInlet.meshVolumes(prm);
        columnInlet.writeFragments(prm->outpath, basefilename + "_inlet",  extension );
    }

    if(prm->periodicOutlet > 0)
    {
        std::cout << "Writing periodic outlet unit..." << std::endl;
        columnOutlet.write(prm->outpath, basefilename + "_outlet",  extension );
        columnInlet.meshVolumes(prm);
        columnOutlet.writeFragments(prm->outpath, basefilename + "_outlet",  extension );
    }

}
