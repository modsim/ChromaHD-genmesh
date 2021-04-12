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
    column = Column(geom->dimTagsFragmented, prm, prm->periodic);

    // remove z from periodic string
    std::string periodicInOut = prm->periodic;
    periodicInOut.erase(std::remove(periodicInOut.begin(), periodicInOut.end(), 'z'), periodicInOut.end());

    if(prm->periodicInlet > 0)
    {
        columnInlet = Column(geom->dimTagsFragmentedPeriodicInlet, prm, periodicInOut);
        std::cout << "Linking inlet and column" << std::endl;
        columnInlet.linkPeriodicZ(column);
    }

    if(prm->periodicOutlet > 0)
    {
        columnOutlet = Column(geom->dimTagsFragmentedPeriodicOutlet, prm, periodicInOut);
        std::cout << "Linking column and outlet" << std::endl;
        column.linkPeriodicZ(columnOutlet);
    }

    // TODO: re-implement
    /* if geometry is imported instead: */
    /*     std::cout << "Importing geometry: " << prm->geomInfile << "... " << std::flush; */
    /*     factory::importShapes("geometries/" + prm->geomInfile, dimTagsFragmented); */
    /*     std::cout << "done!" << std::endl; */

    /*     // Synchronize gmsh model with geometry kernel. */
    /*     std::cout << "synchronizing... " << std::flush; */
    /*     factory::synchronize(); */
    /*     std::cout << "done!" << std::endl; */

    /*     createNamedGroups(dimTagsFragmented, prm->containerShape); */

    //if write geometry
    if(!prm->geomOutfile.empty())
        gmsh::write(prm->outpath + "/"+ prm->geomOutfile );

    // Synchronize gmsh model with geometry kernel.
    std::cout << "synchronizing... " << std::flush;
    factory::synchronize();
    std::cout << "done!" << std::endl;
}


Model::~Model()
{

}


void Model::mesh(Parameters * prm)
{
    if(prm->dryRun)
        return;

    // For debug purposes.
    /* for(int i=1; i<prm->MeshGenerate; i++) */
    /* { */
    /*     model::mesh::generate(i); */
    /*     gmsh::write(prm->outpath + "/"+ std::to_string(i) + "D_" + outfile ); */
    /* } */

    // Generate mesh
    model::mesh::generate(prm->MeshGenerate);


}

void Model::write(std::string outfile, Parameters * prm)
{
    /* // Output Full mesh */
    /* gmsh::write(prm->outpath + "/" + outfile-); */

    std:: string basefilename = remove_extension(outfile);
    std:: string extension = get_extension(outfile);

    column.write(prm->outpath, basefilename + "_column",  extension );
    column.writeFragments(prm->outpath, basefilename + "_column",  extension );

    if(prm->periodicInlet > 0)
    {
        columnInlet.write(prm->outpath, basefilename + "_inlet",  extension );
        columnInlet.writeFragments(prm->outpath, basefilename + "_inlet",  extension );

    }

    if(prm->periodicOutlet > 0)
    {
        columnOutlet.write(prm->outpath, basefilename + "_outlet",  extension );
        columnOutlet.writeFragments(prm->outpath, basefilename + "_outlet",  extension );
    }

}
