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

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>

#define PI 3.1415926535897932

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

Model::Model(Parameters * prm, Geometry * geom)
{
    /* model::add("Model"); */
}


Model::~Model()
{

}

void Model::mesh(std::string outfile, Parameters * prm, Geometry * geom)
{
    if(prm->dryRun)
        return;

    // Store dimTags and dimTagsMap
    std::vector<std::pair<int, int> > dimTagsFragmented  = geom->dimTagsFragmented;

    std::vector<std::pair<int, int> > dimTagsDummy;
    std::vector<std::vector<std::pair<int, int> > > ovv;

    // ============================
    // Named Physical Groups

    // if not reading geometry
    if (prm->geomInfile.empty())
    {
        if (prm->fragment == 1)
            createNamedGroups(dimTagsFragmented, prm->containerShape);

        if ((prm->periodic == "xy") || (prm->periodic == "xyz"))
            setupPeriodicSurfaces(prm, geom);
    }
    else
    {
        std::cout << "Importing geometry: " << prm->geomInfile << "... " << std::flush;
        factory::importShapes("geometries/" + prm->geomInfile, dimTagsFragmented);
        std::cout << "done!" << std::endl;

        // Synchronize gmsh model with geometry kernel.
        std::cout << "synchronizing... " << std::flush;
        factory::synchronize();
        std::cout << "done!" << std::endl;

        createNamedGroups(dimTagsFragmented, prm->containerShape);
    }

    std::cout << std::endl;

    //if write geometry
    if(!prm->geomOutfile.empty())
        gmsh::write(prm->outpath + "/"+ prm->geomOutfile );

    // Synchronize gmsh model with geometry kernel.
    std::cout << "synchronizing... " << std::flush;
    factory::synchronize();
    std::cout << "done!" << std::endl;

    for(int i=1; i<prm->MeshGenerate; i++)
    {
        model::mesh::generate(i);
        gmsh::write(prm->outpath + "/"+ std::to_string(i) + "D_" + outfile );

    }
    // Generate mesh
    model::mesh::generate(prm->MeshGenerate);

    // Output Full mesh
    gmsh::write(prm->outpath + "/"+ outfile);

    //TODO: Calculate surface areas as well
    //TODO: Separate out into a function

    std::cout << std::endl;

    // plugin volume
    std::cout << "Calculating mesh volumes. Note: Multiply outputs by " << pow(prm->MeshScalingFactor, 3) << std::endl;
    meshVolumes();


    if (prm->outputFragments == 1)
    {

        model::removePhysicalGroups();

        if (prm->NamedInterstitialVolume)
        {
            model::addPhysicalGroup(3,tVInt,15);
            model::setPhysicalName(3,15,"interstitialVolume");
            gmsh::write(prm->outpath + "/" + outfile + "_interstitial_volume." + prm->fragmentFormat);
            model::removePhysicalGroups();
        }


        if (prm->NamedOuterSurface)
        {
            model::addPhysicalGroup(2, tSInlet  , 11      );
            model::addPhysicalGroup(2, tSOutlet , 12      );
            model::addPhysicalGroup(2, tSWall   , 13      );
            model::setPhysicalName (2, 11       , "inlet" );
            model::setPhysicalName (2, 12       , "outlet");
            model::setPhysicalName (2, 13       , "wall"  );
            gmsh::write(prm->outpath + "/" + outfile + "_interstitial_surfaces." + prm->fragmentFormat);
            model::removePhysicalGroups();

            model::addPhysicalGroup(2, tSInlet  , 11      );
            model::setPhysicalName (2, 11       , "inlet" );
            gmsh::write(prm->outpath + "/" + outfile + "_inlet_surface." + prm->fragmentFormat);
            model::removePhysicalGroups();

            model::addPhysicalGroup(2, tSOutlet , 12      );
            model::setPhysicalName (2, 12       , "outlet");
            gmsh::write(prm->outpath + "/" + outfile + "_outlet_surface." + prm->fragmentFormat);
            model::removePhysicalGroups();

            model::addPhysicalGroup(2, tSWall   , 13      );
            model::setPhysicalName (2, 13       , "wall"  );
            gmsh::write(prm->outpath + "/" + outfile + "_wall_surface." + prm->fragmentFormat);
            model::removePhysicalGroups();
        }


        if (prm->NamedBeadVolume)
        {
            model::addPhysicalGroup(3,tVBeads,16 );
            model::setPhysicalName(3,16,"beadVolume");
            gmsh::write(prm->outpath + "/" + outfile + "_beads_volumes." + prm->fragmentFormat);
            model::removePhysicalGroups();
        }


        if (prm->NamedBeadSurface)
        {
            model::addPhysicalGroup(2,tSBeads,14);
            model::setPhysicalName(2,14, "beadSurface");
            gmsh::write(prm->outpath + "/" + outfile + "_beads_surfaces." + prm->fragmentFormat);
        }

    }


}

void Model::createNamedGroups(std::vector<std::pair<int,int>> bv, int containerShape)
{
    if (containerShape == -1)
    {
        std::cout << ">>> Not creating Named Groups!!" << std::endl;
        return;
    }

    // bv is the vector storing dimtags after fragmentation.
    // the last entry should be the interstitial volume
    /* std::vector<std::pair<int,int>> cv; */
    std::vector<std::pair<int,int>> dimTagsBeadSurface;
    std::vector<std::pair<int,int>> dimTagsOuterSurface;

    std::cout << "Listing Interstitial Volume... " << std::flush;
    tVInt.push_back(bv.back().second);
    bv.pop_back();
    std::cout << "done!" << std::endl;

    // Extract surfaces from interstitial volume
    std::cout << "Listing Outer Surfaces... " << std::flush;
    model::getBoundary({{3,tVInt[0]}}, dimTagsOuterSurface, false, false, false);
    model::getBoundary(bv, dimTagsBeadSurface, false, false, false);
    if(containerShape == 0)
    {
        tSWall   = {dimTagsOuterSurface[0].second};
        tSOutlet = {dimTagsOuterSurface[1].second};
        tSInlet  = {dimTagsOuterSurface[2].second};

        // NOTE: What would be the normal for a cylinder??
    }
    else if (containerShape == 1)
    {
        std::vector<double> points, parametricCoord, curvatures, normals;
        std::vector<std::pair <int, int>> dimTagsBound;

        // Interstitial surfaces
        for (std::vector<std::pair<int, int>>::iterator iter = dimTagsOuterSurface.begin(); iter != dimTagsOuterSurface.end(); iter++)
        {

            // Get the curvatures of the surfaces,
            // Select planar surfaces
            // Get normals of planes to match surface position to tag

            model::getBoundary({(*iter)}, dimTagsBound, false, false, true);
            for (std::vector<std::pair<int, int>>::iterator iteraa = dimTagsBound.begin(); iteraa != dimTagsBound.end(); iteraa++)
            {
                model::getValue((*iteraa).first, (*iteraa).second, {}, points);
            }
            model::getParametrization((*iter).first, (*iter).second, points, parametricCoord);
            model::getCurvature((*iter).first, (*iter).second, parametricCoord, curvatures);

            //length of curvatures was always 1
            if (curvatures[0] == 0)
            {
                model::getNormal((*iter).second, parametricCoord, normals);

                if (int(normals[0]) == -1)
                    tXLeftWallInt.push_back((*iter).second);
                else if (int(normals[1]) == -1)
                    tYLeftWallInt.push_back((*iter).second);
                else if (int(normals[2]) == -1)
                    tZLeftWallInt.push_back((*iter).second);
                else if (int(normals[0]) == 1)
                    tXRightWallInt.push_back((*iter).second);
                else if (int(normals[1]) == 1)
                    tYRightWallInt.push_back((*iter).second);
                else if (int(normals[2]) == 1)
                    tZRightWallInt.push_back((*iter).second);

            }

        }

        // Bead Surfaces
        // For cut beads, each surface is different. nBeads != nSurfaces (I think)
        for (std::vector<std::pair<int, int>>::iterator iter = dimTagsBeadSurface.begin(); iter != dimTagsBeadSurface.end(); iter++)
        {
            /*
             * Get the curvatures of the surfaces,
             * Select planar surfaces
             * Get normals of planes to match surface position to tag
             */
            model::getBoundary({(*iter)}, dimTagsBound, false, false, true); // Get points dimtags
            /* std::cout << "Getting points for " << (*iter).second << std::endl; */
            for (std::vector<std::pair<int, int>>::iterator iteraa = dimTagsBound.begin(); iteraa != dimTagsBound.end(); iteraa++)
            {
                model::getValue((*iteraa).first, (*iteraa).second, {}, points); // get points coords
                /* std::cout << "[" << (*iteraa).first << ", "<< (*iteraa).second << "]: " << std::flush; */
                /* for (auto itp = points.begin(); itp != points.end(); itp++) */
                /*     std::cout << *itp << "  " << std::flush; */
                /* std::cout << std::endl; */
            }
            model::getParametrization((*iter).first, (*iter).second, points, parametricCoord);
            model::getCurvature((*iter).first, (*iter).second, parametricCoord, curvatures);

            //length of curvatures was always 1
            if (curvatures[0] == 0)
            {
                model::getNormal((*iter).second, parametricCoord, normals);

                if (int(normals[0]) == -1)
                    tXLeftWallBead.push_back((*iter).second);
                else if (int(normals[1]) == -1)
                    tYLeftWallBead.push_back((*iter).second);
                else if (int(normals[2]) == -1)
                    tZLeftWallBead.push_back((*iter).second);
                else if (int(normals[0]) == 1)
                    tXRightWallBead.push_back((*iter).second);
                else if (int(normals[1]) == 1)
                    tYRightWallBead.push_back((*iter).second);
                else if (int(normals[2]) == 1)
                    tZRightWallBead.push_back((*iter).second);
            }

        }

        tSWall.reserve(tSWall.size() + tXLeftWallInt.size() + tXRightWallInt.size() + tYLeftWallInt.size() + tYRightWallInt.size());

        tSWall.insert(tSWall.end(), tXLeftWallInt.begin(), tXLeftWallInt.end());
        tSWall.insert(tSWall.end(), tXRightWallInt.begin(), tXRightWallInt.end());
        tSWall.insert(tSWall.end(), tYLeftWallInt.begin(), tYLeftWallInt.end());
        tSWall.insert(tSWall.end(), tYRightWallInt.begin(), tYRightWallInt.end());

        tSWall.reserve(tSWall.size() + tXLeftWallBead.size() + tXRightWallBead.size() + tYLeftWallBead.size() + tYRightWallBead.size());

        tSWall.insert(tSWall.end(), tXLeftWallBead.begin(), tXLeftWallBead.end());
        tSWall.insert(tSWall.end(), tXRightWallBead.begin(), tXRightWallBead.end());
        tSWall.insert(tSWall.end(), tYLeftWallBead.begin(), tYLeftWallBead.end());
        tSWall.insert(tSWall.end(), tYRightWallBead.begin(), tYRightWallBead.end());

        // TODO: In z-periodic cases, do I want to include the cut bead surfaces as inlet/outlet?
        // Also: How do I convert data from outlet -> inlet and provide it to xns.in?
        // so far, rngdexp... might need to check the option to specify node-by-node dof constraint
        tSInlet =  tZLeftWallInt ;
        tSOutlet=  tZRightWallInt ;
    }
    std::cout << "done!" << std::endl;

    // The remaining entries in bv are bead volumes
    std::cout << "Listing Bead Volumes... " << std::flush;
    for (auto it = bv.begin(); it != bv.end(); it++)
    {
        tVBeads.push_back((*it).second);
    }
    std::cout << "done!" << std::endl;

    std::cout << "Listing Bead Surfaces... " << std::flush;
    for ( std::vector<std::pair< int , int>>::iterator it = dimTagsBeadSurface.begin(); it != dimTagsBeadSurface.end(); it++   )
    {
        tSBeads.push_back((*it).second);
    }
    std::cout << "done!" << std::endl;

    /*
     * Remove bead surfaces that are also part of the wall
     * from being categorized as bead surfaces.
     * These 'tSBeads' will go on to be doubled in `gmsh2mixdv2 -d 4`
     */
    for (auto it : tXLeftWallBead)
    {
        auto it2 = std::find(tSBeads.begin(), tSBeads.end(), it);
        if (it2 != tSBeads.end()) tSBeads.erase(it2);
    }
    for (auto it : tXRightWallBead)
    {
        auto it2 = std::find(tSBeads.begin(), tSBeads.end(), it);
        if (it2 != tSBeads.end()) tSBeads.erase(it2);
    }
    for (auto it : tYLeftWallBead)
    {
        auto it2 = std::find(tSBeads.begin(), tSBeads.end(), it);
        if (it2 != tSBeads.end()) tSBeads.erase(it2);
    }
    for (auto it : tYRightWallBead)
    {
        auto it2 = std::find(tSBeads.begin(), tSBeads.end(), it);
        if (it2 != tSBeads.end()) tSBeads.erase(it2);
    }

    std::cout << "Adding Physical Groups... " << std::flush;
    model::addPhysicalGroup(2, tSInlet , 1 );
    model::addPhysicalGroup(2, tSOutlet, 2 );
    model::addPhysicalGroup(2, tSWall  , 3 );
    model::addPhysicalGroup(2, tSBeads , 4 );
    model::addPhysicalGroup(3, tVInt   , 5 );
    model::addPhysicalGroup(3, tVBeads , 6 );
    std::cout << "done!" << std::endl;

    std::cout << "Setting Physical Names... " << std::flush;
    model::setPhysicalName(2,1,"inlet");
    model::setPhysicalName(2,2,"outlet");
    model::setPhysicalName(2,3,"wall");
    model::setPhysicalName(2,4,"beadSurface");
    model::setPhysicalName(3,5,"interstitialVolume");
    model::setPhysicalName(3,6,"beadVolume");
    std::cout << "done!" << std::endl;

}

void Model::setupPeriodicSurfaces(Parameters * prm, Geometry * geom)
{
    std::cout << "Setting up periodic surfaces..." << std::flush;

    assert(tXLeftWallInt.size()  == tXRightWallInt.size());
    assert(tYLeftWallInt.size()  == tYRightWallInt.size());
    assert(tXLeftWallBead.size() == tXRightWallBead.size());
    assert(tYLeftWallBead.size() == tYRightWallBead.size());

    if (prm->periodic == "xyz")
    {
        assert(tZLeftWallInt.size()  == tZRightWallInt.size());
        assert(tZLeftWallBead.size() == tZRightWallBead.size());
    }

    double dx, dy, dz;

    dx = geom->dx;
    dy = geom->dy;
    dz = geom->dz;

    std::vector<double> affineTranslationX = {1, 0, 0, dx, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    std::vector<double> affineTranslationY = {1, 0, 0, 0, 0, 1, 0, dy, 0, 0, 1, 0, 0, 0, 0, 1};
    std::vector<double> affineTranslationZ = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, dz, 0, 0, 0, 1};

    /* std::cout << "periodic X wall interstitial surfaces: " << tXLeftWallInt.size() << std::endl; */
    /* std::cout << "periodic Y wall interstitial surfaces: " << tYLeftWallInt.size() << std::endl; */
    /* std::cout << "periodic Z wall interstitial surfaces: " << tZLeftWallInt.size() << std::endl; */

    matchPeriodicSurfaces(tXLeftWallInt, tXRightWallInt, 0, affineTranslationX);
    matchPeriodicSurfaces(tYLeftWallInt, tYRightWallInt, 1, affineTranslationY);
    if (prm->periodic == "xyz")
        matchPeriodicSurfaces(tZLeftWallInt, tZRightWallInt, 2, affineTranslationZ);

    /* std::cout << "periodic X wall bead surfaces: " << tXLeftWallBead.size() << std::endl; */
    /* std::cout << "periodic Y wall bead surfaces: " << tYLeftWallBead.size() << std::endl; */
    /* std::cout << "periodic Z wall bead surfaces: " << tZLeftWallBead.size() << std::endl; */

    matchPeriodicSurfaces(tXLeftWallBead, tXRightWallBead, 0, affineTranslationX);
    matchPeriodicSurfaces(tYLeftWallBead, tYRightWallBead, 1, affineTranslationY);
    if (prm->periodic == "xyz")
        matchPeriodicSurfaces(tZLeftWallBead, tZRightWallBead, 2, affineTranslationZ);

    std::cout << "done!" << std::endl;
}

void Model::matchPeriodicSurfaces(std::vector<int>& ltags, std::vector<int>& rtags, int per_dir, std::vector<double> affineTranslation)
{
    /* model::mesh::setPeriodic(2, tXRightWallBead, tXLeftWallBead, affineTranslationX); */
    // Unfortunately the order of the surfaces in the vectors matters when providing them in groups to the setPeriodic function.
    // So we match the surfaces individually instead in this function.
    // For each surface, find center of bbox. Project that along per_dir and find matching centers on either surface.
    // If match is found, setPeriodic

    std::cout << "genmesh: Matching periodic surfaces in " << per_dir << " direction..." << std::endl;

    double eps = 1e-10;

    for (auto itl = ltags.begin(); itl != ltags.end(); itl++)
    {
        double xlmin, ylmin, zlmin, xlmax, ylmax, zlmax;
        model::getBoundingBox(2, *itl, xlmin, ylmin, zlmin, xlmax, ylmax, zlmax);
        std::vector<double> lcenter = { (xlmin + xlmax) / 2, (ylmin + ylmax) / 2, (zlmin + zlmax) / 2};

        for (auto itr = rtags.begin(); itr != rtags.end(); itr++)
        {
            double xrmin, yrmin, zrmin, xrmax, yrmax, zrmax;
            model::getBoundingBox(2, *itr, xrmin, yrmin, zrmin, xrmax, yrmax, zrmax);
            std::vector<double> rcenter = { (xrmin + xrmax) / 2, (yrmin + yrmax) / 2, (zrmin + zrmax) / 2};

            double delta = 0.0;
            for (int dir = 0; dir < 3; dir++)
            {
                if (dir != per_dir)
                {
                    delta += (lcenter.at(dir) - rcenter.at(dir)) * (lcenter.at(dir) - rcenter.at(dir));
                }
            }

            delta = sqrt(delta);
            if (delta < eps)
            {
                /* std::cout << "Matched surfaces: " << *itl << " and " << *itr  << " with delta: " << delta << std::endl; */

                /* for (auto itlc = lcenter.begin(); itlc != lcenter.end(); itlc++) */
                /*     std::cout << *itlc << " " << std::flush; */
                /* std::cout << std::endl; */
                /* for (auto itrc = rcenter.begin(); itrc != rcenter.end(); itrc++) */
                /*     std::cout << *itrc << " " << std::flush; */
                /* std::cout << std::endl; */

                model::mesh::setPeriodic(2, {*itr}, {*itl}, affineTranslation);
                break;
            }
        }

    }

}

void Model::meshVolumes()
{
    std::cout << "[ Column Volume ]" << std::endl;
    gmsh::plugin::setNumber("MeshVolume", "PhysicalGroup", -1);
    gmsh::plugin::setNumber("MeshVolume", "Dimension", 3);
    gmsh::plugin::run("MeshVolume");
    std::cout << std::endl;

    std::cout << "[ Interstitial Volume ]" << std::endl;
    gmsh::plugin::setNumber("MeshVolume", "PhysicalGroup", 5);
    gmsh::plugin::setNumber("MeshVolume", "Dimension", 3);
    gmsh::plugin::run("MeshVolume");
    std::cout << std::endl;

    std::cout << "[ Bead Volume ]" << std::endl;
    gmsh::plugin::setNumber("MeshVolume", "PhysicalGroup", 6);
    gmsh::plugin::setNumber("MeshVolume", "Dimension", 3);
    gmsh::plugin::run("MeshVolume");
    std::cout << std::endl;

}
