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

#define PI 3.1415926535897932

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

Model::Model(Parameters * prm)
{
    model::add("Model");
}


Model::~Model()
{

}

void Model::createGeometry(PackedBed * pb, Parameters * prm)
{

    std::vector<Bead *> beads = pb->beads;

    std::vector<std::pair<int, int> > ov;
    std::vector<std::pair<int, int> > cv;

    std::vector<double> vCount;


    double zCylBot = ( prm->zBot - prm->inlet  );
    double zCylTop = ( prm->zTop + prm->outlet );
    double xCyl    = prm->xCyl;
    double yCyl    = prm->yCyl;
    double rCyl    = prm->rCyl;

    int count = 0;

    int tag, ctag;
    std::cout << "Creating beads... " << std::flush;
    for (std::vector<Bead *>::iterator iter = beads.begin(); iter != beads.end(); iter++ )
    {
        double x = (*iter)->getX();
        double y = (*iter)->getY();
        double z = (*iter)->getZ();
        double r = (*iter)->getR();

        if (prm->geomInfile.empty())
        {
            tag = factory::addSphere(x, y, z, r);
            ctag = factory::addPoint(x, y, z, prm->lc_beads, -1);

            (*iter)->setTag(tag);
            (*iter)->setCTag(ctag);

            dimTagsBeads.push_back({3, tag});
            tBeadCPs.push_back(ctag);
        }

        /* model::mesh::field::add("Ball", ++count); */
        /* model::mesh::field::setNumber(count, "VIn", r/(prm->refBeadRadius) * prm->lc_beads); */
        /* model::mesh::field::setNumber(count, "VOut", prm->lc_out); */
        /* model::mesh::field::setNumber(count, "XCenter", x); */
        /* model::mesh::field::setNumber(count, "YCenter", y); */
        /* model::mesh::field::setNumber(count, "ZCenter", z); */
        /* model::mesh::field::setNumber(count, "Radius", prm->fieldExtensionFactor * r); */
        /* vCount.push_back(count); */

        model::mesh::field::add("Distance", ++count);
        model::mesh::field::setNumbers(count, "NodesList", {double(ctag)});

        model::mesh::field::add("Threshold", ++count);
        model::mesh::field::setNumber(count, "IField", count-1);
        model::mesh::field::setNumber(count, "LcMin", prm->lc_beads);
        model::mesh::field::setNumber(count, "LcMax", prm->lc_out);
        model::mesh::field::setNumber(count, "DistMin", prm->fieldThresholdMinFactor * r);
        model::mesh::field::setNumber(count, "DistMax", prm->fieldThresholdMaxFactor * r);
        vCount.push_back(count);

    }
    std::cout << "done!" << std::endl;

    int tagBridge;

    double dx, dy, dz;
    double x4, y4, z4;
    double x3, y3, z3;
    double x2, y2, z2, r2;
    double x1, y1, z1, r1;
    double dx3, dy3, dz3;
    double factor = prm->bridgeOffsetRatio;
    double dist;

    std::vector<Bead *> remBeads = beads;
    std::cout << "Creating Bridges... " << std::flush;
    for(std::vector<Bead *>::iterator iter = beads.begin(); iter != beads.end(); iter++)
    {
        remBeads.erase(remBeads.begin());
        for(std::vector<Bead *>::iterator riter = remBeads.begin(); riter != remBeads.end(); riter++)
        {

            x1 = (*iter)->getX();
            x2 = (*riter)->getX();

            y1 = (*iter)->getY();
            y2 = (*riter)->getY();

            z1 = (*iter)->getZ();
            z2 = (*riter)->getZ();

            r1 = (*iter)->getR();
            r2 = (*riter)->getR();

            dx = x2-x1;
            dy = y2-y1;
            dz = z2-z1;

            dist = sqrt( pow(dx, 2) + pow(dy, 2) + pow(dz, 2) );

            if (dist < r1 + r2 + prm->bridgeTol)
            {

                //cylinder radius
                double rBeadSmallest = (r1 <= r2? r1 : r2);
                double rBridge = prm->relativeBridgeRadius * rBeadSmallest;
                double r1Bridge = prm->relativeBridgeRadius * r1;
                double r2Bridge = prm->relativeBridgeRadius * r2;


                //Cylinder start point
                x3 = x1 + factor * r1 * dx/dist;
                y3 = y1 + factor * r1 * dy/dist;
                z3 = z1 + factor * r1 * dz/dist;

                x4 = x2 - factor * r2 * dx/dist;
                y4 = y2 - factor * r2 * dy/dist;
                z4 = z2 - factor * r2 * dz/dist;

                //cylinder length
                dx3 = dx * (dist - factor * (r1+r2))/dist;
                dy3 = dy * (dist - factor * (r1+r2))/dist;
                dz3 = dz * (dist - factor * (r1+r2))/dist;

                if (prm->geomInfile.empty())
                {
                    if (r1 == r2)
                        tagBridge = factory::addCylinder(x3, y3, z3, dx3, dy3, dz3, rBridge);
                    else
                        tagBridge = factory::addCone(x3, y3, z3, dx3, dy3, dz3, r1Bridge, r2Bridge);

                    this->dimTagsBridges.push_back({3, tagBridge }) ;
                    /* this->bridgeTagRadiusPairs.push_back({tagBridge, rBridge}); */
                }


                model::mesh::field::add("Frustum", ++count);
                model::mesh::field::setNumber(count, "V1_inner", prm->lc_bridge * rBeadSmallest/(prm->refBeadRadius) );
                model::mesh::field::setNumber(count, "V2_inner", prm->lc_bridge * rBeadSmallest/(prm->refBeadRadius) );
                /* model::mesh::field::setNumber(count, "V1_outer", prm->lc_beads * rBeadSmallest/(radius_max) ); */
                /* model::mesh::field::setNumber(count, "V2_outer", prm->lc_beads * rBeadSmallest/(radius_max) ); */
                model::mesh::field::setNumber(count, "V1_outer", prm->lc_out);
                model::mesh::field::setNumber(count, "V2_outer", prm->lc_out);
                model::mesh::field::setNumber(count, "R1_inner", 0);
                model::mesh::field::setNumber(count, "R2_inner", 0);
                model::mesh::field::setNumber(count, "R1_outer", prm->fieldThresholdMaxFactor * r1Bridge);
                model::mesh::field::setNumber(count, "R2_outer", prm->fieldThresholdMaxFactor * r2Bridge);

                model::mesh::field::setNumber(count, "X1", x3);
                model::mesh::field::setNumber(count, "Y1", y3);
                model::mesh::field::setNumber(count, "Z1", z3);
                model::mesh::field::setNumber(count, "X2", x4);
                model::mesh::field::setNumber(count, "Y2", y4);
                model::mesh::field::setNumber(count, "Z2", z4);
                vCount.push_back(count);




            }

        }
    }

    std::cout << "done!" << std::endl;

    std::cout << "Creating container... " << std::flush;
    if (prm->geomInfile.empty())
    {
        if (prm->containerShape == 0)
            dimTagsCyl.push_back( {3, factory::addCylinder(xCyl,yCyl, zCylBot, 0,0,zCylTop-zCylBot, rCyl) } );
        else if (prm->containerShape == 1)
            dimTagsCyl.push_back( {3, factory::addBox(pb->xMin - prm->rCylDelta,pb->yMin - prm->rCylDelta,zCylBot -prm->rCylDelta, pb->xMax - pb->xMin + prm->rCylDelta,pb->yMax - pb->yMin + prm->rCylDelta, zCylTop-zCylBot + prm->rCylDelta)});
        else
        {
            std::cout << "Invalid Container Shape. Use 0 for Cylinder, 1 for Box." << std::endl;
            exit(-1);
        }

    }
    std::cout << "done!" << std::endl;

    if (prm->meshSizeMethod == 1)
    {
        std::string backgroundField;
        if (prm->lc_beads > prm->lc_out)
            backgroundField="Max";
        else
            backgroundField="Min";

        std::cout << "Setting Background Field: " << backgroundField << std::endl;

        model::mesh::field::add(backgroundField, ++count);
        model::mesh::field::setNumbers(count, "FieldsList", vCount);

        model::mesh::field::setAsBackgroundMesh(count);
    }

}

void Model::mesh(std::string outfile, Parameters * prm)
{
    // Store dimTags and dimTagsMap
    std::vector<std::pair<int, int> > ov;
    std::vector<std::pair<int, int> > bv;
    std::vector<std::pair<int, int> > cv;
    std::vector<std::vector<std::pair<int, int> > > ovv;

    long bool_start = gmsh::logger::getWallTime();

    if (prm->geomInfile.empty())
    {
        /* Fuse beads together (with bridges or without) */
        if (prm->booleanOperation == 1)
        {
            if (dimTagsBridges.size() == 0) dimTagsBridges = {3, dimTagsBeads.back()};

            std::cout << "Fusing Beads and Bridges... " << std::flush;
            factory::fuse(dimTagsBeads, dimTagsBridges, ov, ovv );
            std::cout << "done!" << std::endl;
        }
        else if (prm->booleanOperation == 2)
        {
            if (dimTagsBridges.size() == 0)
            {
                std::cout << "Error: No Bridges to cut from beads" << std::endl;
                exit(-1);
            }

            std::cout << "Capping Beads... " << std::flush;
            factory::cut(dimTagsBeads, dimTagsBridges, ov, ovv );
            std::cout << "done!" << std::endl;

        }
        else
        {
            ov = dimTagsBeads;
        }

        // Fragment cylinder w.r.t. beads
        if (prm->fragment)
        {
            if (ov.size() == 0) ov = dimTagsBeads;

            std::cout << "Fragmenting Volumes... " << std::flush;
            factory::fragment(dimTagsCyl, ov, bv, ovv );
            /* dimTagsInterstitial.push_back(bv.back()); */
            std::cout << "done!" << std::endl;
        }

        // Synchronize gmsh model with geometry kernel.
        std::cout << "synchronizing... " << std::flush;
        factory::synchronize();
        std::cout << "done!" << std::endl;

        std::cout << std::endl;

        std::cout << "Number of Beads: "            << dimTagsBeads.size()   << std::endl;
        std::cout << "Number of Bridges: "          << dimTagsBridges.size() << std::endl;
        std::cout << "Number of internal volumes: " << ov.size()             << std::endl;

    }

    long bool_duration = gmsh::logger::getWallTime() - bool_start;

    gmsh::logger::write("Boolean time: " + std::to_string(bool_duration) + " s", "info");
    std::cout << std::endl;

    // ============================
    // Named Physical Groups

    // if not reading geometry
    if (prm->geomInfile.empty())
        createNamedGroups(bv, prm->containerShape);

    if (prm->meshSizeMethod == 0)
    {
        //Set mesh size globally
        std::cout << "Setting global mesh size... ";
        model::getEntities(cv, 0);
        model::mesh::setSize(cv, prm->lc_beads);
        std:: cout << "done!" << std::endl;
    }

    std::cout << std::endl;

    //if write geometry
    if(!prm->geomOutfile.empty())
        gmsh::write("geometries/" + prm->geomOutfile);

    //if read geometry
    if (!prm->geomInfile.empty())
    {
        std::cout << "Importing geometry: " << prm->geomInfile << "... " << std::flush;
        factory::importShapes("geometries/" + prm->geomInfile, bv);
        std::cout << "done!" << std::endl;

        // Synchronize gmsh model with geometry kernel.
        std::cout << "synchronizing... " << std::flush;
        factory::synchronize();
        std::cout << "done!" << std::endl;

        createNamedGroups(bv, prm->containerShape);
    }

    // Synchronize gmsh model with geometry kernel.
    std::cout << "synchronizing... " << std::flush;
    factory::synchronize();
    std::cout << "done!" << std::endl;

    if(!prm->dryRun)
    {
        for(int i=1; i<prm->MeshGenerate; i++)
        {
            model::mesh::generate(i);
            gmsh::write(prm->outpath + "/"+ std::to_string(i) + "D_" + outfile );

        }
        // Generate mesh
        model::mesh::generate(prm->MeshGenerate);

        // Output Full mesh
        gmsh::write(prm->outpath + "/"+ outfile);

        std::cout << std::endl;
        std::cout << "Calculating mesh volumes. Note: Multiply outputs by " << pow(prm->MeshScalingFactor, 3) << std::endl;
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

        if (prm->outputFragments == 1)
        {

            model::removePhysicalGroups();

            if (prm->NamedInterstitialVolume)
            {
                model::addPhysicalGroup(3,tVInt,15);
                model::setPhysicalName(3,15,"interstitialVolume");
            }

            if (prm->NamedOuterSurface)
            {
                model::addPhysicalGroup(2,tSInlet,11);
                model::setPhysicalName(2,11, "inlet");
                model::addPhysicalGroup(2,tSOutlet,12);
                model::setPhysicalName(2,12, "outlet");
                model::addPhysicalGroup(2,tSWall,13);
                model::setPhysicalName(2,13, "wall");
            }

            gmsh::write(prm->outpath + "/" + outfile + "_interstitial." + prm->fragmentFormat);

            model::removePhysicalGroups();

            if (prm->NamedBeadVolume)
            {
                model::addPhysicalGroup(3,tVBeads,16 );
                model::setPhysicalName(3,16,"beadVolume");
            }

            if (prm->NamedBeadSurface)
            {
                model::addPhysicalGroup(2,tSBeads,14);
                model::setPhysicalName(2,14, "beadSurface");
            }

            gmsh::write(prm->outpath + "/" + outfile + "_beads." + prm->fragmentFormat);

        }


    }


}

void Model::createNamedGroups(std::vector<std::pair<int,int>> bv, int containerShape)
{
    // bv is the vector storing dimtags after fragmentation.
    // the last entry should be the interstitial volume
    std::vector<std::pair<int,int>> cv;
    std::cout << "Listing Interstitial Volume... " << std::flush;
    tVInt.push_back(bv.back().second);
    bv.pop_back();
    std::cout << "done!" << std::endl;

    // Extract surfaces from interstitial volume
    std::cout << "Listing Outer Surfaces... " << std::flush;
    model::getBoundary({{3,tVInt[0]}}, cv, false, false, false);
    if(containerShape == 0)
    {
        tSWall = {cv[0].second};
        tSOutlet= {cv[1].second};
        tSInlet = {cv[2].second};
    }
    else if (containerShape == 1)
    {
        tSWall = {cv[0].second, cv[1].second, cv[3].second, cv[5].second};
        tSOutlet= {cv[2].second};
        tSInlet = {cv[4].second};
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
    model::getBoundary(bv, cv, false, false, false);
    for ( std::vector<std::pair< int , int>>::iterator it = cv.begin(); it != cv.end(); it++   )
    {
        tSBeads.push_back((*it).second);
    }
    std::cout << "done!" << std::endl;

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
