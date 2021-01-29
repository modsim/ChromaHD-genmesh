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

void printVecPair(std::vector<std::pair<int, int>> vp)
{
    for (std::vector<std::pair<int,int>>::iterator iter = vp.begin(); iter != vp.end(); iter++)
    {
        std::cout << (*iter).first << "    " << (*iter).second << std::endl;
    }
}

/* Function to create bead geometry.
 * xoff, yoff are for periodic cases, to stack the bead geometries in the lateral dimensions
 */
void createBeadGeometry(std::vector<Bead *> beads, Parameters* prm, std::vector<std::pair<int, int>> &dimTagsBeads, std::vector<int> &tBeadCPs, int &count, std::vector<double> &vCount, double xoff, double yoff)
{
    std::cout << "Creating beads... " << std::flush;
    int tag, ctag;

    for (std::vector<Bead *>::iterator iter = beads.begin(); iter != beads.end(); iter++ )
    {
        double x = (*iter)->getX() + xoff;
        double y = (*iter)->getY() + yoff;
        double z = (*iter)->getZ();
        double r = (*iter)->getR();

        if (prm->geomInfile.empty())
        {
            tag = factory::addSphere(x, y, z, r);
            /* std::cout << tag << std::endl; */
            ctag = factory::addPoint(x, y, z, prm->lc_beads, -1);

            (*iter)->setTag(tag);
            (*iter)->setCTag(ctag);

            dimTagsBeads.push_back({3, tag});
            tBeadCPs.push_back(ctag);
        }

        /*
         * Distance/Threshold allows creating gradient mesh sizes within the beads
         */
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

}

/*
 * Create container geometry
 */
void createContainerGeometry(PackedBed * pb, Parameters * prm, double xCyl, double yCyl, double rCyl, double zCylTop, double zCylBot, std::vector<std::pair<int, int>> &dimTagsCyl)
{
    std::cout << "Creating container... " << std::flush;
    if (prm->geomInfile.empty())
    {
        if (prm->autoContainment == 1)
        {

            if (prm->containerShape == 0)
            {
                std::cout << "Creating cylindrical container automatically." << std::endl;
                dimTagsCyl.push_back( {3, factory::addCylinder(xCyl,yCyl, zCylBot, 0,0,zCylTop-zCylBot, rCyl) } );

            }
            else if (prm->containerShape == 1)
            {
                std::cout << "Creating rectangular container automatically." << std::endl;
                dimTagsCyl.push_back( {3, factory::addBox(pb->xMin - prm->rCylDelta,pb->yMin - prm->rCylDelta,zCylBot, pb->xMax - pb->xMin + 2 * prm->rCylDelta,pb->yMax - pb->yMin + 2 * prm->rCylDelta, zCylTop-zCylBot)});
            }
            else
            {
                std::cout << "Skipping container creation." << std::endl;
            }

        }
        else if (prm->autoContainment == 0)
        {
            if (prm->containerShape == 0)
            {
                std::cout << "Creating cylindrical container manually." << std::endl;
                //TODO: Ensure that xCyl yCyl rCyl are immutable in prm
                dimTagsCyl.push_back( {3, factory::addCylinder(prm->xCyl,prm->yCyl, zCylBot, 0,0,zCylTop-zCylBot, prm->rCyl) } );

            }
            else if (prm->containerShape == 1)
            {
                std::cout << "Creating rectangular container manually." << std::endl;
                dimTagsCyl.push_back( {3, factory::addBox(prm->x0, prm->y0, prm->z0, prm->dx, prm->dy, prm->dz)});
            }
            else
            {
                std::cout << "Skipping container creation manually." << std::endl;
            }

        }
        else
        {
            std::cout << "Invalid autoContainment value!" << std::endl;
            exit(-1);
        }

    }

}

void Model::createGeometry(PackedBed * pb, Parameters * prm)
{
    if(prm->dryRun)
        return;

    std::vector<Bead *> beads = pb->beads;

    std::vector<std::pair<int, int> > ov;
    std::vector<std::pair<int, int> > cv;

    std::vector<double> vCount;

    zCylBot = ( prm->zBot - prm->inlet  );
    zCylTop = ( prm->zTop + prm->outlet );
    xMax = pb->xMax;
    yMax = pb->yMax;
    xMin = pb->xMin;
    yMin = pb->yMin;

    double xCyl    = prm->xCyl;
    double yCyl    = prm->yCyl;
    double rCyl    = prm->rCyl;

    int count = 0;

    createBeadGeometry(beads, prm, dimTagsBeads, tBeadCPs, count, vCount, 0, 0);

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

    createContainerGeometry(pb, prm, xCyl, yCyl, rCyl, zCylTop, zCylBot, dimTagsCyl);

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
    if(prm->dryRun)
        return;

    // Store dimTags and dimTagsMap
    std::vector<std::pair<int, int> > ov;
    std::vector<std::pair<int, int> > bv;
    std::vector<std::pair<int, int> > cv;
    std::vector<std::vector<std::pair<int, int> > > ovv;

    /* cv = dimTagsBeads; */

    // TODO: Use ov, bv, cv consistently.
    // TODO: Rename ov, bv, cv to be more descriptive

    std::cout << "Intersecting Volumes... " << std::flush;
    factory::intersect(dimTagsBeads, dimTagsCyl, cv, ovv, -1, true, false);
    std::cout << "done!" << std::endl;

    long bool_start = gmsh::logger::getWallTime();

    if (prm->geomInfile.empty())
    {
        /* Fuse beads together (with bridges or without) */
        if (prm->booleanOperation == 1)
        {
            if (dimTagsBridges.size() == 0) dimTagsBridges = {3, dimTagsBeads.back()};

            std::cout << "Fusing Beads and Bridges... " << std::flush;
            /* factory::fuse(dimTagsBeads, dimTagsBridges, ov, ovv ); */
            factory::fuse(cv, dimTagsBridges, ov, ovv );
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
            /* factory::cut(dimTagsBeads, dimTagsBridges, ov, ovv ); */
            factory::cut(cv, dimTagsBridges, ov, ovv );
            std::cout << "done!" << std::endl;

        }
        else
        {
            /* ov = dimTagsBeads; */
            ov = cv;
        }

        // Fragment cylinder w.r.t. beads
        if (prm->fragment)
        {
            /* if (ov.size() == 0) ov = dimTagsBeads; */
            if (ov.size() == 0) ov = cv;

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

        /* std::cout << "Number of Beads: "            << dimTagsBeads.size()   << std::endl; */
        std::cout << "Number of Beads: "            << cv.size()   << std::endl;
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
    {
        if (prm->fragment == 1)
            createNamedGroups(bv, prm->containerShape);

        if (prm->periodic == 1)
            setupPeriodicSurfaces(prm);
    }
    else
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
            model::addPhysicalGroup(2, tSInlet  , 11      );
            model::addPhysicalGroup(2, tSOutlet , 12      );
            model::addPhysicalGroup(2, tSWall   , 13      );
            model::setPhysicalName (2, 11       , "inlet" );
            model::setPhysicalName (2, 12       , "outlet");
            model::setPhysicalName (2, 13       , "wall"  );
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

void Model::createNamedGroups(std::vector<std::pair<int,int>> bv, int containerShape)
{
    if (containerShape == -1)
    {
        std::cout << ">>> Not creating Named Groups!!" << std::endl;
        return;
    }

    // bv is the vector storing dimtags after fragmentation.
    // the last entry should be the interstitial volume
    std::vector<std::pair<int,int>> cv;
    std::vector<std::pair<int,int>> beadSurfaceDimTags;

    std::cout << "Listing Interstitial Volume... " << std::flush;
    tVInt.push_back(bv.back().second);
    bv.pop_back();
    std::cout << "done!" << std::endl;

    // Extract surfaces from interstitial volume
    std::cout << "Listing Outer Surfaces... " << std::flush;
    model::getBoundary({{3,tVInt[0]}}, cv, false, false, false);
    model::getBoundary(bv, beadSurfaceDimTags, false, false, false);
    if(containerShape == 0)
    {
        tSWall = {cv[0].second};
        tSOutlet= {cv[1].second};
        tSInlet = {cv[2].second};

        // NOTE: What would be the normal for a cylinder??
    }
    else if (containerShape == 1)
    {

        std::vector<double> points, parametricCoord, curvatures, normals;
        std::vector<std::pair <int, int>> dimTagsBound;

        // Interstitial surfaces
        for (std::vector<std::pair<int, int>>::iterator iter = cv.begin(); iter != cv.end(); iter++)
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

        std::cout << "SIZE INTER: "<< tXLeftWallInt.size() << " | "<< tXRightWallInt.size() << std::endl;

        // Bead Surfaces
        for (std::vector<std::pair<int, int>>::iterator iter = beadSurfaceDimTags.begin(); iter != beadSurfaceDimTags.end(); iter++)
        {
            /*
             * Get the curvatures of the surfaces,
             * Select planar surfaces
             * Get normals of planes to match surface position to tag
             */
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

        /* std::cout << "SIZE BOTH X : "<< tXLeftWallInt.size() << " | "<< tXRightWallInt.size() << std::endl; */
        /* std::cout << "SIZE BOTH Y : "<< tYLeftWallInt.size() << " | "<< tYRightWallInt.size() << std::endl; */


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

    //TODO: How would this work for periodic meshes?
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

void Model::setupPeriodicSurfaces(Parameters * prm)
{
    std::cout << "Setting up periodic surfaces..." << std::flush;

    if (tXLeftWallInt.size() != tXRightWallInt.size())
    {
        std::cout << "Mismatch in number of surfaces in X planes." << std::endl;
        exit(-1);
    }
    if (tYLeftWallInt.size() != tYRightWallInt.size())
    {
        std::cout << "Mismatch in number of surfaces in Y planes." << std::endl;
        exit(-1);
    }

    double dx, dy, dz;

    if (prm->autoContainment == 1)
    {
        dx = xMax - xMin + 2 * prm->rCylDelta;
        dy = yMax - yMin + 2 * prm->rCylDelta;
        dz = zCylTop - zCylBot;
    }
    else if (prm->autoContainment == 0)
    {
        dx = prm->dx;
        dy = prm->dy;
        dz = prm->dz;
    }

    std::vector<double> affineTranslationX = {1, 0, 0, dx, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    std::vector<double> affineTranslationY = {1, 0, 0, 0, 0, 1, 0, dy, 0, 0, 1, 0, 0, 0, 0, 1};
    std::vector<double> affineTranslationZ = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, dz, 0, 0, 0, 1};

    std::cout << tXLeftWallInt.size() << "   " << tXRightWallInt.size() << std::endl;
    std::cout << tYLeftWallInt.size() << "   " << tYRightWallInt.size() << std::endl;

    std::cout << tXLeftWallBead.size() << "   " << tXRightWallBead.size() << std::endl;
    std::cout << tYLeftWallBead.size() << "   " << tYRightWallBead.size() << std::endl;

    /* std::cout << tZLeftWallInt.size() << "   " << tZRightWallInt.size() << std::endl; */

    std::cout << tXLeftWallInt[0] << " | " << tXRightWallInt[0] << std::endl;
    std::cout << tYLeftWallInt[0] << " | " << tYRightWallInt[0] << std::endl;

    std::cout << tXLeftWallBead[0] << " | " << tXRightWallBead[0] << std::endl;
    std::cout << tYLeftWallBead[0] << " | " << tYRightWallBead[0] << std::endl;

    model::mesh::setPeriodic(2, tXRightWallInt, tXLeftWallInt, affineTranslationX);
    model::mesh::setPeriodic(2, tYRightWallInt, tYLeftWallInt, affineTranslationY);
    /* model::mesh::setPeriodic(2, tZLeftWallInt, tZRightWallInt, affineTranslationZ); */

    std::cout << "done!" << std::endl;

}

