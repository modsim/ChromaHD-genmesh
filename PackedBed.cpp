/*
 * File: PackedBed.cpp
 * Created by: Rao, Jayghosh Subodh
 * Created on: Sat 06 Apr 2019 01:03:26 PM CEST
 *
 * This class is modeled after and copied from
 * BeadsInside from the old mesher.
 *
 * TODO: Save volumes differently. Surfaces (beads and int); volumes (beads and int)
 */

#include "PackedBed.h"
#include <iostream>
#include "Files.h"
#include <gmsh.h>
#include <iterator>

#include <float.h>

#define PI 3.1415926535897932

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

PackedBed::PackedBed(Parameters * prm)
{
    this->prm = prm;
    this->getBeads(this->prm->packfile);

    /* beads.clear(); */
    /* beads.push_back(new Bead(5, 5, 710.5, this->prm->rFactor * 1)); */
    /* beads.push_back(new Bead(6.5, 5,710.5, this->prm->rFactor * 0.50000)); */

    this->transformBeads();
    model::add("PackedBed");
}

PackedBed::~PackedBed()
{

}

void PackedBed::createGeometry()
{
    std::vector<std::pair<int, int> > ov;
    std::vector<std::pair<int, int> > cv;

    std::vector<double> vCount;

    double zCylBot = 0, zCylTop = 0;

    zCylBot     = ( this->prm->zBot - this->prm->inlet );
    zCylTop     = ( this->prm->zTop + this->prm->outlet );
    double xCyl = this->prm->xCyl;
    double yCyl = this->prm->yCyl;
    double rCyl = this->prm->rCyl;
    int count = 0;


    int tag, ctag;
    std::cout << "Creating beads... " << std::flush;
    for (std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++ )
    {
        double x = (*iter)->getX();
        double y = (*iter)->getY();
        double z = (*iter)->getZ();
        double r = (*iter)->getR();

        tag = factory::addSphere(x, y, z, r);
        /* ctag = factory::addPoint(x, y, z, this->prm->lc_beads, -1); */

        (*iter)->setTag(tag);
        /* (*iter)->setCTag(ctag); */

        dimTagsBeads.push_back({3, tag});
        /* tBeadCPs.push_back(ctag); */

        model::mesh::field::add("Ball", ++count);
        model::mesh::field::setNumber(count, "VIn", r/(radius_max) * this->prm->lc_beads);
        model::mesh::field::setNumber(count, "VOut", this->prm->lc_out);
        model::mesh::field::setNumber(count, "XCenter", x);
        model::mesh::field::setNumber(count, "YCenter", y);
        model::mesh::field::setNumber(count, "ZCenter", z);
        model::mesh::field::setNumber(count, "Radius", this->prm->fieldExtensionFactor * r);
        vCount.push_back(count);

    }
    std::cout << "done!" << std::endl;

    int tagBridge;

    std::vector<Bead *> remBeads = this->beads;
    std::cout << "Creating Bridges... " << std::flush;
    for(std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++)
    {
        remBeads.erase(remBeads.begin());
        for(std::vector<Bead *>::iterator riter = remBeads.begin(); riter != remBeads.end(); riter++)
        {
            //check pair with tol
            double dx, dy, dz;
            double x3, y3, z3;
            double x4, y4, z4;
            double dx3, dy3, dz3;

            double x1 = (*iter)->getX();
            double x2 = (*riter)->getX();

            double y1 = (*iter)->getY();
            double y2 = (*riter)->getY();

            double z1 = (*iter)->getZ();
            double z2 = (*riter)->getZ();

            double r1 = (*iter)->getR();
            double r2 = (*riter)->getR();

            dx = x2-x1;
            dy = y2-y1;
            dz = z2-z1;

            double dist = sqrt( pow(dx, 2) + pow(dy, 2) + pow(dz, 2) );

            if (dist < r1 + r2 + this->prm->bridgeTol)
            {

                //cylinder radius
                double rBeadSmallest = (r1 <= r2? r1 : r2);
                double rBridge = this->prm->relativeBridgeRadius * rBeadSmallest;
                double r1Bridge = this->prm->relativeBridgeRadius * r1;
                double r2Bridge = this->prm->relativeBridgeRadius * r2;

                double factor = this->prm->bridgeOffsetRatio;

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

                if (r1 == r2)
                    tagBridge = factory::addCylinder(x3, y3, z3, dx3, dy3, dz3, rBridge);
                else
                    tagBridge = factory::addCone(x3, y3, z3, dx3, dy3, dz3, r1Bridge, r2Bridge);

                this->dimTagsBridges.push_back({3, tagBridge }) ;
                /* this->bridgeTagRadiusPairs.push_back({tagBridge, rBridge}); */


/*                     model::mesh::field::add("Cylinder", ++count); */
/*                     model::mesh::field::setNumber(count, "VIn", this->prm->lc_beads * rBeadSmallest/(radius_max) ); */
/*                     model::mesh::field::setNumber(count, "VOut", this->prm->lc_out); */
/*                     model::mesh::field::setNumber(count, "XCenter", x3); */
/*                     model::mesh::field::setNumber(count, "YCenter", y3); */
/*                     model::mesh::field::setNumber(count, "ZCenter", z3); */
/*                     model::mesh::field::setNumber(count, "XAxis", dx3); */
/*                     model::mesh::field::setNumber(count, "YAxis", dy3); */
/*                     model::mesh::field::setNumber(count, "ZAxis", dz3); */
/*                     model::mesh::field::setNumber(count, "Radius", this->prm->fieldExtensionFactor * rBridge); */
/*                     vCount.push_back(count); */

                model::mesh::field::add("Frustum", ++count);
                model::mesh::field::setNumber(count, "V1_inner", this->prm->lc_beads * rBeadSmallest/(radius_max) );
                model::mesh::field::setNumber(count, "V2_inner", this->prm->lc_beads * rBeadSmallest/(radius_max) );
                /* model::mesh::field::setNumber(count, "V1_outer", this->prm->lc_beads * rBeadSmallest/(radius_max) ); */
                /* model::mesh::field::setNumber(count, "V2_outer", this->prm->lc_beads * rBeadSmallest/(radius_max) ); */
                model::mesh::field::setNumber(count, "V1_outer", this->prm->lc_out);
                model::mesh::field::setNumber(count, "V2_outer", this->prm->lc_out);
                model::mesh::field::setNumber(count, "R1_inner", 0);
                model::mesh::field::setNumber(count, "R2_inner", 0);
                model::mesh::field::setNumber(count, "R1_outer", this->prm->fieldExtensionFactor * r1Bridge);
                model::mesh::field::setNumber(count, "R2_outer", this->prm->fieldExtensionFactor * r2Bridge);

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

    std::cout << "Creating cylinder... " << std::flush;
    dimTagsCyl.push_back( {3, factory::addCylinder(xCyl,yCyl, zCylBot, 0,0,zCylTop-zCylBot, rCyl) } );
    std::cout << "done!" << std::endl;

    if (this->prm->meshSizeMethod == 1)
    {
        std::string backgroundField;
        if (this->prm->lc_beads > this->prm->lc_out)
            backgroundField="Max";
        else
            backgroundField="Min";

        std::cout << "Setting Background Field: " << backgroundField << std::endl;

        model::mesh::field::add(backgroundField, ++count);
        model::mesh::field::setNumbers(count, "FieldsList", vCount);

        model::mesh::field::setAsBackgroundMesh(count);
    }

}

void PackedBed::mesh(std::string outfile)
{
    // Store dimTags and dimTagsMap
    std::vector<std::pair<int, int> > ov;
    std::vector<std::pair<int, int> > bv;
    std::vector<std::pair<int, int> > cv;
    std::vector<std::vector<std::pair<int, int> > > ovv;


    /* // Synchronize gmsh model with geometry kernel. */
    /* std::cout << "synchronizing... " << std::flush; */
    /* factory::synchronize(); */
    /* std::cout << "done!" << std::endl; */

    /* // Set scaled mesh size on beads */
    /* std::cout << "Setting bead sizes..." << std::endl; */
    /* std::vector<Bead *>::iterator biter = beads.begin(); */
    /* for (std::vector<std::pair <int, int >>::iterator it = this->dimTagsBeads.begin(); */
    /*         it != this->dimTagsBeads.end(); it++) */
    /* { */
    /*     model::getBoundary({(*it)}, cv, false, false, true ); */
    /*     model::mesh::setSize(cv, this->prm->lc * (*biter)->getR()/radius_max); */
    /*     /1* model::mesh::setSize({(*it)}, this->prm->lc * (*biter)->getR()/radius_avg); *1/ */
    /*     biter++; */
    /* } */

    /* std::cout << "Setting bridge sizes..." << std::endl; */
    /* /1* for (std::vector<std::pair <int, int >>::iterator it = this->dimTagsBridges.begin(); *1/ */
    /* /1*         it != this->dimTagsBridges.end(); it++) *1/ */
    /* for (std::vector<std::pair <int, double>>::iterator it = this->bridgeTagRadiusPairs.begin(); */
    /*         it != this->bridgeTagRadiusPairs.end(); it++) */
    /* { */
    /*     model::getBoundary({{3,(*it).first}}, cv, false, false, true ); */
    /*     model::mesh::setSize(cv, this->prm->lc * (*it).second/(this->prm->relativeBridgeRadius * radius_max) ); */
    /*     model::mesh::setSize(cv, this->prm->lc * (*it).second/(radius_max) ); */
    /*     /1* model::mesh::setSize({{3,(*it).first}}, this->prm->lc * (*it).second/(this->prm->relativeBridgeRadius * radius_avg) ); *1/ */
    /* } */

    long bool_start = gmsh::logger::time();

    /* Fuse beads together (with bridges or without) */
    if (this->prm->booleanOperation == 1)
    {
        if (dimTagsBridges.size() == 0) dimTagsBridges = {3, dimTagsBeads.back()};

        std::cout << "Fusing Beads and Bridges... " << std::flush;
        factory::fuse(dimTagsBeads, dimTagsBridges, ov, ovv );
        std::cout << "done!" << std::endl;
    }
    else if (this->prm->booleanOperation == 2)
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
    if (this->prm->fragment)
    {
        if (ov.size() == 0) ov = dimTagsBeads;

        std::cout << "Fragmenting Volumes... " << std::flush;
        factory::fragment(dimTagsCyl, ov, bv, ovv );
        dimTagsInterstitial.push_back(bv.back());
        std::cout << "done!" << std::endl;
    }

    long bool_duration = gmsh::logger::time() - bool_start;

    // Synchronize gmsh model with geometry kernel.
    std::cout << "synchronizing... " << std::flush;
    factory::synchronize();
    std::cout << "done!" << std::endl;

    std::cout << std::endl;

    std::cout << "Number of Beads: "            << dimTagsBeads.size()   << std::endl;
    std::cout << "Number of Bridges: "          << dimTagsBridges.size() << std::endl;
    std::cout << "Number of internal volumes: " << ov.size()             << std::endl;

    gmsh::logger::write("Boolean time: " + std::to_string(bool_duration) + " s", "info");
    std::cout << std::endl;

    /* std:: cout << "Embedding control points within spheres..." << std::flush; */
    /* model::mesh::embed(0,tBeadCPs, 3, ov[0].second); */

    /* std:: cout << "Embedding control points within spheres..." << std::flush; */
    /* //Embed points in beads to control element size within them */
    /* if (this->prm->booleanOperation == 1) */
    /* { */
    /*     // if fuse: embed all bead centers into the bead fragment */
    /*     /1* std::cout << bv.size() << std::endl; *1/ */
    /*     if (bv.size() == 0) bv = ov; */
    /*     std::cout << ov.size() << std::endl; */
    /*     std::cout << bv.size() << std::endl; */


    /*     for (std::vector<std::pair<int,int>>::iterator iter = bv.begin(); iter != bv.end(); iter++ ) */
    /*     { */
    /*         std::cout << "inside loop.." << std::endl; */
    /*         model::mesh::embed(0,tBeadCPs, 3, (*iter).second); */
    /*     } */
    /* } */
    /* else */
    /* { */
    /*     //if cap or reduce, embed into beads themeselves. */
    /*     for (std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++ ) */
    /*     { */
    /*         model::mesh::embed(0,{(*iter)->getCTag()}, 3, (*iter)->getTag()); */
    /*     } */
    /* } */
    /* std:: cout << "done!" << std::endl; */

    /* /1* // Set mesh size on cylinder surface *1/ */
    /* model::getBoundary(dimTagsInterstitial, cv, false, false, true ); */
    /* cv.erase(cv.begin()+3, cv.end()); */
    /* model::mesh::setSize(cv, this->prm->lc); */


    // ============================
    // Named Physical Groups
    if (this->prm->NamedInterstitialVolume)
    {
        std::cout << "Naming Interstitial Volumes... ";
        model::addPhysicalGroup(3,{bv.back().second},5);
        model::setPhysicalName(3,5,"interstitialVolume");
        std::cout << "done!" << std::endl;
    }

    if (this->prm->NamedOuterSurface)
    {
        std::cout << "Naming Outer Surfaces... ";
        model::getBoundary({{3, bv.back().second}}, cv, false, false,false);
        model::addPhysicalGroup(2,{cv[0].second},3);
        model::addPhysicalGroup(2,{cv[1].second},2);
        model::addPhysicalGroup(2,{cv[2].second},1);
        model::setPhysicalName(2,1, "inlet");
        model::setPhysicalName(2,2, "outlet");
        model::setPhysicalName(2,3, "wall");
        std::cout << "done!" << std::endl;
    }


    bv.pop_back();

    if (this->prm->NamedBeadVolume)
    {
        tBeads.clear();
        std::cout << "Naming Bead Volumes... ";
        for ( std::vector<std::pair< int , int>>::iterator it = bv.begin(); it != bv.end(); it++   )
        {
            tBeads.push_back((*it).second);
        }
        model::addPhysicalGroup(3,tBeads,6 );
        model::setPhysicalName(3,6,"beadVolume");
        std::cout << "done!" << std::endl;

    }

    if (this->prm->NamedBeadSurface)
    {
        std::cout << "Naming Bead Surfaces... ";
        model::getBoundary(bv, cv, false, false,false);
        tBeads.clear();
        for ( std::vector<std::pair< int , int>>::iterator it = cv.begin(); it != cv.end(); it++   )
        {
            tBeads.push_back((*it).second);
        }
        model::addPhysicalGroup(2,tBeads,4);
        model::setPhysicalName(2,4, "beadSurface");
        std::cout << "done!" << std::endl;
    }

    if (this->prm->meshSizeMethod == 0)
    {
        //Set mesh size globally
        std::cout << "Setting global mesh size... ";
        model::getEntities(cv, 0);
        model::mesh::setSize(cv, this->prm->lc_beads);
        std:: cout << "done!" << std::endl;
    }

    std::cout << std::endl;


    // Set mesh size for beads
    /* model::getBoundary(bv, cv, false, false, true); */
    /* model::mesh::setSize(cv, this->prm->lc_beads); */

    /* // Set mesh size for interstitial */
    /* std:: cout << "Setting mesh size for surfaces..."; */
    /* model::getBoundary(dimTagsInterstitial, cv, false, false, true); */
    /* model::mesh::setSize(cv, this->prm->lc); */
    /* std:: cout << "done!" << std::endl; */

    //set mesh size on bead surface
    // maybe use the inside of the interstitial fragment instead?
    /* model::getBoundary(dimTagsInterstitial, cv, false, false, false); */
    /* model::getBoundary(bv, cv, false, false, false); */
    /* model::getBoundary(cv, ov, false, false, true); */
    /* model::mesh::setSize(ov, this->prm->lc); */

    /* gmsh::write(this->prm->outpath + "geometry.brep" ); */


    if(!this->prm->dryRun)
    {
        // Generate mesh
        model::mesh::generate(this->prm->MeshGenerate);

        // Output Full mesh
        gmsh::write(this->prm->outpath + outfile);

        model::removePhysicalGroups();

        if (this->prm->NamedInterstitialVolume)
        {
            model::addPhysicalGroup(3,{dimTagsInterstitial.back().second},15);
            model::setPhysicalName(3,15,"interstitialVolume");
        }

        if (this->prm->NamedOuterSurface)
        {
            model::getBoundary(dimTagsInterstitial, cv, false, false,false);
            model::addPhysicalGroup(2,{cv[0].second},13);
            model::addPhysicalGroup(2,{cv[1].second},12);
            model::addPhysicalGroup(2,{cv[2].second},11);
            model::setPhysicalName(2,11, "inlet");
            model::setPhysicalName(2,12, "outlet");
            model::setPhysicalName(2,13, "wall");
        }

        gmsh::write(this->prm->outpath + outfile + "_interstitial.vtk");

        model::removePhysicalGroups();

        if (this->prm->NamedBeadVolume)
        {
            tBeads.clear();
            for ( std::vector<std::pair< int , int>>::iterator it = bv.begin(); it != bv.end(); it++   )
            {
                tBeads.push_back((*it).second);
            }
            model::addPhysicalGroup(3,tBeads,16 );
            model::setPhysicalName(3,16,"beadVolume");

        }

        /* if (this->prm->NamedBeadSurface) */
        /* { */
        /*     model::getBoundary(bv, cv, false, false,false); */
        /*     tBeads.clear(); */
        /*     for ( std::vector<std::pair< int , int>>::iterator it = cv.begin(); it != cv.end(); it++   ) */
        /*     { */
        /*         tBeads.push_back((*it).second); */
        /*     } */
        /*     model::addPhysicalGroup(2,tBeads,14); */
        /*     model::setPhysicalName(2,14, "beadSurface"); */
        /* } */

        gmsh::write(this->prm->outpath + outfile + "_beads.msh2");

        gmsh::plugin::run("NewView");
        gmsh::plugin::setString("ModifyComponents", "Expression0", "1");
        gmsh::plugin::run("ModifyComponents");
        gmsh::plugin::setNumber("Integrate", "Dimension", 3);
        gmsh::plugin::run("Integrate");

    }




}

void PackedBed::printPacking()
{
    std::cout << "Printing packing geometry..." << std::endl;
    for (std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++ )
    {
        std::cout << std::setw(20) << ( *iter )->getX()
            << std::setw(20) << ( *iter )->getY()
            << std::setw(20) << ( *iter )->getZ()
            << std::setw(20) << ( *iter )->getR()
            << std::setw(20) << std::endl;
    }

}

void PackedBed::getBeads(std::string packingFilename)
{
    // open the file:
    std::ifstream file(packingFilename.c_str(), std::ios::binary);

    // Stop eating new lines in binary mode!!!
    file.unsetf(std::ios::skipws);

    file.seekg(0, std::ios::end);
    const size_t num_elements = file.tellg() / sizeof(float);
    file.seekg(0, std::ios::beg);

    std::cout << "reading packing ... " << std::flush;
    std::vector<float> data(num_elements);
    file.read(reinterpret_cast<char*>(&data[0]), num_elements*sizeof(float));

    if (isBigEndian())
        swapbytes(reinterpret_cast<char*>(&data[0]), data.size(), sizeof (float));

    std::cout << "done!" << std::endl << std::endl;

    std::cout << "Note: Bead data stored in the vector is already modified by rFactor." << std::endl;
    std::cout << "Selecting beads in range..." << std::endl;

    double cz1 = this->prm->zBot;
    double cz2 = this->prm->zTop;

    nBeadsMax = data.size()/4;

    //NOTE: Bead data stored in the vector is already modified by rFactor
    for(size_t i = 0; i < nBeadsMax; ++i)
    {
        double x = data[i * 4 ];
        double y = data[i * 4 + 1];
        double z = data[i * 4 + 2];
        double r = data[i * 4 + 3] * 0.5;

        beads.push_back(new Bead(x, y, z, this->prm->rFactor * r));
    }
    std::cout << "done!" << std::endl;

    if (beads.size() == 0) {
        std::cerr << "ERROR: No beads found!" << std::endl;
        exit(1);
    }

    // sort using a lambda expression
    std::cout << "Sorting beads according to z-value..." << std::flush;
    std::sort(beads.begin(), beads.end(), [](const Bead* b1, const Bead* b2) {
            return b1->getZ() < b2->getZ();
    });
    std::cout << "done!" << std::endl;

    //Only store nBeads
    double nBeads = this->prm->nBeads > nBeadsMax? nBeadsMax: this->prm->nBeads;
    beads.erase(beads.begin() + nBeads, beads.end());

    std::cout << this->beads.size() << "/" << nBeadsMax << " beads in the selected range." << std::endl;

}

void PackedBed::transformBeads()
{
    //Calculate bounding box
    double xCyl = 0, yCyl = 0, rCyl = 0;
    double xMax = -DBL_MAX, yMax = -DBL_MAX, zMax = -DBL_MAX;
    double xMin = DBL_MAX, yMin = DBL_MAX, zMin = DBL_MAX;
    double x=0, y=0, r=0, z=0;
    double zBot=0, zTop=0;
    radius_avg=0;

    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
    {
        x = (*it)->getX();
        y = (*it)->getY();
        r = (*it)->getR();

        radius_avg +=  r;

        if ( (x + r) > xMax ) xMax = x + r;
        if ( (x - r) < xMin ) xMin = x - r;
        if ( (y + r) > yMax ) yMax = y + r;
        if ( (y - r) < yMin ) yMin = y - r;

        if (r > radius_max) radius_max = r;
        if (r < radius_min) radius_min = r;

    }
    radius_avg /= beads.size();

    xCyl = (xMax + xMin) / 2;
    yCyl = (yMax + yMin) / 2;
    rCyl = std::max( ((xMax - xMin) / 2), ((yMax - yMin)/2) ) + this->prm->rCylDelta;
    /* rCyl = std::max(std::max(xMax, std::abs(xMin)), std::max(yMax, std::abs(yMin) ) ) + this->prm->rCylDelta; */

    //Adjust zBot and zTop
    zBot=beads.front()->getZ();
    zTop=beads.back()->getZ();

    std::cout << std::endl;
    std::cout << "== Before Transform ==" << std::endl;
    std::cout << "xCyl: " << xCyl << std::endl;
    std::cout << "yCyl: " << yCyl << std::endl;
    std::cout << "rCyl: " << rCyl << std::endl;
    std::cout << "xMax: " << xMax << std::endl;
    std::cout << "xMin: " << xMin << std::endl;
    std::cout << "yMax: " << yMax << std::endl;
    std::cout << "yMin: " << yMin << std::endl;
    std::cout << "zBot: " << zBot << std::endl;
    std::cout << "zTop: " << zTop << std::endl;
    std::cout << std::endl;

    double offsetx = -xCyl;
    double offsety = -yCyl;
    double offsetz = -zBot;

    std::cout << "Translating beads by (" << offsetx << ", "<< offsety << ", " << offsetz << ")... "<<  std::flush;
    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
        (*it)->translate(offsetx, offsety, offsetz);
    std::cout << "done!" << std::endl;

    std::cout << "Scaling beads by " << this->prm->preScalingFactor << "... " << std::flush;
    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
        (*it)->scale(this->prm->preScalingFactor);
    std::cout << "done!" << std::endl << std::endl;

    //Adjust zBot and zTop
    zBot=beads.front()->getZ();
    zTop=beads.back()->getZ();

    xCyl = 0, yCyl = 0, rCyl = 0;
    xMax = -DBL_MAX, yMax = -DBL_MAX, zMax = -DBL_MAX;
    xMin = DBL_MAX, yMin = DBL_MAX, zMin = DBL_MAX;
    x=0, y=0, r=0, z=0;
    radius_avg=0;

    double vol_real_beads=0;
    double vol_geom_beads=0;
    double vol_real_int=0;
    double vol_mesh_int=0;
    double vol_cylinder=0;

    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
    {
        x = (*it)->getX();
        y = (*it)->getY();
        z = (*it)->getZ();
        r = (*it)->getR();

        radius_avg +=  r;
        vol_real_beads += PI * 4/3 * pow(r * this->prm->MeshScalingFactor * 1/(this->prm->rFactor), 3);
        vol_geom_beads += PI * 4/3 * pow(r * this->prm->MeshScalingFactor, 3);

        if ( (x + r) > xMax ) xMax = x + r;
        if ( (x - r) < xMin ) xMin = x - r;
        if ( (y + r) > yMax ) yMax = y + r;
        if ( (y - r) < yMin ) yMin = y - r;
        if ( (z + r) > zMax ) zMax = z + r;
        if ( (z - r) < zMin ) zMin = z - r;

        if (r > radius_max) radius_max = r;
        if (r < radius_min) radius_min = r;

    }
    radius_avg /= beads.size();

    xCyl = (xMax + xMin) / 2;
    yCyl = (yMax + yMin) / 2;
    rCyl = std::max( ((xMax - xMin) / 2), ((yMax - yMin)/2) ) + this->prm->rCylDelta;
    /* rCyl = std::max(std::max(xMax, std::abs(xMin)), std::max(yMax, std::abs(yMin) ) ) + this->prm->rCylDelta; */

    vol_cylinder = PI * pow(rCyl * this->prm->MeshScalingFactor, 2) * (zTop - zBot + this->prm->inlet + this->prm->outlet) * this->prm->MeshScalingFactor;
    vol_real_int = vol_cylinder - vol_real_beads;

    this->prm->xCyl = xCyl;
    this->prm->yCyl = yCyl;
    this->prm->rCyl = rCyl;
    this->prm->zBot = zBot;
    this->prm->zTop = zTop;

    std::cout << "== After Transform ==" << std::endl;
    std::cout << "xCyl: " << xCyl << std::endl;
    std::cout << "yCyl: " << yCyl << std::endl;
    std::cout << "rCyl: " << rCyl << std::endl;
    std::cout << "xMax: " << xMax << std::endl;
    std::cout << "xMin: " << xMin << std::endl;
    std::cout << "yMax: " << yMax << std::endl;
    std::cout << "yMin: " << yMin << std::endl;
    std::cout << "zBot: " << zBot << std::endl;
    std::cout << "zTop: " << zTop << std::endl << std::endl;

    std::cout << "average bead radius: " << std::scientific << std::setprecision(10) << radius_avg << std::endl;
    std::cout << "maximum bead radius: " << std::scientific << std::setprecision(10) << radius_max << std::endl;
    std::cout << "minimum bead radius: " << std::scientific << std::setprecision(10) << radius_min << std::endl << std::endl;

    std::cout << "=== After Mesh Scaling ===" << std::endl;
    std::cout << "Cylinder Volume: "      << vol_cylinder   << std::endl;
    std::cout << "Real Bead Volume: "     << vol_real_beads << std::endl;
    /* std::cout << "Real Int Volume: "      << vol_real_int   << std::endl; */
    std::cout << "Modified Bead Volume: " << vol_geom_beads << std::endl<< std::endl;

    std::cout << "Bed Length: " << (zMax-zMin)*this->prm->MeshScalingFactor << std::endl;
    std::cout << "Column Length: " << (zTop - zBot + this->prm->inlet + this->prm->outlet) * this->prm->MeshScalingFactor << std::endl<<std::endl;

    std::cout << "Note: Modified Bead Volume is currently only accurate for reduced beads." << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Bead Geometry Volume Error: " << (vol_real_beads - vol_geom_beads)/vol_real_beads*100 << "%" << std::endl << std::endl;

}
