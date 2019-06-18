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

#define PI 3.1415926535897932

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

PackedBed::PackedBed(Parameters * prm)
{
    this->prm = prm;
    this->readFile(this->prm->packfile);

    /* beads.push_back(new Bead(0, 0, 0, this->prm->rFactor * 1)); */
    /* beads.push_back(new Bead(0, 0, 2, this->prm->rFactor * 1)); */

    model::add("PackedBed");
}

PackedBed::~PackedBed()
{

}

void PackedBed::createGeometry()
{

    double zCylBot = 0, zCylTop = 0;

    double psf  = this->prm->preScalingFactor;
    zCylBot     = psf * ( this->prm->zBot - this->prm->inlet );
    zCylTop     = psf * ( this->prm->zTop + this->prm->outlet );
    double xCyl = psf * this->prm->xCyl;
    double yCyl = psf * this->prm->yCyl;
    double rCyl = psf * this->prm->rCyl;

    std::cout << "Creating cylinder... " << std::flush;
    /* dimTagsCyl.push_back( {3, factory::addCylinder(this->prm->xCyl,this->prm->yCyl, zCylBot, 0,0,zCylTop-zCylBot, this->prm->rCyl) } ); */
    dimTagsCyl.push_back( {3, factory::addCylinder(xCyl,yCyl, zCylBot, 0,0,zCylTop-zCylBot, rCyl) } );
    std::cout << "done!" << std::endl;

    std::vector<std::pair<int, int> > ov;

    int tag, ctag;
    std::cout << "Creating bead geometries... " << std::flush;
    for (std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++ )
    {
        double x = psf * (*iter)->getX();
        double y = psf * (*iter)->getY();
        double z = psf * (*iter)->getZ();
        double r = psf * (*iter)->getR();

        tag = factory::addSphere(x, y, z, this->prm->rFactor * r);
        ctag = factory::addPoint(x, y, z, this->prm->lc_beads, -1);

        (*iter)->setTag(tag);
        (*iter)->setCTag(ctag);

        dimTagsBeads.push_back({3, tag});
        tBeadCPs.push_back(ctag);

    }
    std::cout << "done!" << std::endl;

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
            double dx3, dy3, dz3;

            double x1 = psf * (*iter)->getX();
            double x2 = psf * (*riter)->getX();

            double y1 = psf * (*iter)->getY();
            double y2 = psf * (*riter)->getY();

            double z1 = psf * (*iter)->getZ();
            double z2 = psf * (*riter)->getZ();

            double r1 = psf * (*iter)->getR();
            double r2 = psf * (*riter)->getR();

            dx = x2-x1;
            dy = y2-y1;
            dz = z2-z1;

            double dist = sqrt( pow(dx, 2) + pow(dy, 2) + pow(dz, 2) );

            if (dist < r1 + r2 + this->prm->bridgeTol)
            {

                //cylinder radius
                double rBeadSmallest = (r1 <= r2? r1 : r2);
                double rBridge = this->prm->relativeBridgeRadius * rBeadSmallest;

                /* this->dimTagsBridges.push_back({3, factory::addCylinder((*iter)->getX(), (*iter)->getY(), (*iter)->getZ(), dx, dy, dz, rBridge) }) ; */
                /* this->dimTagsBeads.push_back({3, factory::addCylinder((*iter)->getX(), (*iter)->getY(), (*iter)->getZ(), dx, dy, dz, rBridge) }) ; */

                double factor = this->prm->bridgeOffsetRatio;

                //Cylinder start point
                x3 = x1 + factor * r1 * dx/dist;
                y3 = y1 + factor * r1 * dy/dist;
                z3 = z1 + factor * r1 * dz/dist;

                //cylinder length
                dx3 = dx * (dist - factor * (r1+r2))/dist;
                dy3 = dy * (dist - factor * (r1+r2))/dist;
                dz3 = dz * (dist - factor * (r1+r2))/dist;

                /* this->dimTagsBeads.push_back({3, factory::addCylinder(x3, y3, z3, dx3, dy3, dz3, rBridge) }) ; */
                this->dimTagsBridges.push_back({3, factory::addCylinder(x3, y3, z3, dx3, dy3, dz3, rBridge) }) ;

            }

        }
    }

    std::cout << "done!" << std::endl;

}

void PackedBed::mesh(std::string outfile)
{
    // Store dimTags and dimTagsMap
    std::vector<std::pair<int, int> > ov;
    std::vector<std::pair<int, int> > bv;
    std::vector<std::pair<int, int> > cv;
    std::vector<std::vector<std::pair<int, int> > > ovv;

    /* //TODO:extract wall, cut, */

    /* //Fuse beads together (with bridges or without) */
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

    // Fragment cylinder w.r.t. beads
    if (this->prm->fragment)
    {
        if (ov.size() == 0) ov = dimTagsBeads;

        std::cout << "Fragmenting Volumes... " << std::flush;
        factory::fragment(dimTagsCyl, ov, bv, ovv );
        std::cout << "done!" << std::endl;
    }

    // Synchronize gmsh model with geometry kernel.
    std::cout << "synchronizing... " << std::flush;
    factory::synchronize();
    std::cout << "done!" << std::endl;

    std::cout << std::endl;

    // ============================
    // Named Physical Groups

    dimTagsInterstitial.push_back(bv.back());

    if (this->prm->NamedInterstitialVolume)
    {
        model::addPhysicalGroup(3,{bv.back().second},5);
        model::setPhysicalName(3,5,"interstitialVolume");
    }

    if (this->prm->NamedOuterSurface)
    {
        model::getBoundary({{3, bv.back().second}}, cv, false, false,false);
        model::addPhysicalGroup(2,{cv[0].second},3);
        model::addPhysicalGroup(2,{cv[1].second},2);
        model::addPhysicalGroup(2,{cv[2].second},1);
        model::setPhysicalName(2,1, "inlet");
        model::setPhysicalName(2,2, "outlet");
        model::setPhysicalName(2,3, "wall");
    }


    bv.pop_back();

    if (this->prm->NamedBeadVolume)
    {
        tBeads.clear();
        for ( std::vector<std::pair< int , int>>::iterator it = bv.begin(); it != bv.end(); it++   )
        {
            tBeads.push_back((*it).second);
        }
        model::addPhysicalGroup(3,tBeads,6 );
        model::setPhysicalName(3,6,"beadVolume");

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

    //Set mesh size globally
    /* model::getEntities(cv, 0); */
    /* model::mesh::setSize(cv, this->prm->lc); */

    /* model::getEntitiesInBoundingBox(-0.1,-0.1,7.2, 0.1,0.1,7.8, cv, 0); */
    /* model::mesh::setSize(cv, this->prm->lc_beads); */

    // Set mesh size for beads
    /* model::getBoundary(bv, cv, false, false, true); */
    /* model::mesh::setSize(cv, this->prm->lc_beads); */

    /* model::mesh::embed(0,{tagCenter}, 3,dimTagsBeads.back().second); */

    //Embed points in beads to control element size within them
    if (this->prm->booleanOperation == 1)
    {
        for (std::vector<std::pair<int,int>>::iterator iter = bv.begin(); iter != bv.end(); iter++ )
        {
            model::mesh::embed(0,tBeadCPs, 3, (*iter).second);
        }
    }
    else
    {
        for (std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++ )
        {
            model::mesh::embed(0,{(*iter)->getCTag()}, 3, (*iter)->getTag());
        }

    }

    // Set mesh size for interstitial
    model::getBoundary(dimTagsInterstitial, cv, false, false, true);
    model::mesh::setSize(cv, this->prm->lc);

    //set mesh size on bead surface
    // maybe use the inside of the interstitial fragment instead?
    /* model::getBoundary(dimTagsInterstitial, cv, false, false, false); */
    /* model::getBoundary(bv, cv, false, false, false); */
    /* model::getBoundary(cv, ov, false, false, true); */
    /* model::mesh::setSize(ov, this->prm->lc); */

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

        if (this->prm->NamedBeadSurface)
        {
            std::cout << "Naming Bead Surfaces... ";
            model::getBoundary(bv, cv, false, false,false);
            tBeads.clear();
            for ( std::vector<std::pair< int , int>>::iterator it = cv.begin(); it != cv.end(); it++   )
            {
                tBeads.push_back((*it).second);
            }
            model::addPhysicalGroup(2,tBeads,14);
            model::setPhysicalName(2,14, "beadSurface");
            std::cout << "done!" << std::endl;
        }

        gmsh::write(this->prm->outpath + outfile + "_beads.vtk");


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

void PackedBed::readFile(std::string packingFilename)
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
    std::cout << "Selecting beads in range..." << std::endl;

    double cz1 = this->prm->zBot;
    double cz2 = this->prm->zTop;

    double nBeadsMax = data.size()/4;

    for(size_t i = 0; i < nBeadsMax; ++i)
    {
        double x = data[i * 4 ];
        double y = data[i * 4 + 1];
        double z = data[i * 4 + 2];
        double r = data[i * 4 + 3] * 0.5;

        beads.push_back(new Bead(x, y, z, this->prm->rFactor * r));
    }

    // sort using a lambda expression
    std::cout << "Sorting beads according to z-value..." << std::flush;
    std::sort(beads.begin(), beads.end(), [](const Bead* b1, const Bead* b2) {
            return b1->getZ() < b2->getZ();
    });
    std::cout << "done!" << std::endl;

    //Only store nBeads
    // testing nBeads
    double nBeads = this->prm->nBeads > nBeadsMax? nBeadsMax: this->prm->nBeads;
    beads.erase(beads.begin() + nBeads, beads.end());

    //Adjust zBot and zTop
    this->prm->zBot=beads.front()->getZ();
    this->prm->zTop=beads.back()->getZ();

    if (beads.size() == 0) {
        std::cerr << "ERROR: no beads in selection!!!" << std::endl;
        exit(1);
    }

    double radius_total = 0.0;
    double radius_max = -1.0;

    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
    {
        // remember: this is already the shrunken radius!!!
        double r = (*it)->getR();
        radius_total +=  r;
        if(r > radius_max)
            radius_max = r;
    }

    std::cout << this->beads.size() << "/" << data.size()/4 << " beads in the selected range." << std::endl;
    std::cout << "average bead radius: " << std::scientific << std::setprecision(10) << radius_total/beads.size() << std::endl;
    std::cout << "maximum bead radius: " << std::scientific << std::setprecision(10) << radius_max << std::endl << std::endl;

}
