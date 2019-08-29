/*
 * File: PackedBed.cpp
 * Created by: Rao, Jayghosh Subodh
 * Created on: Thu 08 Aug 2019 12:57:40 PM CEST
 *
 */

#include "PackedBed.h"
#include "Parameters.h"
#include "Files.h"

#include <gmsh.h>
#include <assert.h>

#define PI 3.1415926535897932

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;


PackedBed::PackedBed(Parameters * prm)
{
    this->getBeads(prm);

    /* beads.clear(); */
    /* beads.push_back(new Bead(5, 5, 710.5, prm->rFactor * 1)); */
    /* beads.push_back(new Bead(6.5, 5,710.5, prm->rFactor * 0.50000)); */

    this->transformBeads(prm);
    model::add("PackedBed");
}

PackedBed::~PackedBed()
{

}


void PackedBed::getBeads(Parameters * prm)
{
    // open the file:
    std::ifstream file(prm->packfile.c_str(), std::ios::binary);

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

    double cz1 = prm->zBot;
    double cz2 = prm->zTop;

    nBeadsMax = data.size()/4;

    //Read and store all beads if nbeads > 0
    //else use zBot and zTop limits
    for(size_t i = 0; i < nBeadsMax; ++i)
    {
        double x = data[i * 4 ];
        double y = data[i * 4 + 1];
        double z = data[i * 4 + 2];
        double r = data[i * 4 + 3] * 0.5;

        if(prm->nBeads < 0)
        {
            if (z >= cz1 && z <= cz2)
                beads.push_back(new Bead(x, y, z, r));
        }
        else
            beads.push_back(new Bead(x, y, z, r));

    }
    std::cout << "done!" << std::endl;


    if (beads.size() == 0) {
        std::cerr << "ERROR: No beads found!" << std::endl;
        exit(1);
    }

    // sort using a lambda expression
    std::cout << "Sorting beads according to z-value... " << std::flush;
    std::sort(beads.begin(), beads.end(), [](const Bead* b1, const Bead* b2) {
            return b1->getZ() < b2->getZ();
    });
    std::cout << "done!" << std::endl;

    //Only store nBeads
    if(prm->nBeads >= 0)
    {
        nBeads = prm->nBeads > nBeadsMax? nBeadsMax: prm->nBeads;
        beads.erase(beads.begin() + nBeads, beads.end());
    }

    std::cout << this->beads.size() << "/" << nBeadsMax << " beads in the selected range." << std::endl;

    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
    {
        x = (*it)->getX();
        y = (*it)->getY();
        r = (*it)->getR();

        if ( (x + r) > xMax ) xMax = x + r;
        if ( (x - r) < xMin ) xMin = x - r;
        if ( (y + r) > yMax ) yMax = y + r;
        if ( (y - r) < yMin ) yMin = y - r;

    }
    prm->rCyl = std::max( ((xMax - xMin) / 2), ((yMax - yMin)/2) );

    std::cout << "Scaling bead radii in place... " << std::flush;
    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
        (*it)->scaleRadius(prm->rFactor);
    std::cout << "done!" << std::endl;

}

void PackedBed::transformBeads(Parameters * prm)
{
    xCyl = 0, yCyl = 0;
    xMax = -DBL_MAX, yMax = -DBL_MAX, zMax = -DBL_MAX;
    xMin = DBL_MAX, yMin = DBL_MAX, zMin = DBL_MAX;
    x=0, y=0, r=0, z=0;
    zBot=0, zTop=0;
    radius_avg=0;

    rCyl = prm->rCyl;

    //Calculate bounding box
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
    /* rCyl = std::max( ((xMax - xMin) / 2), ((yMax - yMin)/2) ) + prm->rCylDelta; */
    /* rCyl = std::max(std::max(xMax, std::abs(xMin)), std::max(yMax, std::abs(yMin) ) ) + prm->rCylDelta; */

    //Adjust zBot and zTop
    zBot=beads.front()->getZ();
    zTop=beads.back()->getZ();

    /* std::cout << std::setprecision(4); */
    std::cout << std::endl;
    std::cout << "== Before Transform ==" << std::endl;
    std::cout << "xCyl: " << xCyl << std::endl;
    std::cout << "yCyl: " << yCyl << std::endl;
    std::cout << "rCyl: " << rCyl << std::endl;
    std::cout << "xMax: " << xMax << std::endl;
    std::cout << "xMin: " << xMin << std::endl;
    std::cout << "yMax: " << yMax << std::endl;
    std::cout << "yMin: " << yMin << std::endl;
    /* std::cout << std::setprecision(10); */
    std::cout << "zBot: " << zBot << std::endl;
    std::cout << "zTop: " << zTop << std::endl;
    std::cout << std::endl;

    double offsetx = -xCyl;
    double offsety = -yCyl;
    double offsetz = -zBot;

    //Use computed values to autoscale and translate packed bed
    std::cout << "Translating beads by (" << offsetx << ", "<< offsety << ", " << offsetz << ")... "<<  std::flush;
    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
        (*it)->translate(offsetx, offsety, offsetz);
    std::cout << "done!" << std::endl;

    std::cout << "Scaling beads by " << prm->preScalingFactor << "... " << std::flush;
    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
        (*it)->scale(prm->preScalingFactor);
    std::cout << "done!" << std::endl << std::endl;

    //Get new zBot and zTop
    prm->rCyl *= prm->preScalingFactor;
    rCyl = rCyl * prm->preScalingFactor + prm->rCylDelta;;
    zBot=beads.front()->getZ();
    zTop=beads.back()->getZ();

    //reset variables
    xCyl = 0, yCyl = 0;
    xMax = -DBL_MAX, yMax = -DBL_MAX, zMax = -DBL_MAX;
    xMin = DBL_MAX, yMin = DBL_MAX, zMin = DBL_MAX;
    x=0, y=0, r=0, z=0;
    radius_avg=0;

    double vol_real_beads=0;
    double vol_geom_beads=0;
    double vol_real_int=0;
    double vol_mesh_int=0;
    double vol_cylinder=0;

    //Get bounding parameters
    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
    {
        x = (*it)->getX();
        y = (*it)->getY();
        z = (*it)->getZ();
        r = (*it)->getR();

        radius_avg +=  r;
        vol_real_beads += PI * 4/3 * pow(r * prm->MeshScalingFactor * 1/(prm->rFactor), 3);
        vol_geom_beads += PI * 4/3 * pow(r * prm->MeshScalingFactor, 3);

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

    assert ( rCyl >= std::max( ((xMax - xMin) / 2), ((yMax - yMin)/2) ));


    /* rCyl = std::max( ((xMax - xMin) / 2), ((yMax - yMin)/2) ) + prm->rCylDelta; */
    /* rCyl = std::max(std::max(xMax, std::abs(xMin)), std::max(yMax, std::abs(yMin) ) ) + prm->rCylDelta; */

    vol_cylinder = PI * pow(rCyl * prm->MeshScalingFactor, 2) * (zTop - zBot + prm->inlet + prm->outlet) * prm->MeshScalingFactor;
    vol_real_int = vol_cylinder - vol_real_beads;

    prm->xCyl = xCyl;
    prm->yCyl = yCyl;
    prm->rCyl = rCyl;

    // If we choose to use zBot and zTop as inputs
    // make sure that the column length is the same
    // regardless of the preScalingFactor used for
    // different packings.
    if (prm->nBeads < 0)
    {
        prm->zBot += offsetz;
        prm->zTop += offsetz;
        prm->zBot *= prm->preScalingFactor;
        prm->zTop *= prm->preScalingFactor;
        zBot = prm->zBot;
        zTop = prm->zTop;
    }
    else
    {
        prm->zBot = zBot;
        prm->zTop = zTop;
    }

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

    /* std::cout << std::scientific << std::setprecision(4); */
    std::cout << "average bead radius: " << radius_avg << std::endl;
    std::cout << "maximum bead radius: " << radius_max << std::endl;
    std::cout << "minimum bead radius: " << radius_min << std::endl << std::endl;

    std::cout << "=== After Mesh Scaling ===" << std::endl;
    std::cout << "Mesh Scaling Factor: " << prm->MeshScalingFactor        << std::endl;
    std::cout << "Cylinder Radius: "     << rCyl * prm->MeshScalingFactor << std::endl;
    std::cout << "Cylinder Volume: "     << vol_cylinder                  << std::endl;
    std::cout << "Real Int Volume: "     << vol_real_int                  << std::endl;
    std::cout << "Real Bead Volume: "    << vol_real_beads                << std::endl;

    std::cout << "Modified Bead Volume: (without bridges) " << vol_geom_beads << std::endl<< std::endl;
    /* std::cout << std::fixed << std::setprecision(2); */
    std::cout << "Bead Geometry Volume Error: " << (vol_real_beads - vol_geom_beads)/vol_real_beads*100 << "%" << std::endl << std::endl;

    bedLength = (zMax-zMin)*prm->MeshScalingFactor;
    vol_bed_cyl = PI * pow(rCyl * prm->MeshScalingFactor, 2) * bedLength;
    por_real_bed  =  ((vol_bed_cyl - vol_real_beads))  / vol_bed_cyl;
    por_geom_bed  =  ((vol_bed_cyl - vol_geom_beads)) / vol_bed_cyl;

    /* std::cout << std::scientific << std::setprecision(4); */
    std::cout << "Column Length: " << (zTop - zBot + prm->inlet + prm->outlet) * prm->MeshScalingFactor << std::endl<<std::endl;
    std::cout << "Bed Length: " << (zMax-zMin)*prm->MeshScalingFactor << std::endl;
    std::cout << "zTop - zBot: " << (zTop-zBot)*prm->MeshScalingFactor << std::endl;
    std::cout << "Bed Cylinder Volume: " << vol_bed_cyl << std::endl;

    /* std::cout << std::fixed << std::setprecision(2); */
    std::cout << "Real Porosity: " << por_real_bed << std::endl;
    std::cout << "Modified Porosity (without bridges): " << por_geom_bed << std::endl << std::endl;


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
