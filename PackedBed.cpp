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
#include <bits/stdc++.h>

#define PI 3.1415926535897932

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;


PackedBed::PackedBed(Parameters * prm)
{
    // extract bead data from packing file.
    this->getBeads(prm);

    // transform packing data by scaling and offsetting the bed.
    if (this->beads.size() != 0)
        this->transformBeads(prm);

    // fix or calculate porosity
    if (prm->por_eps != DBL_MAX && prm->por_target != 0.0)
        this->fixPorosity(prm);
    else
        this->calcPorosity(prm);

    std::ofstream beadsfile;

    beadsfile.open(prm->outpath + "/beads.double.xyzd", std::ios::out | std::ios::binary);
    for (std::vector<Bead *>::iterator it = beads.begin(); it!=beads.end(); it++)
    {
        auto x = (*it)->x;
        auto y = (*it)->y;
        auto z = (*it)->z;
        auto d = (*it)->r * 2;

        beadsfile.write( reinterpret_cast<char *>(&x), sizeof(double));
        beadsfile.write( reinterpret_cast<char *>(&y), sizeof(double));
        beadsfile.write( reinterpret_cast<char *>(&z), sizeof(double));
        beadsfile.write( reinterpret_cast<char *>(&d), sizeof(double));
    }

    /* writeVecToBin(beads, beadsfile); */
    beadsfile.close();

    updateBounds(prm);
    // print data to stdout
    geometryStats(prm);

    if ((prm->periodic == "xy") || (prm->periodic == "xyz"))
        stackPeriodicPacking(prm);

    model::add("PackedBed");
}

PackedBed::~PackedBed()
{

}

template<typename T>
std::vector<double> PackedBed::readPacking(Parameters * prm)
{
    /* Read packing data */
    std::ifstream file(prm->packfile.c_str(), std::ios::binary);
    file.unsetf(std::ios::skipws);

    file.seekg(0, std::ios::end);
    const size_t num_elements = file.tellg() / prm->packingPrecision;
    file.seekg(0, std::ios::beg);
    std::vector<T> data(num_elements);

    std::cout << "reading packing ... " << std::flush;
    file.read(reinterpret_cast<char*>(&data[0]), num_elements*sizeof(T));
    if (isBigEndian())
        swapbytes(reinterpret_cast<char*>(&data[0]), data.size(), sizeof (T));

    std::vector<double> doubleVec(data.begin(), data.end());

    std::cout << "done!" << std::endl << std::endl;
    return doubleVec;
}


void PackedBed::getBeads(Parameters * prm)
{
    std::vector<double> data;
    if (prm->packingPrecision == 4)
        data = readPacking<float>(prm);
    else if (prm->packingPrecision == 8)
        data = readPacking<double>(prm);
    else
    {
        std::cout << "Invalid Packing Precision. Should be 4 or 8." << std::endl;
        exit(-1);
    }

    /* Select beads by number or range [zBot,zTop] */
    std::cout << "Selecting beads in range..." << std::flush;

    double cz1 = prm->zBot;
    double cz2 = prm->zTop;

    nBeadsMax = data.size()/4;

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

    // sort by z coordinate using a lambda expression
    std::cout << "Sorting beads according to z-value... " << std::flush;
    std::sort(beads.begin(), beads.end(), [](const Bead* b1, const Bead* b2) {
            return b1->getZ() < b2->getZ();
    });
    std::cout << "done!" << std::endl;

    // Only store nBeads
    if(prm->nBeads >= 0)
    {
        nBeads = prm->nBeads > nBeadsMax? nBeadsMax: prm->nBeads;
        beads.erase(beads.begin() + nBeads, beads.end());
    }

    std::cout << this->beads.size() << "/" << nBeadsMax << " beads in the selected range." << std::endl;

    if (this->beads.size() < 0)
        return

    updateBounds(prm);
    printBounds();

    std::cout << "Scaling bead radii in place... " << std::flush;
    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
        (*it)->scaleRadius(prm->rFactor);
    std::cout << "done!" << std::endl;

}

void PackedBed::transformBeads(Parameters * prm)
{
    updateBounds(prm);

    std::cout << std::endl;
    std::cout << "== Before Transform ==" << std::endl;
    printBounds();

    double offsetx;
    double offsety;
    double offsetz;

    if (prm->translateOffsets == "auto")
    {
        std::cout << "Calculating translation offsets automatically from bounding box." << std::endl;
        offsetx = -xCyl;
        offsety = -yCyl;
        offsetz = -zBot;
    }
    else
    {
        std::cout << "Using provided translation offsets." << std::endl;
        offsetx = prm->tOffX;
        offsety = prm->tOffY;
        offsetz = prm->tOffZ;
    }

    // Use computed offset values to autoscale and translate packed bed
    // This way the bottom of the bed coincides with the x-y plane
    // And the bead diameter should correspond to 1 unit.
    std::cout << "Translating beads by (" << offsetx << ", "<< offsety << ", " << offsetz << ")... "<<  std::flush;
    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
        (*it)->translate(offsetx, offsety, offsetz);
    std::cout << "done!" << std::endl;

    std::cout << "Scaling beads by " << prm->preScalingFactor << "... " << std::flush;
    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
        (*it)->scale(prm->preScalingFactor);
    std::cout << "done!" << std::endl << std::endl;

    updateBounds(prm);

    prm->xCyl = xCyl;
    prm->yCyl = yCyl;
    prm->rCyl = rCyl;

    if ( prm->refBeadSize == "avg")
        prm->refBeadRadius = radius_avg;
    else if ( prm->refBeadSize == "max" )
        prm->refBeadRadius = radius_max;
    else if ( prm->refBeadSize == "min" )
        prm->refBeadRadius = radius_min;

    assert ( rCyl >= std::max( ((xMax - xMin) / 2), ((yMax - yMin)/2) ));

    //TODO: Can be moved to updateBounds
    if (prm->nBeads < 0)
    {
        // if packed bed is sliced based on zTop and zBot instead of a limit on number of beads,
        // apply transforms on the zTop, zBot inputs rather than calculating them from the post-transform geometry
        // this way, the zTop-zBot (and hence column length) is pre-determined by the inputs even though the packed bed's
        // actual size between different models may vary. This is helpful with generating "identical" columns of mono and poly beads.

        zCylBot = (prm->zBot + offsetz) * prm->preScalingFactor - prm->inlet;
        zCylTop = (prm->zTop + offsetz) * prm->preScalingFactor + prm->outlet;
    }
    else
    {
        zCylBot = zBot - prm->inlet;
        zCylTop = zTop + prm->outlet;
    }

    std::cout << "== After Transform ==" << std::endl;
    printBounds();

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

void PackedBed::updateBounds(Parameters * prm)
{
    // NOTE: Changes zTop and zBot, will affect BED porosity calculations

    //reset variables
    xCyl = 0, yCyl = 0;
    xMax = -DBL_MAX, yMax = -DBL_MAX, zMax = -DBL_MAX;
    xMin = DBL_MAX, yMin = DBL_MAX, zMin = DBL_MAX;
    radius_avg=0;
    radius_max= -DBL_MAX;
    radius_min= DBL_MAX;

    zBot = DBL_MAX;
    zTop = -DBL_MAX;

    /* vol_real_beads=0; */
    /* vol_geom_beads=0; */

    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
    {
        double x=0, y=0, r=0, z=0;
        x = (*it)->getX();
        y = (*it)->getY();
        z = (*it)->getZ();
        r = (*it)->getR();

        radius_avg +=  r;
        /* vol_real_beads += PI * 4/3 * pow(r * prm->MeshScalingFactor * 1/(prm->rFactor), 3); */
        /* vol_geom_beads += PI * 4/3 * pow(r * prm->MeshScalingFactor, 3); */

        if ( (x + r) > xMax ) xMax = x + r;
        if ( (x - r) < xMin ) xMin = x - r;
        if ( (y + r) > yMax ) yMax = y + r;
        if ( (y - r) < yMin ) yMin = y - r;
        if ( (z + r) > zMax ) zMax = z + r;
        if ( (z - r) < zMin ) zMin = z - r;

        // Note: These are local values of zTop and zBot.
        // Different from prm->zBot,zTop, which are used as slice coordinates and part of the zCylBot/Top calculations.
        // Local zBot,Top are used to find bed lengths.
        // They can be different because of porosity manipulation/control.
        // They can be different (slightly) if nBeads > 0.
        if ( z < zBot ) zBot = z;
        if ( z > zTop ) zTop = z;

        if (r > radius_max) radius_max = r;
        if (r < radius_min) radius_min = r;

    }
    radius_avg /= beads.size();

    xCyl = (xMax + xMin) / 2;
    yCyl = (yMax + yMin) / 2;
    rCyl = std::max( ((xMax - xMin) / 2), ((yMax - yMin)/2) ) + prm->rCylDelta;

}

void PackedBed::printBounds()
{
    std::cout << "xCyl: " << xCyl << std::endl;
    std::cout << "yCyl: " << yCyl << std::endl;
    std::cout << "rCyl: " << rCyl << std::endl;
    std::cout << "xMax: " << xMax << std::endl;
    std::cout << "xMin: " << xMin << std::endl;
    std::cout << "yMax: " << yMax << std::endl;
    std::cout << "yMin: " << yMin << std::endl;
    std::cout << "zMax: " << zMax << std::endl;
    std::cout << "zMin: " << zMin << std::endl;
    std::cout << "zBot: " << zBot << std::endl;
    std::cout << "zTop: " << zTop << std::endl;
    std::cout << std::endl;

    std::cout << "average bead radius: " << radius_avg << std::endl;
    std::cout << "maximum bead radius: " << radius_max << std::endl;
    std::cout << "minimum bead radius: " << radius_min << std::endl << std::endl;
}

void PackedBed::calcPorosity(Parameters * prm)
{
    std::cout << "Calculating porosity..." << std::endl;

    vol_real_beads = 0;
    vol_geom_beads = 0;

    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
    {
        double r = (*it)->getR();
        vol_real_beads += PI * 4/3 * pow(r * prm->MeshScalingFactor * 1/(prm->rFactor), 3);
        vol_geom_beads += PI * 4/3 * pow(r * prm->MeshScalingFactor, 3);
    }

    vol_cylinder = PI * pow(rCyl * prm->MeshScalingFactor, 2) * (zCylTop - zCylBot) * prm->MeshScalingFactor;
    vol_real_int = vol_cylinder - vol_real_beads;

    por_real_col = vol_real_int / vol_cylinder;                         // real (ideal) packing porosity of the full column.
    por_geom_col = (vol_cylinder - vol_geom_beads) / vol_cylinder;      // packing porosity of the geometry of the column after bead shrinking (without bridges)

    bedLength = (zMax-zMin)*prm->MeshScalingFactor;                             // length of packed bed region (Tip to tip)
    vol_bed_cyl = PI * pow(rCyl * prm->MeshScalingFactor, 2) * bedLength;       // volume of the cylinder corresponding to the length of the packed bed.

    por_real_bed  =  ((vol_bed_cyl - vol_real_beads))  / vol_bed_cyl;           // ideal porosity of the packed bed region
    por_geom_bed  =  ((vol_bed_cyl - vol_geom_beads)) / vol_bed_cyl;            // porosity of the packed bed region after bead shrinking

    std::cout << "Real Column Porosity: " << por_real_col << std::endl;         // ideal porosity of the full column

    /* std::cout << "Modified Column Porosity (without bridges): " << por_geom_col << std::endl << std::endl; */

    /* NOTE: zMax and zMin might change in the scenario where porosity control happens */
    std::cout << "Bed Length (zMax - zMin): " << (zMax-zMin)*prm->MeshScalingFactor << std::endl;
    std::cout << "Bed Cylinder Volume: " << vol_bed_cyl << std::endl;
    std::cout << "Real Bed Porosity: " << por_real_bed << std::endl;
    std::cout << "Modified Bed Porosity (without bridges): " << por_geom_bed << std::endl << std::endl;


}

void PackedBed::geometryStats(Parameters * prm)
{
    std::cout << std::endl;
    std::cout << "=== After Mesh Scaling ===" << std::endl;

    std::cout << "Mesh Scaling Factor: " << prm->MeshScalingFactor              << std::endl;
    std::cout << "Avg Bead Radius: "     << radius_avg * prm->MeshScalingFactor << std::endl;
    std::cout << "Cylinder Radius: "     << rCyl * prm->MeshScalingFactor       << std::endl;
    std::cout << "Cylinder Volume: "     << vol_cylinder                        << std::endl;
    std::cout << "Real Int Volume: "     << vol_real_int                        << std::endl;
    std::cout << "Real Bead Volume: "    << vol_real_beads                      << std::endl;

    std::cout << "Modified Bead Volume: (without bridges) " << vol_geom_beads << std::endl<< std::endl;
    std::cout << "Bead Geometry Volume Error: " << (vol_real_beads - vol_geom_beads)/vol_real_beads*100 << "%" << std::endl << std::endl;

    std::cout << "Column Length: " << (zCylTop - zCylBot) * prm->MeshScalingFactor << std::endl<<std::endl;
    std::cout << "zTop - zBot (given): " << (prm->zTop - prm->zBot)*prm->preScalingFactor * prm->MeshScalingFactor << std::endl;
    std::cout << "zTop - zBot: " << (zTop - zBot)*prm->MeshScalingFactor << std::endl;
    std::cout << "zMax - zMin: " << (zMax - zMin)*prm->MeshScalingFactor << std::endl<< std::endl;;

    // TODO: Consider adding porosity values here instead of in fix/calcporosity?

}

void PackedBed::fixPorosity(Parameters * prm)
{
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>" <<std::endl;
    std::cout << "Modifying column porosity..." <<std::endl;
    calcPorosity(prm);

    std::vector<std::pair<int,Bead*>> endZoneBeads ;
    std::vector<double> endZoneBeadRads;

    if (prm->fixPorosityMethod == 0)
    {
        // useful for removing beads based on volume. But changes the bulk of the packed bed.
        std::cout << "Sorting beads according to r-value... " << std::flush;
        std::sort(beads.begin(), beads.end(), [](const Bead* b1, const Bead* b2) {
                return b1->getR() < b2->getR();
                });
        std::cout << "done!" << std::endl;
    }
    else if (prm->fixPorosityMethod == 2)
    {
        std::cout << "Finding End Zone Beads..." << std::flush;
        for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
            if (((*it)->getZ() >= zTop - radius_max) || ((*it)->getZ() <= zBot + radius_max) )
                endZoneBeads.push_back(std::pair<int,Bead*>(std::distance(beads.begin(),it), *it));
        std::cout << "done!" << std::endl;

        std::cout << "Sorting endZoneBeads..." << std::flush;
        std::sort(endZoneBeads.begin(), endZoneBeads.end(), [](const std::pair<int,Bead*> b1, const std::pair<int,Bead*> b2) {
                return b1.second->getR() < b2.second->getR();
                });
        std::cout << "done!" << std::endl;

        for (std::vector<std::pair<int,Bead*>>::iterator it = endZoneBeads.begin(); it != endZoneBeads.end(); it++)
            endZoneBeadRads.push_back((*it).second->getR() * prm->MeshScalingFactor);

    }

    std::vector<double> vBeadRads;

    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
    {
        double radscaled= (*it)->getR() * prm->MeshScalingFactor;
        vBeadRads.push_back(radscaled);
    }

    while ((prm->por_target - por_real_col) > prm->por_eps)
    {

        double vol_target_beads;

        vol_target_beads = - vol_cylinder * (prm->por_target - 1);
        double vol_bead_removal = vol_real_beads - vol_target_beads;
        double rad_bead_removal = pow( (3*fabs(vol_bead_removal))/(4*PI), 1.0/3.0);

        std::cout << "Volume to be removed: " << vol_bead_removal << std::endl;
        std::cout << "Radius to be removed: " << rad_bead_removal << std::endl;

        if (prm->fixPorosityMethod == 0)
        {
            // Requires sorting!!
            // useful for removing beads based on volume. But changes the bulk of the packed bed.
            int index = findBeadWithRadius(rad_bead_removal, vBeadRads);
            std::cout << "==> Removing bead " << index << " with radius " << vBeadRads[index] << std::endl;
            vBeadRads.erase(vBeadRads.begin() + index);
            beads.erase(beads.begin() + index);
        }
        else if (prm->fixPorosityMethod == 1)
        {
            std::cout << "Removing last bead..." << std::endl;
            vBeadRads.pop_back();
            beads.pop_back();
        }
        else if (prm->fixPorosityMethod == 2)
        {
            std::cout << "Radius Max is " << radius_max << std::endl;


            int indexEndZoneRads = findBeadWithRadius(rad_bead_removal, endZoneBeadRads);
            int indexvBeadRads = endZoneBeads[indexEndZoneRads].first;

            endZoneBeadRads.erase(endZoneBeadRads.begin() + indexEndZoneRads);
            endZoneBeads.erase(endZoneBeads.begin() + indexEndZoneRads);
            vBeadRads.erase(vBeadRads.begin() + indexvBeadRads);
            beads.erase(beads.begin() + indexvBeadRads);

        }
        else
        {
            std::cout << "Invalid fixPorosityMethod input!" << std::endl;
            exit(-1);
        }

        calcPorosity(prm);
        std::cout << "Target porosity: " << prm->por_target << std::endl;
    }

    //unsort??
    /* if (prm->fixPorosityMethod == 0) */
    /* { */
    /*     // useful for removing beads based on volume. But changes the bulk of the packed bed. */
    /*     std::cout << "Sorting Beads by z-value..." << std::flush; */
    /*     std::sort(beads.begin(), beads.end(), [](const Bead* b1, const Bead* b2) { */
    /*             return b1->getR() < b2->getR(); */
    /*             }); */
    /*     std::cout << "done!" << std::endl; */
    /* } */


    std::cout << beads.size() << " beads remaining after porosity control." << std::endl;

    std::cout << "Porosity modification complete!" <<std::endl;
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" <<std::endl << std::endl;

    updateBounds(prm);
    calcPorosity(prm);
    printBounds();


}

int PackedBed::findBeadWithRadius(double value, std::vector<double> vBeadRads)
{
    auto it = std::lower_bound(vBeadRads.begin(), vBeadRads.end(), value, [](double a, double b){return a < b;});
    int index = std::distance(vBeadRads.begin(), it);

    //if value > radius of last bead, index = last bead
    if (it == vBeadRads.end())
        index -= 1;

    return index;
}

void PackedBed::stackPeriodicPacking(Parameters * prm)
{
    std::cout << "Stacking periodic packings in " << prm->periodic << " directions... " << std::endl;

    /* double eps = 1e-3;  // value of tolerated overlap */
    double eps = 0;  // value of tolerated overlap

    std::vector<Bead *> xmBeads;
    std::vector<Bead *> xpBeads;
    std::vector<Bead *> ymBeads;
    std::vector<Bead *> ypBeads;

    std::vector<Bead *> xpypBeads;
    std::vector<Bead *> xpymBeads;
    std::vector<Bead *> xmymBeads;
    std::vector<Bead *> xmypBeads;

    std::vector<Bead *> stackedBeads;

    double xoff, yoff, zoff;

    if (prm->periodicOffsets == "auto")
    {
        if (prm->autoContainment == 1)
        {
            std::cout << "WARNING: Using autoContainment generated bounds for stacking periodic packed beds. NOT RECOMMENDED FOR PERIODICITY." << std::endl;
            xoff = (this->xMax-this->xMin); // Offset it completely... effectively stack the outer containers
            yoff = (this->yMax-this->yMin);
            zoff = (this->zMax-this->zMin);
        }
        else if (prm->autoContainment == 0) // manual containment
        {
            xoff = (prm->dx);               // Offset it completely using box deltas... effectively stack the outer containers
            yoff = (prm->dy);
            zoff = (prm->dz);
        }

    }
    else
    {
        xoff = prm->pOffX;
        yoff = prm->pOffY;
        zoff = prm->pOffZ;
    }

    std::vector<int> xOffsetMultiplier = {0};
    std::vector<int> yOffsetMultiplier = {0};
    std::vector<int> zOffsetMultiplier = {0};

    /* if (prm->periodic == "xyz") */
    /*     zOffsetMultiplier = {-1, 0, 1}; */
    /* else */
    /*     zOffsetMultiplier = {0}; */

    std::size_t found;
    found = prm->periodic.find('x');
    if (found != std::string::npos) xOffsetMultiplier = {-1, 0, 1};

    found = prm->periodic.find('y');
    if (found != std::string::npos) yOffsetMultiplier = {-1, 0, 1};

    found = prm->periodic.find('z');
    if (found != std::string::npos) zOffsetMultiplier = {-1, 0, 1};

    for (std::vector<int>::iterator zom = zOffsetMultiplier.begin(); zom != zOffsetMultiplier.end(); zom++)
    {
        for (std::vector<int>::iterator yom = yOffsetMultiplier.begin(); yom != yOffsetMultiplier.end(); yom++)
        {
            for (std::vector<int>::iterator xom = xOffsetMultiplier.begin(); xom != xOffsetMultiplier.end(); xom++)
            {
                if ( (*xom == 0) && (*yom == 0) && (*zom == 0)) continue;
                for (std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++)
                    stackedBeads.push_back(new Bead((*iter)->getX() + (*xom) * xoff, (*iter)->getY() + (*yom) * yoff, (*iter)->getZ() + (*zom) * zoff, (*iter)->getR()));
            }
        }
    }

    this->beads.reserve(this->beads.size() + stackedBeads.size());
    this->beads.insert(this->beads.end(), stackedBeads.begin(), stackedBeads.end());


    /* // Full copy beads with offset in the x-y directions */
    /* for (std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++) */
    /*     xpBeads.push_back( new Bead((*iter)->getX() + xoff, (*iter)->getY(), (*iter)->getZ(), (*iter)->getR()) ); */

    /* for (std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++) */
    /*     xmBeads.push_back( new Bead((*iter)->getX() - xoff, (*iter)->getY(), (*iter)->getZ(), (*iter)->getR()) ); */

    /* for (std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++) */
    /*     ypBeads.push_back( new Bead((*iter)->getX(), (*iter)->getY() + yoff, (*iter)->getZ(), (*iter)->getR()) ); */

    /* for (std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++) */
    /*     ymBeads.push_back( new Bead((*iter)->getX(), (*iter)->getY() - yoff, (*iter)->getZ(), (*iter)->getR()) ); */

    /* for (std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++) */
    /*     xpypBeads.push_back( new Bead((*iter)->getX() + xoff, (*iter)->getY() + yoff, (*iter)->getZ(), (*iter)->getR()) ); */

    /* for (std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++) */
    /*     xmypBeads.push_back( new Bead((*iter)->getX() - xoff, (*iter)->getY() + yoff, (*iter)->getZ(), (*iter)->getR()) ); */

    /* for (std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++) */
    /*     xmymBeads.push_back( new Bead((*iter)->getX() - xoff, (*iter)->getY() - yoff, (*iter)->getZ(), (*iter)->getR()) ); */

    /* for (std::vector<Bead *>::iterator iter = this->beads.begin(); iter != this->beads.end(); iter++) */
    /*     xpymBeads.push_back( new Bead((*iter)->getX() + xoff, (*iter)->getY() - yoff, (*iter)->getZ(), (*iter)->getR()) ); */

    /* // Reserve space to extend this->beads vector */
    /* this->beads.reserve(this->beads.size() * 9); */

    /* // Extend this->beads vector */
    /* this->beads.insert(this->beads.end(), xpBeads.begin(), xpBeads.end()); */
    /* this->beads.insert(this->beads.end(), xmBeads.begin(), xmBeads.end()); */
    /* this->beads.insert(this->beads.end(), ypBeads.begin(), ypBeads.end()); */
    /* this->beads.insert(this->beads.end(), ymBeads.begin(), ymBeads.end()); */

    /* this->beads.insert(this->beads.end(), xpypBeads.begin(), xpypBeads.end()); */
    /* this->beads.insert(this->beads.end(), xmypBeads.begin(), xmypBeads.end()); */
    /* this->beads.insert(this->beads.end(), xmymBeads.begin(), xmymBeads.end()); */
    /* this->beads.insert(this->beads.end(), xpymBeads.begin(), xpymBeads.end()); */

}

double PackedBed::calculateMinDistance(std::vector<Bead *> _beads)
{
    double minDistance = INT_MAX;
    // FIXME: Inefficient
    for (std::vector<Bead *>::iterator iteri = this->beads.begin(); iteri != this->beads.end(); iteri++)
        for (std::vector<Bead *>::iterator iterj = _beads.begin(); iterj != _beads.end(); iterj++)
        {
            double distance = sqrt(
                    ( (*iteri)->getX() - (*iterj)->getX() ) * ( (*iteri)->getX() - (*iterj)->getX() ) +
                    ( (*iteri)->getY() - (*iterj)->getY() ) * ( (*iteri)->getY() - (*iterj)->getY() ) +
                    ( (*iteri)->getZ() - (*iterj)->getZ() ) * ( (*iteri)->getZ() - (*iterj)->getZ() )
                    );

            /* std::cout << iteri - this->beads.begin() << ", " << iterj - _beads.begin() << ": " << distance << std::endl; */

            if ( minDistance > (distance - ( (*iteri)->getR() + (*iterj)->getR() )) )
                minDistance = distance - ( (*iteri)->getR() + (*iterj)->getR() );

        }
    return minDistance;

}
