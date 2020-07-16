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
    // extract bead data from packing file.
    this->getBeads(prm);

    // transform packing data by scaling and offsetting the bed.
    this->transformBeads(prm);

    if (prm->por_eps != DBL_MAX && prm->por_target != 0.0)
        this->fixPorosity(prm);
    else
        this->calcPorosity(prm);

    geometryStats(prm);

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

    //Only store nBeads
    if(prm->nBeads >= 0)
    {
        nBeads = prm->nBeads > nBeadsMax? nBeadsMax: prm->nBeads;
        beads.erase(beads.begin() + nBeads, beads.end());
    }

    std::cout << this->beads.size() << "/" << nBeadsMax << " beads in the selected range." << std::endl;

    updateBounds();

    std::cout << "Scaling bead radii in place... " << std::flush;
    for(std::vector<Bead*>::iterator it = beads.begin(); it != beads.end(); it++)
        (*it)->scaleRadius(prm->rFactor);
    std::cout << "done!" << std::endl;

}

void PackedBed::transformBeads(Parameters * prm)
{
    updateBounds();

    /* std::cout << std::setprecision(4); */
    std::cout << std::endl;
    std::cout << "== Before Transform ==" << std::endl;
    printBounds();

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


    updateBounds();

    rCyl += prm->rCylDelta;;

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

    if (prm->nBeads < 0)
    {
        // if packed bed is sliced based on zTop and zBot instead of a limit on number of beads,
        // apply transforms on the zTop, zBot inputs rather than calculating them from the post-transform geometry
        // this way, the zTop-zBot (and hence column length) is pre-determined by the inputs even though the packed bed's
        // actual size between different models may vary. This is helpful with generating "identical" columns of mono and poly beads.

        prm->zBot += offsetz;
        prm->zTop += offsetz;
        prm->zBot *= prm->preScalingFactor;
        prm->zTop *= prm->preScalingFactor;
            /* std::cout << "Scaled: " << prm->zTop << "\t" << prm->zBot << prm->zTop - prm->zBot << std::endl; */
            /* std::cout << "Calced: " << zTop << "\t" << zBot << zTop - zBot << std::endl; */
        zBot = prm->zBot;
        zTop = prm->zTop;
    }
    else
    {
        // update parameters object. This is used in the Model class to create the cylinder.
        prm->zBot = zBot;
        prm->zTop = zTop;
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

void PackedBed::updateBounds()
{
    //reset variables
    xCyl = 0, yCyl = 0;
    xMax = -DBL_MAX, yMax = -DBL_MAX, zMax = -DBL_MAX;
    xMin = DBL_MAX, yMin = DBL_MAX, zMin = DBL_MAX;
    radius_avg=0;
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

        if (r > radius_max) radius_max = r;
        if (r < radius_min) radius_min = r;

    }
    radius_avg /= beads.size();

    xCyl = (xMax + xMin) / 2;
    yCyl = (yMax + yMin) / 2;
    rCyl = std::max( ((xMax - xMin) / 2), ((yMax - yMin)/2) );

    //Adjust zBot and zTop
    zBot=beads.front()->getZ();
    zTop=beads.back()->getZ();


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

    vol_cylinder = PI * pow(rCyl * prm->MeshScalingFactor, 2) * (zTop - zBot + prm->inlet + prm->outlet) * prm->MeshScalingFactor;
    vol_real_int = vol_cylinder - vol_real_beads;

    por_real_col = vol_real_int / vol_cylinder;                         // real (ideal) packing porosity of the full column.
    por_geom_col = (vol_cylinder - vol_geom_beads) / vol_cylinder;      // packing porosity of the geometry of the column after bead shrinking (without bridges)

    bedLength = (zMax-zMin)*prm->MeshScalingFactor;                             // length of packed bed region (Tip to tip)
    vol_bed_cyl = PI * pow(rCyl * prm->MeshScalingFactor, 2) * bedLength;       // volume of the cylinder corresponding to the length of the packed bed.

    por_real_bed  =  ((vol_bed_cyl - vol_real_beads))  / vol_bed_cyl;           // ideal porosity of the packed bed region
    por_geom_bed  =  ((vol_bed_cyl - vol_geom_beads)) / vol_bed_cyl;            // porosity of the packed bed region after bead shrinking

    std::cout << "Real Column Porosity: " << por_real_col << std::endl;         // ideal porosity of the full column

        /* std::cout << "Bed Length (zMax - zMin): " << (zMax-zMin)*prm->MeshScalingFactor << std::endl; */
        /* std::cout << "Bed Cylinder Volume: " << vol_bed_cyl << std::endl; */
        /* /1* std::cout << std::fixed << std::setprecision(2); *1/ */
        /* std::cout << "Real Bed Porosity: " << por_real_bed << std::endl; */
        /* std::cout << "Modified Bed Porosity (without bridges): " << por_geom_bed << std::endl << std::endl; */

    /* std::cout << "Modified Column Porosity (without bridges): " << por_geom_col << std::endl << std::endl; */

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
    /* std::cout << std::fixed << std::setprecision(2); */
    std::cout << "Bead Geometry Volume Error: " << (vol_real_beads - vol_geom_beads)/vol_real_beads*100 << "%" << std::endl << std::endl;

    /* std::cout << std::scientific << std::setprecision(4); */
    std::cout << "Column Length: " << (zTop - zBot + prm->inlet + prm->outlet) * prm->MeshScalingFactor << std::endl<<std::endl;
    std::cout << "zTop - zBot: " << (zTop-zBot)*prm->MeshScalingFactor << std::endl;

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
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" <<std::endl;

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
