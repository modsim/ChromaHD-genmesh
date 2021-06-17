#include "Tools.h"

/* #include <iterator> */
#include <algorithm>
#include <iostream>
#include <gmsh.h>
#include <unordered_set>

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

void subtractTags(std::vector<int>& mainVec, std::vector<int>& subVec)
{
    /*
     * Remove bead surfaces that are also part of the wall
     * from being categorized as bead surfaces.
     * These 'tSBeads' will go on to be doubled in `gmsh2mixdv2 -d 4`
     */
    for (auto it : subVec)
    {
        auto it2 = std::find(mainVec.begin(), mainVec.end(), it);
        if (it2 != mainVec.end()) mainVec.erase(it2);
    }

}

void addTags(std::vector<int>& mainVec, const std::vector<int> subVec)
{
    // Add subvec to mainvec
    mainVec.reserve(mainVec.size() + subVec.size());
    mainVec.insert(mainVec.end(), subVec.begin(), subVec.end());
}

void printDimTags(std::vector<std::pair<int,int>> dimTags)
{
    for (auto it : dimTags)
        printf("%d\t%4d\n", it.first, it.second);
}

void subtractDimTags(std::vector<std::pair<int,int>>& mainVec, std::vector<std::pair<int,int>> subVec)
{
    /* std::cout << "Subtracting dimtags from " << mainVec.size(); */
    for (auto it : subVec)
    {
        auto it2 = std::find(mainVec.begin(), mainVec.end(), it);
        if (it2 != mainVec.end()) mainVec.erase(it2);
    }
    /* std::cout << " to " << mainVec.size() << std::endl; */

}

void addDimTags(std::vector<std::pair<int,int>>& mainVec, const std::vector<std::pair<int,int>> subVec)
{
    // Add subvec to mainvec
    /* std::cout << "Adding dimtags from " << mainVec.size(); */
    /* mainVec.reserve(mainVec.size() + subVec.size()); */
    /* mainVec.insert(mainVec.end(), subVec.begin(), subVec.end()); */
    /* std::cout << " to " << mainVec.size() << std::endl; */

    for (auto it:subVec)
    {
        auto it2 = std::find(mainVec.begin(), mainVec.end(), it);
        if (it2 == mainVec.end()) mainVec.push_back(it);
    }


}

void findSurfacesWithNormal(std::vector<std::pair<int,int>> dimTagsInput, std::vector<double> _normals, std::vector<std::pair<int,int>>& dimTagsOutput)
{
    //Extracts surfaces with a given normal out of the 3D object

    factory::synchronize();
    std::vector<std::pair<int,int>> dimTagsInputBoundaries;
    model::getBoundary(dimTagsInput, dimTagsInputBoundaries, false, false, false);

    for (std::vector<std::pair<int, int>>::iterator iter = dimTagsInputBoundaries.begin(); iter != dimTagsInputBoundaries.end(); iter++)
    {

        // Get the curvatures of the surfaces,
        // Select planar surfaces
        // Get normals of planes to match surface position to tag

        std::vector<std::pair<int,int>> dimTagsBound;
        std::vector<double> points, parametricCoord, curvatures, normals;

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

            // NOTE: equating doubles. Might wanna be wary here
            if (normals == _normals) dimTagsOutput.push_back({(*iter).first, (*iter).second});

        }
    }
}

void findIfSurfaceWithNormal(std::vector<std::pair<int,int>> dimTagsInput, std::vector<double> _normals, std::vector<std::pair<int,int>>& dimTagsOutput)
{
    // NOTE: returns list of 3d objects from input containing surfaces that point in _normals direction.
    // Useful to calculate beads with cut surfaces, for eg.

    factory::synchronize();

    for (auto it:dimTagsInput) // for every bead
    {
        std::vector<std::pair<int,int>> dimTagsBoundaries;

        model::getBoundary({it}, dimTagsBoundaries, false, false, false); // get boundaries of each bead

        for (auto it2:dimTagsBoundaries) // for every boundary surface
        {
            std::vector<std::pair<int,int>> dimTagsBound;
            std::vector<double> points, _points, parametricCoord, curvatures, normals;

            model::getBoundary({it2}, dimTagsBound, false, false, true);
            for (auto it3:dimTagsBound) // for every recursed "point" in the boundary surface representation
            {
                model::getValue(it3.first, it3.second, {}, _points);

                model::getParametrization(it2.first, it2.second, _points, parametricCoord);

                /* model::getCurvature(it2.first, it2.second, parametricCoord, curvatures); */
                /* std::cout << curvatures.size() << ": "; */
                /* for(auto ic: curvatures) */
                /*     std::cout << ic << " "; */
                /* std::cout << std::endl; */

                model::getNormal(it2.second, parametricCoord, normals);

                // if normal matches provided normal, push bead (3d) onto output vector
                if (normals == _normals)
                {
                    dimTagsOutput.push_back(it);
                    // break out of surface loop and check next bead
                    break;
                }
            }
        }
    }
}

void findCutBeads(std::vector<std::pair<int,int>> dimTagsBeads, std::vector<std::pair<int,int>> dimTagsOutput)
{
    // synch
    // for each
        // getsurfaces
        // getcurvature
        // if curvature is not 0, push dimtag
        //

    factory::synchronize();

    /* std::cout << dimTagsBeads.size() << std::endl; */

    for (auto it:dimTagsBeads)
    {
        std::vector<std::pair<int,int>>  dtbound;
        model::getBoundary({it}, dtbound, false, false, false);

        /* std::cout << "  " << dtbound.size() << std::endl; */
        /* std::cout << it.first << "  " << it.second << std::endl; */

        /* double xmax, ymax, zmax, xmin, ymin, zmin; */
        /* model::getBoundingBox(dtbound, xmin, ymin,zmin, xmax,ymax,zmax); */


        if (dtbound.size() > 1)
            dimTagsOutput.push_back(it);

    }

}

void findUncutBeads(std::vector<std::pair<int,int>> dimTagsBeads, std::vector<std::pair<int,int>> dimTagsOutput)
{
    // synch
    // for each
        // getsurfaces
        // getcurvature
        // if curvature is not 0, push dimtag
        //

    factory::synchronize();

    /* std::cout << dimTagsBeads.size() << std::endl; */

    for (auto it:dimTagsBeads)
    {
        std::vector<std::pair<int,int>>  dtbound;
        model::getBoundary({it}, dtbound, false, false, false);

        /* std::cout << "  " << dtbound.size() << std::endl; */
        /* std::cout << it.first << "  " << it.second << std::endl; */

        /* double xmax, ymax, zmax, xmin, ymin, zmin; */
        /* model::getBoundingBox(dtbound, xmin, ymin,zmin, xmax,ymax,zmax); */


        if (dtbound.size() == 1)
            dimTagsOutput.push_back(it);

    }

}


void printDimTagsMap(const std::vector<std::vector<std::pair<int,int>>>& ovv)
{
    int index = 0;
    for (auto iv:ovv)
    {
        index++;
        if (iv.size() == 0 ) continue;
        std::cout <<  index << ": " << iv.size() << ": ";
        for (auto i:iv)
            std::cout << i.second << " ";
        std::cout << std::endl;
    }

}

void extractTags(const std::vector<std::pair<int,int>>& dt, std::vector<int>& tags)
{
    for (auto it : dt)
    {
        tags.push_back(it.second);
    }
}

