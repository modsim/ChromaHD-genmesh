#include "Column.h"
#include "Tools.h"

#include <gmsh.h>

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>

namespace model = gmsh::model;
namespace factory = gmsh::model::occ;

Column::Column() {};
Column::~Column() {};

Column::Column(std::vector<std::pair<int,int>> dimTagsFragmentedColumn, Parameters * prm, std::string _periodic)
{
    std::cout << "Creating Column... " << std::endl;

    std::vector<std::pair<int,int>> dimTagsBeadSurface;
    std::vector<std::pair<int,int>> dimTagsOuterSurface;

    periodic = _periodic;
    double xMax, xMin, yMax, yMin, zMax, zMin;

    std::cout << "  > Listing Interstitial Volume... " << std::flush;
    volumes.interstitial.push_back(dimTagsFragmentedColumn.back().second);
    dimTagsFragmentedColumn.pop_back();
    std::cout << "done!" << std::endl;

    // Extract surfaces from interstitial volume
    std::cout << "  > Extracting Surfaces... " << std::flush;
    model::getBoundary({{3,volumes.interstitial[0]}}, dimTagsOuterSurface, false, false, false);
    model::getBoundary(dimTagsFragmentedColumn, dimTagsBeadSurface, false, false, false);
    std::cout << "done!" << std::endl;

    std::cout << "  > Listing Bounding Surfaces... " << std::flush;
    if(prm->containerShape == 0)
    {
        surfaces.walls   = {dimTagsOuterSurface[0].second};
        surfaces.outlet  = {dimTagsOuterSurface[1].second};
        surfaces.inlet   = {dimTagsOuterSurface[2].second};
        // NOTE: What would be the normal for a cylinder wall??
    }
    else if (prm->containerShape == 1)
    {
        separateBoundingSurfaces(dimTagsOuterSurface, outerWalls);
        separateBoundingSurfaces(dimTagsBeadSurface, beadWalls);

        generateBoxSurfaces();
    }
    std::cout << "done!" << std::endl;

    // The remaining entries in bv are bead volumes
    std::cout << "  > Listing Bead Volumes... " << std::flush;
    for (auto it = dimTagsFragmentedColumn.begin(); it != dimTagsFragmentedColumn.end(); it++)
    {
        volumes.beads.push_back((*it).second);
    }
    std::cout << "done!" << std::endl;

    std::cout << "  > Listing Bead Surfaces... " << std::flush;
    for ( std::vector<std::pair< int , int>>::iterator it = dimTagsBeadSurface.begin(); it != dimTagsBeadSurface.end(); it++   )
    {
        surfaces.beads.push_back((*it).second);
    }
    std::cout << "done!" << std::endl;

    // Remove bead surfaces that are also part of the wall
    // from being categorized as bead surfaces.
    // These 'tSBeads' will go on to be doubled in `gmsh2mixdv2 -d 4`
    std::cout << "  > Cleaning Bead Surfaces... " << std::flush;
    subtractTags(surfaces.beads, beadWalls.xleft);
    subtractTags(surfaces.beads, beadWalls.xright);
    subtractTags(surfaces.beads, beadWalls.yleft);
    subtractTags(surfaces.beads, beadWalls.yright);
    subtractTags(surfaces.beads, beadWalls.zleft); // if not xyz, there should be no beadWalls zleft and zright, so it's fine
    subtractTags(surfaces.beads, beadWalls.zright);
    std::cout << "done!" << std::endl;

    model::getBoundingBox(
            3, volumes.interstitial[0],
            xMin, yMin, zMin,
            xMax, yMax, zMax
            );

    dx = xMax - xMin;
    dy = yMax - yMin;
    dz = zMax - zMin;

    stats(); // print out info, good for debugging on dryRun

    if (periodic == "xyz" || periodic == "xy")
    {
        setupPeriodicSurfaces(outerWalls);
        setupPeriodicSurfaces(beadWalls);
    }


}

void Column::AddPhysicalGroups()
{
    model::removePhysicalGroups();

    std::cout << "  > Adding Physical Groups... " << std::flush;
    model::addPhysicalGroup(2, surfaces.inlet      , 1 );
    model::addPhysicalGroup(2, surfaces.outlet     , 2 );
    model::addPhysicalGroup(2, surfaces.walls      , 3 );
    model::addPhysicalGroup(2, surfaces.beads      , 4 );
    model::addPhysicalGroup(3, volumes.interstitial, 5 );
    model::addPhysicalGroup(3, volumes.beads       , 6 );
    std::cout << "done!" << std::endl;

    std::cout << "  > Setting Physical Names... " << std::flush;
    model::setPhysicalName(2,1,"inlet");
    model::setPhysicalName(2,2,"outlet");
    model::setPhysicalName(2,3,"wall");
    model::setPhysicalName(2,4,"beadSurface");
    model::setPhysicalName(3,5,"interstitialVolume");
    model::setPhysicalName(3,6,"beadVolume");
    std::cout << "done!" << std::endl;

}


void Column::separateBoundingSurfaces(std::vector<std::pair<int, int>> dimTagsSurfaces, Walls& tWall)
{

    std::vector<double> points, parametricCoord, curvatures, normals;
    std::vector<std::pair <int, int>> dimTagsBound;

    for (std::vector<std::pair<int, int>>::iterator iter = dimTagsSurfaces.begin(); iter != dimTagsSurfaces.end(); iter++)
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
                tWall.xleft.push_back((*iter).second);
            else if (int(normals[1]) == -1)
                tWall.yleft.push_back((*iter).second);
            else if (int(normals[2]) == -1)
                tWall.zleft.push_back((*iter).second);
            else if (int(normals[0]) == 1)
                tWall.xright.push_back((*iter).second);
            else if (int(normals[1]) == 1)
                tWall.yright.push_back((*iter).second);
            else if (int(normals[2]) == 1)
                tWall.zright.push_back((*iter).second);

        }

    }

}


void Column::setupPeriodicSurfaces(Walls& tWall)
{

    std::vector<double> affineTranslationX = {1, 0, 0, dx, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};
    std::vector<double> affineTranslationY = {1, 0, 0, 0, 0, 1, 0, dy, 0, 0, 1, 0, 0, 0, 0, 1};
    std::vector<double> affineTranslationZ = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, dz, 0, 0, 0, 1};

    std::size_t found;
    found = periodic.find('x');
    if (found != std::string::npos)
    {
        assert(tWall.xleft.size()  == tWall.xright.size());
        matchPeriodicSurfaces(tWall.xleft, tWall.xright, 0, affineTranslationX);
    }

    found = periodic.find('y');
    if (found != std::string::npos)
    {
        assert(tWall.yleft.size()  == tWall.yright.size());
        matchPeriodicSurfaces(tWall.yleft, tWall.yright, 1, affineTranslationY);
    }

    found = periodic.find('z');
    if (found != std::string::npos)
    {
        assert(tWall.zleft.size()  == tWall.zright.size());
        matchPeriodicSurfaces(tWall.zleft, tWall.zright, 2, affineTranslationZ);
    }

}


// TODO: Might need to move to tools in case i have to link outside
void Column::matchPeriodicSurfaces(std::vector<int>& ltags, std::vector<int>& rtags, int per_dir, std::vector<double> affineTranslation)
{
    /* model::mesh::setPeriodic(2, tXRightWallBead, tXLeftWallBead, affineTranslationX); */
    // Unfortunately the order of the surfaces in the vectors matters when providing them in groups to the setPeriodic function.
    // So we match the surfaces individually instead in this function.
    // For each surface, find center of bbox. Project that along per_dir and find matching centers on either surface.
    // If match is found, setPeriodic

    std::cout << "  > Matching periodic surfaces in " << per_dir << " direction... " << std::endl;

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

void Column::generateBoxSurfaces()
{
    surfaces.walls.reserve(
            surfaces.walls.size() +
            outerWalls.xleft.size() +
            outerWalls.yleft.size() +
            outerWalls.xright.size() +
            outerWalls.yright.size() +
            beadWalls.xleft.size() +
            beadWalls.yleft.size() +
            beadWalls.xright.size() +
            beadWalls.yright.size()
            );

    surfaces.walls.insert( surfaces.walls.end(), outerWalls.xleft.begin() , outerWalls.xleft.end()  );
    surfaces.walls.insert( surfaces.walls.end(), outerWalls.xright.begin(), outerWalls.xright.end() );
    surfaces.walls.insert( surfaces.walls.end(), outerWalls.yleft.begin() , outerWalls.yleft.end()  );
    surfaces.walls.insert( surfaces.walls.end(), outerWalls.yright.begin(), outerWalls.yright.end() );

    surfaces.walls.insert( surfaces.walls.end(), beadWalls.xleft.begin() , beadWalls.xleft.end()  );
    surfaces.walls.insert( surfaces.walls.end(), beadWalls.xright.begin(), beadWalls.xright.end() );
    surfaces.walls.insert( surfaces.walls.end(), beadWalls.yleft.begin() , beadWalls.yleft.end()  );
    surfaces.walls.insert( surfaces.walls.end(), beadWalls.yright.begin(), beadWalls.yright.end() );

    // TODO: In z-periodic cases, do I want to include the cut bead surfaces as inlet/outlet?
    // Also: How do I convert data from outlet -> inlet and provide it to xns.in?
    // so far, rngdexp... might need to check the option to specify node-by-node dof constraint

    surfaces.inlet = outerWalls.zleft;
    surfaces.outlet = outerWalls.zright;

}

void Column::write(std::string outpath, std::string outfile, std::string extension)
{
    AddPhysicalGroups();
    gmsh::write(outpath + "/" + outfile + extension);
}


void Column::writeFragments(std::string outpath, std::string outfile, std::string extension)
{
    //TODO: modularize

    model::removePhysicalGroups();

    model::addPhysicalGroup(3,volumes.interstitial,15);
    model::setPhysicalName(3,15,"interstitialVolume");
    gmsh::write(outpath + "/" + outfile + "_volumes_interstitial" + extension);
    model::removePhysicalGroups();

    model::addPhysicalGroup(2, surfaces.inlet  , 11      );
    model::addPhysicalGroup(2, surfaces.outlet , 12      );
    model::addPhysicalGroup(2, surfaces.walls   , 13      );
    model::setPhysicalName (2, 11       , "inlet" );
    model::setPhysicalName (2, 12       , "outlet");
    model::setPhysicalName (2, 13       , "wall"  );
    gmsh::write(outpath + "/" + outfile + "_surfaces_interstitial" + extension);
    model::removePhysicalGroups();

    model::addPhysicalGroup(2, surfaces.inlet  , 11      );
    model::setPhysicalName (2, 11       , "inlet" );
    gmsh::write(outpath + "/" + outfile + "_surfaces_inlet" + extension);
    model::removePhysicalGroups();

    model::addPhysicalGroup(2, surfaces.outlet , 12      );
    model::setPhysicalName (2, 12       , "outlet");
    gmsh::write(outpath + "/" + outfile + "_surfaces_outlet" + extension);
    model::removePhysicalGroups();

    model::addPhysicalGroup(2, surfaces.walls   , 13      );
    model::setPhysicalName (2, 13       , "wall"  );
    gmsh::write(outpath + "/" + outfile + "_surfaces_wall" + extension);
    model::removePhysicalGroups();


    model::addPhysicalGroup(3,volumes.beads,16 );
    model::setPhysicalName(3,16,"beadVolume");
    gmsh::write(outpath + "/" + outfile + "_volumes_beads" + extension);
    model::removePhysicalGroups();

    model::addPhysicalGroup(2,surfaces.beads,14);
    model::setPhysicalName(2,14, "beadSurface");
    gmsh::write(outpath + "/" + outfile + "_surfaces_beads" + extension);

}

void Column::meshVolumes(Parameters * prm)
{
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

}

// Links current Column to a second Column that's right ahead of it in Z
void Column::linkPeriodicZ(Column& second)
{
    assert(outerWalls.zright.size()  == second.outerWalls.zleft.size());
    assert(beadWalls.zright.size()  == second.beadWalls.zleft.size());

    double dz = 0;
    std::vector<double> affineTranslationZ = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, dz, 0, 0, 0, 1};

    matchPeriodicSurfaces(outerWalls.zright, second.outerWalls.zleft, 2, affineTranslationZ);
    matchPeriodicSurfaces(beadWalls.zright, second.beadWalls.zleft, 2, affineTranslationZ);
}

void Column::stats()
{
    std::cout << "========== Column stats ==========" << std::endl
        << "Periodic: " << periodic << std::endl
        << "----- Volumes -----" << std::endl
        << "      beads       : " << volumes.beads.size() << std       ::endl
        << "      interstitial: " << volumes.interstitial.size() << std::endl
        << "----- Surfaces -----" << std::endl
        << "      beads : " << surfaces.beads.size() << std ::endl
        << "      inlet : " << surfaces.inlet.size() << std ::endl
        << "      outlet: " << surfaces.outlet.size() << std::endl
        << "      walls : " << surfaces.walls.size() << std ::endl
        << "----- outerWalls -----" << std::endl
        << "      x: " << outerWalls.xleft.size() << "  " << outerWalls.xright.size() << std::endl
        << "      y: " << outerWalls.yleft.size() << "  " << outerWalls.yright.size() << std::endl
        << "      z: " << outerWalls.zleft.size() << "  " << outerWalls.zright.size() << std::endl
        << "----- beadWalls -----" << std::endl
        << "      x: " << beadWalls.xleft.size() << "  " << beadWalls.xright.size() << std::endl
        << "      y: " << beadWalls.yleft.size() << "  " << beadWalls.yright.size() << std::endl
        << "      z: " << beadWalls.zleft.size() << "  " << beadWalls.zright.size() << std::endl;
    std::cout << "==================================" << std::endl;
}
