#include "Geometry.h"
#include "Tools.h"
#include <gmsh.h>

namespace factory = gmsh::model::occ;
namespace model = gmsh::model;

#include "Data.h"
#include <numeric> // ::accumulate

Geometry::Geometry(Parameters * prm, PackedBed * pb)
{
    createPackedBed(pb, prm, dt_beads);
    createBridges(pb, prm, dt_bridges);
    createContainer(pb, prm, dt_containers);
    operate(prm);

    // New method
    /* generate(prm); */

    // TODO: Move to Model
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
    else if (prm->meshSizeMethod == 0)
    {
        //Set mesh size globally
        std::cout << "Setting global mesh size... ";
        model::getEntities(dt_dummy, 0);
        model::mesh::setSize(dt_dummy, prm->lc_beads);
        std:: cout << "done!" << std::endl;
    }


}

Geometry::~Geometry()
{
}

// Create packed bed geometry. Basically all the beads.
// Store dimTags in last argument
void Geometry::createPackedBed(PackedBed * pb, Parameters * prm, std::vector<std::pair<int, int>> &dt_beads)
{
    std::cout << "Creating beads... " << std::flush;
    int tag, ctag;

    for (std::vector<Bead *>::iterator iter = pb->beads.begin(); iter != pb->beads.end(); iter++ )
    {
        double x = (*iter)->getX();
        double y = (*iter)->getY();
        double z = (*iter)->getZ();
        double r = (*iter)->getR();

        if (prm->geomInfile.empty())
        {
            tag = factory::addSphere(x, y, z, r);
            /* std::cout << tag << std::endl; */
            ctag = factory::addPoint(x, y, z, prm->lc_beads, -1);

            (*iter)->setTag(tag);
            (*iter)->setCTag(ctag);

            dt_beads.push_back({3, tag});
        }

        /*
         * Distance/Threshold allows creating gradient mesh sizes within the beads
         */
        //TODO: Update to new GMSH changes
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


void Geometry::createBridges(PackedBed * pb, Parameters * prm, std::vector<std::pair<int, int>> &dt_bridges)
{
    int tagBridge;

    double dx, dy, dz;
    double x4, y4, z4;
    double x3, y3, z3;
    double x2, y2, z2, r2;
    double x1, y1, z1, r1;
    double dx3, dy3, dz3;
    double factor = prm->bridgeOffsetRatio;
    double dist;

    std::vector<Bead *> remBeads = pb->beads;
    std::cout << "Creating Bridges... " << std::flush;
    for(std::vector<Bead *>::iterator iter = pb->beads.begin(); iter != pb->beads.end(); iter++)
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

                    this->dt_bridges.push_back({3, tagBridge }) ;
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

}



// Creates a container (cylinder or box) for the specified packed bed
// either automatically or using the specified "cyl" or "box" parameters
// Adds the dimtag of the new geometry to the last argument.
// Generating the additional container for the inlet or outlet periodic special case should be done outside

void Geometry::createContainer(PackedBed * pb, Parameters * prm, std::vector<std::pair<int,int>> &dt_containers)
{

    if (prm->containerShape == 0)
    {
        if(prm->autoContainment == 1)
        {
            x = pb->xCyl;         dx = 0;
            y = pb->yCyl;         dy = 0;
            z = pb->zCylBot;      dz = pb->zCylTop - pb->zCylBot; // includes inlet/outlet
            R = pb->rCyl;       // should include delta
        }
        else
        {
            // Strictly specified dimensions.
            x = prm->x0;         dx = prm->dx;
            y = prm->y0;         dy = prm->dy;
            z = prm->z0;         dz = prm->dz;
            R = prm->rCyl;
        }

        std::cout << "Creating cylinder: " <<
            x << " " << y << " " << z << " " <<
            dx << " " << dy << " " << dz << " " <<
            R << std::endl;;

        dt_containers.push_back({3, factory::addCylinder(x, y, z, dx, dy, dz, R)});

        // TODO: if periodic == "z"
    }
    else if (prm->containerShape == 1)
    {
        if(prm->autoContainment == 1)
        {
            x = pb->xMin;         dx = pb->xMax - pb->xMin;
            y = pb->yMin;         dy = pb->yMax - pb->yMin;
            z = pb->zCylBot;      dz = pb->zCylTop - pb->zCylBot; // includes inlet/outlet?

            // NOTE: Not acccounting for rcyldelta here, as opposed to cylindrical containers.
            // TODO: Figure this out
        }
        else
        {
            // Strictly specified dimensions.
            x = prm->x0;         dx = prm->dx;
            y = prm->y0;         dy = prm->dy;
            z = prm->z0;         dz = prm->dz;
        }

        std::cout << std::setprecision(2);
        std::cout << "Creating box: " <<
            x << " " << y << " " << z << " " <<
            dx << " " << dy << " " << dz << std::endl;
        std::cout << std::setprecision(10);

        dt_containers.push_back( {3, factory::addBox(x, y, z, dx, dy, dz)});

        // NOTE: You can have periodicInlet and Outlet sections without ANY periodicity.
        // PeriodicInlet and Outlet only necessitate periodic linking which will happen no matter what
        //TODO: Rename these? Make the first ones a plane instead?
        if (prm->periodicInlet > 0)
        {
            dt_containers.push_back({3, factory::addBox(x - dx, y - dy, z - prm->periodicInlet, 3*dx, 3*dy, prm->periodicInlet)});
            dt_containers.push_back({3, factory::addBox(x, y, z - prm->periodicInlet, dx, dy, prm->periodicInlet)});
        }
        if (prm->periodicOutlet > 0)
        {
            dt_containers.push_back({3, factory::addBox(x - dx, y - dy, z+dz, 3*dx, 3*dy, prm->periodicOutlet)});
            dt_containers.push_back({3, factory::addBox(x, y, z+dz, dx, dy, prm->periodicInlet)});
        }

    }
}

void Geometry::operate(Parameters * prm)
{
    if (!prm->geomInfile.empty())
        return;

    long bool_start = gmsh::logger::getWallTime();

    // dimTagsMap
    std::vector<std::vector<std::pair<int, int> > > ovv;

    std::cout << "Intersecting Volumes... " << std::endl;

    if (dt_containers.size() == 1)
        factory::intersect(dt_beads, dt_containers, dt_beadsInside, ovv, -1, true, false);
    else
    {
        std::cout << "  > Preparing Column Geometries... " << std::flush;;
        factory::intersect(dt_beads, {dt_containers[0]}, dt_beadsInside, ovv, -1, false, false);
        std::cout << "done!" << std::endl;;

        //assuming both periodicInlet and periodicOutlet are given!!!
        std::vector<std::pair<int,int>> dt_beadsInlet, dt_beadsOutlet, dt_dummy;

        // NOTE:
        // Possible fixes to this mess:
        // 1. Intersect with plane, extract dimtags out of ovv
        // 2. getEntities in Bounding Box after intersect? Not enough.. need to fragment cleanly.
        // 3. Move to a fragment only based workflow. Fragment, Remove Outside bounding box.
            // if so, I might not even have to do all this manual handling. All Containers X All beads.
            // Still leaves the issue of ensuring xy periodicity in inlet-outlet
            //

        // Periodic Inlet and Outlet geometries are tricky to properly capture when dealing with full periodicity.
        // Things to note:
        //      1. We need to capture beads that are within the (in/out) volume, but not all of them
        //      2. With XYZ periodicity, Z-direction stacking is necessary, filling the whole volume with unnecessary beads.
        //      3. In case of XYZ periodicity, Inlet and Outlet will have XY periodicity. If so, finding the beads to keep in the
        //          volume becomes harder. The criterion is not just intersection of the bead with the z0 plane of the column; this
        //          doesn't account for XY projections of beads. Instead, we find intersection of beads with the full z0 plane, and then
        //          intersect again to keep them within the specified volume.


        std::cout << "  > Preparing Inlet Geometries... " << std::flush;;
        factory::intersect(dt_beads, {dt_containers[1]}, dt_beadsInPeriodicInlet, ovv, -1, false, false);
        filterIfSurfaceWithNormal(dt_beadsInPeriodicInlet, {0,0,1}, dt_beadsInlet);
        factory::intersect(dt_beadsInlet, {dt_containers[2]}, dt_beadsInPeriodicInlet, ovv, -1, false, false);
        std::cout << "done!" << std::endl;

        std::cout << "  > Preparing Outlet Geometries... " << std::flush;;
        factory::intersect(dt_beads, {dt_containers[3]}, dt_beadsInPeriodicOutlet, ovv, -1, false, false);
        filterIfSurfaceWithNormal(dt_beadsInPeriodicOutlet, {0,0,-1}, dt_beadsOutlet);
        factory::intersect(dt_beadsOutlet, {dt_containers[4]}, dt_beadsInPeriodicOutlet, ovv, -1, false, false);
        std::cout << "done!" << std::endl;;

        // TODO: Do I need this both here and after fragmentation?
        // TODO: Parallellize
        std::cout << "  > Cleaning Model... " << std::flush;;
        factory::getEntities(dt_dummy);
        subtractDimTags(dt_dummy, dt_beadsInside);
        subtractDimTags(dt_dummy, dt_beadsInPeriodicInlet);
        subtractDimTags(dt_dummy, dt_beadsInPeriodicOutlet);
        subtractDimTags(dt_dummy, {dt_containers[0]});
        subtractDimTags(dt_dummy, {dt_containers[2]});
        subtractDimTags(dt_dummy, {dt_containers[4]});
        factory::remove(dt_dummy, true);
        std::cout << "done!" << std::endl;

    }

    // Synchronize gmsh model with geometry kernel.
    std::cout << "Synchronizing... " << std::flush;
    factory::synchronize();
    std::cout << "done!" << std::endl;

    /* Fuse beads together (with bridges or without) */
    if (prm->booleanOperation == 1)
    {
        if (dt_bridges.size() == 0) dt_bridges = {3, dt_beads.back()};

        std::cout << "Fusing Beads and Bridges... " << std::flush;
        factory::fuse(dt_beadsInside, dt_bridges, dt_fused, ovv );
        std::cout << "done!" << std::endl;
    }
    else if (prm->booleanOperation == 2)
    {
        if (dt_bridges.size() == 0)
        {
            std::cout << "Error: No Bridges to cut from beads" << std::endl;
            exit(-1);
        }

        std::cout << "Capping Beads... " << std::flush;
        /* factory::cut(dimTagsBeads, dimTagsBridges, ov, ovv ); */
        factory::cut(dt_beadsInside, dt_bridges, dt_fused, ovv );
        std::cout << "done!" << std::endl;

    }
    else
    {
        /* ov = dimTagsBeads; */
        dt_fused = dt_beadsInside;
    }

    // Fragment cylinder w.r.t. beads
    if (prm->fragment)
    {
        if (dt_fused.size() == 0) dt_fused = dt_beadsInside;

        std::cout << "Fragmenting Volumes..." << std::endl;
        /* factory::fragment(dimTagsContainers, dt_fused, dimTagsFragmented, ovv ); */

        //TODO: don't use dimtagsfused here?
        if (dt_containers.size() == 1)
            factory::fragment(dt_containers, dt_fused, dt_fragmented, ovv );
        else
        {
            std::cout << "  > Fragmenting Column... " << std::flush;
            factory::fragment({dt_containers[0]}, dt_fused, dt_fragmented, ovv );
            std::cout << "done!" << std::endl;

            //NOTE: assuming both periodicInlet and periodicOutlet are given!!!
            std::cout << "  > Fragmenting Inlet..." << std::flush;
            factory::fragment({dt_containers[2]}, dt_beadsInPeriodicInlet, dt_fragmentedPeriodicInlet, ovv );
            std::cout << "done!" << std::endl;
            std::cout << "  > Fragmenting Outlet..." << std::flush;
            factory::fragment({dt_containers[4]}, dt_beadsInPeriodicOutlet, dt_fragmentedPeriodicOutlet, ovv );
            std::cout << "done!" << std::endl;

            // have to do this because I can't delete the beads directly during/after intersection as with simpler cases.
            // But there might be simpler ways of achieving the same thing: I won't need intersect if I'm just fragmenting and deleting by bounding box
            std::cout << "  > Cleaning Model..." << std::flush;;
            std::vector<std::pair<int,int>> dtdummy;
            factory::getEntities(dtdummy);
            subtractDimTags(dtdummy, dt_fragmented);
            subtractDimTags(dtdummy, dt_fragmentedPeriodicInlet);
            subtractDimTags(dtdummy, dt_fragmentedPeriodicOutlet);
            factory::remove(dtdummy,true);
            std::cout << "done!" << std::endl;

        }

    }

    // Synchronize gmsh model with geometry kernel.
    std::cout << "Synchronizing model... " << std::flush;
    factory::synchronize();
    std::cout << "done!" << std::endl;

    std::cout << "Number of Geometrical Beads:            " << dt_beadsInside.size() << std::endl;
    std::cout << "Number of Geometrical Bridges:          " << dt_bridges.size()     << std::endl;
    std::cout << "Number of Geometrical Internal Volumes: " << dt_fused.size()       << std::endl;

    long bool_duration = gmsh::logger::getWallTime() - bool_start;

    gmsh::logger::write("Boolean time: " + std::to_string(bool_duration) + " s", "info");

}

void Geometry::generate(Parameters * prm)
{
    if (!prm->geomInfile.empty())
        return;

    factory::synchronize();

    std::vector<std::pair<int,int>> dt_containerSurfaces, dt_beadSurfaces;
    std::vector<int> t_containerSurfaces, t_beadSurfaces, t_shells;
    std::vector<int> t_beadVolumes;
    int t_interstitialVolume;

    model::getBoundary(dt_containers, dt_containerSurfaces);
    model::getBoundary(dt_beads, dt_beadSurfaces);

    //TODO: Check bounding boxes for collisions and other stuff

    if (prm->periodic == "off")
    {
        extractTags(dt_containerSurfaces, t_containerSurfaces);
        extractTags(dt_beadSurfaces, t_beadSurfaces);

        t_shells.push_back(model::geo::addSurfaceLoop(t_containerSurfaces));
        for (auto it : t_beadSurfaces)
        {
            int t_beadShell = model::geo::addSurfaceLoop({it});
            t_shells.push_back(t_beadShell);
        }

        t_interstitialVolume = model::geo::addVolume(t_shells);

        t_shells.erase(t_shells.begin());

        for (auto it : t_shells)
            t_beadVolumes.push_back(model::geo::addVolume({it}));

        factory::remove(dt_containers);
        factory::remove(dt_beads);

        factory::synchronize();
        model::geo::synchronize();

        // ------------

        Surfaces surfaces;
        Volumes volumes;

        std::vector<std::pair<int,int>> dt_points;
        std::vector<double> points, parametricCoord, curvatures, normals;

        std::vector<double> nzleft = {0, 0, -1};
        std::vector<double> nzright= {0, 0, 1};

        for (auto iSurface: dt_containerSurfaces)
        {
            model::getBoundary({iSurface}, dt_points, false, false, true);
            for (auto iPoint: dt_points)
            {
                model::getValue(iPoint.first, iPoint.second, {}, points);
                model::getParametrization(iSurface.first, iSurface.second, points, parametricCoord);
                model::getCurvature(iSurface.first, iSurface.second, parametricCoord, curvatures);

                if (std::accumulate(curvatures.begin(), curvatures.end(), 0) > 0)
                {
                    surfaces.walls.push_back(iSurface.second);
                    break;
                }

                model::getNormal(iSurface.second, parametricCoord, normals);

                if (normals == nzleft)
                {
                    surfaces.inlet.push_back(iSurface.second);
                    break;
                }
                else if (normals == nzright)
                {
                    surfaces.outlet.push_back(iSurface.second);
                    break;
                }
            }
        }

        surfaces.beads = t_beadSurfaces;
        volumes.beads = t_beadVolumes;
        volumes.interstitial.push_back(t_interstitialVolume);

        model::addPhysicalGroup(2, surfaces.inlet      , 1 );
        model::addPhysicalGroup(2, surfaces.outlet     , 2 );
        model::addPhysicalGroup(2, surfaces.walls      , 3 );
        model::addPhysicalGroup(2, surfaces.beads      , 4 );
        model::addPhysicalGroup(3, volumes.interstitial, 5 );
        model::addPhysicalGroup(3, volumes.beads       , 6 );

        model::setPhysicalName(2,1,"inlet");
        model::setPhysicalName(2,2,"outlet");
        model::setPhysicalName(2,3,"wall");
        model::setPhysicalName(2,4,"beadSurface");
        model::setPhysicalName(3,5,"interstitialVolume");
        model::setPhysicalName(3,6,"beadVolume");

        model::mesh::generate(3);

        gmsh::write("output_new.vtk");

        exit(0);

    }
    else
    {

        std::vector<std::vector<std::pair<int,int>>> dtm_fragmented;
        std::vector<std::pair<int,int>> dt_fragmented, dt_EntitiesInBox;
        std::vector<int> t_EntitiesInBox;

        double xmin, ymin, zmin;
        double xmax, ymax, zmax;

        extractTags(dt_containerSurfaces, t_containerSurfaces);
        extractTags(dt_beadSurfaces, t_beadSurfaces);

        factory::fragment(dt_containerSurfaces, dt_beadSurfaces, dt_fragmented, dtm_fragmented);

        factory::getBoundingBox(3, dt_containers[0].second, xmin, ymin, zmin, xmax,ymax,zmax);

        double eps = 1e-3;
        Surfaces surfaces;
        Volumes volumes;

        factory::getEntitiesInBoundingBox(xmin-eps, ymin-eps, zmin-eps, xmax+eps, ymax+eps, zmax+eps, dt_EntitiesInBox, 2);

        std::vector<std::pair<int,int>> dt_points;
        std::vector<double> points, parametricCoord, curvatures, normals;

        std::vector<double> nzleft = {0, 0, -1};
        std::vector<double> nzright= {0, 0, 1};

        for (auto iSurface : dt_EntitiesInBox)
        {
            model::getBoundary({iSurface}, dt_points, false, false, true);
            for (auto iPoint: dt_points)
            {
                model::getValue(iPoint.first, iPoint.second, {}, points);
                model::getParametrization(iSurface.first, iSurface.second, points, parametricCoord);
                model::getCurvature(iSurface.first, iSurface.second, parametricCoord, curvatures);

                if (std::accumulate(curvatures.begin(), curvatures.end(), 0) > 0)
                {
                    surfaces.beads.push_back(iSurface.second);
                    // TODO:
                    // if complete sphere: continue. How to check? (Positive curvature at all points?)
                        // gettype, getprinciplecurvature, parameterization, parametrizationbounds, derivative, secondderivative, classifysurfaces, computehomology
                        // getadjacencies
                    // if not, put it in list(cut-bead-inner) (later use this for surface loop)
                    break;
                }

                model::getNormal(iSurface.second, parametricCoord, normals);

                // TODO:
                // If bbox is same as any wall bbox: list(interstitial-bound)
                // If bbox is smaller: Add to list(cut-bead-walls)

                if      (normals == nzleft)  { surfaces.inlet.push_back(iSurface.second); break; }
                else if (normals == nzright) { surfaces.outlet.push_back(iSurface.second); break; }

            }

        }

        //TODO:
        // loopOuter = cut-bead-walls + interstitial-bounds
        // wholeBeads = whole beads surfaces
        // cut beads = cut beads-inner + cut beads-walls //NOTE: How to match these?, just a fuse?, assuming you get closed indiv. surfaces, each can be surflooped


    }

}
