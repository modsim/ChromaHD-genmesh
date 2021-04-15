#include "Geometry.h"
#include "Tools.h"
#include <gmsh.h>

namespace factory = gmsh::model::occ;
namespace model = gmsh::model;

Geometry::Geometry(Parameters * prm, PackedBed * pb)
{
    createPackedBed(pb, prm, dimTagsBeads);
    createBridges(pb, prm, dimTagsBridges);
    createContainer(pb, prm, dimTagsContainers);
    operate(prm);

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
        model::getEntities(dimTagsDummy, 0);
        model::mesh::setSize(dimTagsDummy, prm->lc_beads);
        std:: cout << "done!" << std::endl;
    }


}

Geometry::~Geometry()
{
}

// Create packed bed geometry. Basically all the beads.
// Store dimTags in last argument
void Geometry::createPackedBed(PackedBed * pb, Parameters * prm, std::vector<std::pair<int, int>> &dimTagsBeads)
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

            dimTagsBeads.push_back({3, tag});
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


void Geometry::createBridges(PackedBed * pb, Parameters * prm, std::vector<std::pair<int, int>> &dimTagsBridges)
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

}



// Creates a container (cylinder or box) for the specified packed bed
// either automatically or using the specified "cyl" or "box" parameters
// Adds the dimtag of the new geometry to the last argument.
// Generating the additional container for the inlet or outlet periodic special case should be done outside

void Geometry::createContainer(PackedBed * pb, Parameters * prm, std::vector<std::pair<int,int>> &dimTagsContainers)
{

    /* double x = 0, dx = 0; */
    /* double y = 0, dy = 0; */
    /* double z = 0, dz = 0; */
    /* double r = 0; */

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

        dimTagsContainers.push_back({3, factory::addCylinder(x, y, z, dx, dy, dz, R)});

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

        dimTagsContainers.push_back( {3, factory::addBox(x, y, z, dx, dy, dz)});

        // NOTE: You can have periodicInlet and Outlet sections without ANY periodicity.
        // PeriodicInlet and Outlet only necessitate periodic linking which will happen no matter what
        //TODO: Rename these? Make the first ones a plane instead?
        if (prm->periodicInlet > 0)
        {
            dimTagsContainers.push_back({3, factory::addBox(x - dx, y - dy, z - prm->periodicInlet, 3*dx, 3*dy, prm->periodicInlet)});
            dimTagsContainers.push_back({3, factory::addBox(x, y, z - prm->periodicInlet, dx, dy, prm->periodicInlet)});
        }
        if (prm->periodicOutlet > 0)
        {
            dimTagsContainers.push_back({3, factory::addBox(x - dx, y - dy, z+dz, 3*dx, 3*dy, prm->periodicOutlet)});
            dimTagsContainers.push_back({3, factory::addBox(x, y, z+dz, dx, dy, prm->periodicInlet)});
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
    /* factory::intersect(dimTagsBeads, dimTagsContainers, dimTagsBeadsInside, ovv, -1, true, false); */

    if (dimTagsContainers.size() == 1)
        factory::intersect(dimTagsBeads, dimTagsContainers, dimTagsBeadsInside, ovv, -1, true, false);
    else
    {
        std::cout << "  > Preparing Column Geometries... " << std::flush;;
        factory::intersect(dimTagsBeads, {dimTagsContainers[0]}, dimTagsBeadsInside, ovv, -1, false, false);
        std::cout << "done!" << std::endl;;

        //assuming both periodicInlet and periodicOutlet are given!!!
        std::vector<std::pair<int,int>> dimTagsBeadsInlet, dimTagsBeadsOutlet, dtdummy;

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
        factory::intersect(dimTagsBeads, {dimTagsContainers[1]}, dimTagsBeadsInPeriodicInlet, ovv, -1, false, false);
        findIfSurfaceWithNormal(dimTagsBeadsInPeriodicInlet, {0,0,1}, dimTagsBeadsInlet);
        factory::intersect(dimTagsBeadsInlet, {dimTagsContainers[2]}, dimTagsBeadsInPeriodicInlet, ovv, -1, false, false);
        std::cout << "done!" << std::endl;

        std::cout << "  > Preparing Outlet Geometries... " << std::flush;;
        factory::intersect(dimTagsBeads, {dimTagsContainers[3]}, dimTagsBeadsInPeriodicOutlet, ovv, -1, false, false);
        findIfSurfaceWithNormal(dimTagsBeadsInPeriodicOutlet, {0,0,-1}, dimTagsBeadsOutlet);
        factory::intersect(dimTagsBeadsOutlet, {dimTagsContainers[4]}, dimTagsBeadsInPeriodicOutlet, ovv, -1, false, false);
        std::cout << "done!" << std::endl;;

        // TODO: Do I need this both here and after fragmentation?
        std::cout << "  > Cleaning Model... " << std::flush;;
        factory::getEntities(dtdummy);
        subtractDimTags(dtdummy, dimTagsBeadsInside);
        subtractDimTags(dtdummy, dimTagsBeadsInPeriodicInlet);
        subtractDimTags(dtdummy, dimTagsBeadsInPeriodicOutlet);
        subtractDimTags(dtdummy, {dimTagsContainers[0]});
        subtractDimTags(dtdummy, {dimTagsContainers[2]});
        subtractDimTags(dtdummy, {dimTagsContainers[4]});
        factory::remove(dtdummy, true);
        std::cout << "done!" << std::endl;

    }

    // Synchronize gmsh model with geometry kernel.
    std::cout << "Synchronizing... " << std::flush;
    factory::synchronize();
    std::cout << "done!" << std::endl;

    /* Fuse beads together (with bridges or without) */
    if (prm->booleanOperation == 1)
    {
        if (dimTagsBridges.size() == 0) dimTagsBridges = {3, dimTagsBeads.back()};

        std::cout << "Fusing Beads and Bridges... " << std::flush;
        factory::fuse(dimTagsBeadsInside, dimTagsBridges, dimTagsFused, ovv );
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
        factory::cut(dimTagsBeadsInside, dimTagsBridges, dimTagsFused, ovv );
        std::cout << "done!" << std::endl;

    }
    else
    {
        /* ov = dimTagsBeads; */
        dimTagsFused = dimTagsBeadsInside;
    }

    // Fragment cylinder w.r.t. beads
    if (prm->fragment)
    {
        if (dimTagsFused.size() == 0) dimTagsFused = dimTagsBeadsInside;

        std::cout << "Fragmenting Volumes..." << std::endl;
        /* factory::fragment(dimTagsContainers, dimTagsFused, dimTagsFragmented, ovv ); */

        //TODO: don't use dimtagsfused here?
        if (dimTagsContainers.size() == 1)
            factory::fragment(dimTagsContainers, dimTagsFused, dimTagsFragmented, ovv );
        else
        {
            std::cout << "  > Fragmenting Column... " << std::flush;
            factory::fragment({dimTagsContainers[0]}, dimTagsFused, dimTagsFragmented, ovv );
            std::cout << "done!" << std::endl;

            //NOTE: assuming both periodicInlet and periodicOutlet are given!!!
            std::cout << "  > Fragmenting Inlet..." << std::flush;
            factory::fragment({dimTagsContainers[2]}, dimTagsBeadsInPeriodicInlet, dimTagsFragmentedPeriodicInlet, ovv );
            std::cout << "done!" << std::endl;
            std::cout << "  > Fragmenting Outlet..." << std::flush;
            factory::fragment({dimTagsContainers[4]}, dimTagsBeadsInPeriodicOutlet, dimTagsFragmentedPeriodicOutlet, ovv );
            std::cout << "done!" << std::endl;

            // have to do this because I can't delete the beads directly during/after intersection as with simpler cases.
            // But there might be simpler ways of achieving the same thing: I won't need intersect if I'm just fragmenting and deleting by bounding box
            std::cout << "  > Cleaning Model..." << std::flush;;
            std::vector<std::pair<int,int>> dtdummy;
            factory::getEntities(dtdummy);
            subtractDimTags(dtdummy, dimTagsFragmented);
            subtractDimTags(dtdummy, dimTagsFragmentedPeriodicInlet);
            subtractDimTags(dtdummy, dimTagsFragmentedPeriodicOutlet);
            factory::remove(dtdummy,true);
            std::cout << "done!" << std::endl;

        }

    }

    // Synchronize gmsh model with geometry kernel.
    std::cout << "Synchronizing model... " << std::flush;
    factory::synchronize();
    std::cout << "done!" << std::endl;

    std::cout << "Number of Geometrical Beads:            " << dimTagsBeadsInside.size() << std::endl;
    std::cout << "Number of Geometrical Bridges:          " << dimTagsBridges.size()     << std::endl;
    std::cout << "Number of Geometrical Internal Volumes: " << dimTagsFused.size()       << std::endl;

    long bool_duration = gmsh::logger::getWallTime() - bool_start;

    gmsh::logger::write("Boolean time: " + std::to_string(bool_duration) + " s", "info");

}
