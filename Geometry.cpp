#include "Geometry.h"
#include <gmsh.h>

namespace factory = gmsh::model::occ;
namespace model = gmsh::model;

Geometry::Geometry(Parameters * prm, PackedBed * pb)
{

    if(prm->dryRun)
        return;

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

void Geometry::operate(Parameters * prm)
{
    if (!prm->geomInfile.empty())
        return;

    long bool_start = gmsh::logger::getWallTime();

    // dimTagsMap
    std::vector<std::vector<std::pair<int, int> > > ovv;

    std::cout << "Intersecting Volumes... " << std::flush;
    //TODO: How would this react with multiple things in containers?
    factory::intersect(dimTagsBeads, dimTagsContainers, dimTagsBeadsInside, ovv, -1, true, false);
    std::cout << "done!" << std::endl;

    /* Fuse beads together (with bridges or without) */
    if (prm->booleanOperation == 1)
    {
        if (dimTagsBridges.size() == 0) dimTagsBridges = {3, dimTagsBeads.back()};

        std::cout << "Fusing Beads and Bridges... " << std::flush;
        /* factory::fuse(dimTagsBeads, dimTagsBridges, ov, ovv ); */
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
        /* if (ov.size() == 0) ov = dimTagsBeads; */
        if (dimTagsFused.size() == 0) dimTagsFused = dimTagsBeadsInside;

        std::cout << "Fragmenting Volumes... " << std::flush;
        factory::fragment(dimTagsContainers, dimTagsFused, dimTagsFragmented, ovv );
        /* dimTagsInterstitial.push_back(bv.back()); */
        std::cout << "done!" << std::endl;
    }

    // Synchronize gmsh model with geometry kernel.
    std::cout << "synchronizing... " << std::flush;
    factory::synchronize();
    std::cout << "done!" << std::endl;

    std::cout << std::endl;
    std::cout << "Number of Geometrical Beads:            " << dimTagsBeadsInside.size() << std::endl;
    std::cout << "Number of Geometrical Bridges:          " << dimTagsBridges.size()     << std::endl;
    std::cout << "Number of Geometrical Internal Volumes: " << dimTagsFused.size()       << std::endl;

    long bool_duration = gmsh::logger::getWallTime() - bool_start;

    gmsh::logger::write("Boolean time: " + std::to_string(bool_duration) + " s", "info");
    std::cout << std::endl;

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
            /* tBeadCPs.push_back(ctag); */
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
            x << " " <<
            y << " " <<
            z << " " <<
            dx << " " <<
            dy << " " <<
            dz << " " <<
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

        std::cout << "Creating cylinder: " <<
            x << " " <<
            y << " " <<
            z << " " <<
            dx << " " <<
            dy << " " <<
            dz << std::endl;

        dimTagsContainers.push_back( {3, factory::addBox(x, y, z, dx, dy, dz)});

        if ((prm->periodic == "xyz") || prm->periodic == "z")
        {
            if (prm->periodicInlet > 0)
                dimTagsContainers.push_back({3, factory::addBox(x, y, z - prm->periodicInlet, dx, dy, prm->periodicInlet)});
            if (prm->periodicOutlet > 0)
                dimTagsContainers.push_back({3, factory::addBox(x, y, z+dz, dx, dy, prm->periodicInlet)});
        }
    }
}
