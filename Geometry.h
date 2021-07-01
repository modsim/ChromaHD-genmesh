/*
 * File:Geometry.h
 * Desc:Header for Geometry class.
 * Created By: Rao, Jayghosh Subodh
 * Created On:
 *
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "Parameters.h"
#include "PackedBed.h"

class Geometry{
    public:
        Geometry();
        Geometry(Parameters * prm, PackedBed * pb);
        virtual ~Geometry();

        int count = 0;
        std::vector<double> vCount;

        double x, y, z, dx, dy, dz, R;

        std::vector<int> tBeads;
        std::vector<int> tBeadCPs;
        std::vector<std::pair<int, int>> dt_beads;
        std::vector<std::pair<int, int>> dt_containers;
        std::vector<std::pair<int, int>> dt_bridges;

        // Store dimTags and dimTagsMap
        std::vector<std::pair<int, int> > dt_beadsInside;
        std::vector<std::pair<int, int> > dt_beadsInPeriodicInlet;
        std::vector<std::pair<int, int> > dt_beadsInPeriodicOutlet;
        std::vector<std::pair<int, int> > dt_fused;
        std::vector<std::pair<int, int> > dt_fragmented;
        std::vector<std::pair<int, int> > dt_fragmentedPeriodicInlet;
        std::vector<std::pair<int, int> > dt_fragmentedPeriodicOutlet;
        std::vector<std::pair<int, int> > dt_dummy;

        void createContainer(PackedBed * pb, Parameters * prm, std::vector<std::pair<int,int>> &dt_containers);
        void createPackedBed(PackedBed * pb, Parameters * prm, std::vector<std::pair<int, int>> &dt_beads);
        void createBridges(PackedBed * pb, Parameters * prm, std::vector<std::pair<int, int>> &dt_bridges);
        void operate(Parameters * prm);

};

#endif /* GEOMETRY_H */
