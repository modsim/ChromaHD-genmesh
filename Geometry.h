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
        std::vector<std::pair<int, int>> dimTagsBeads;
        std::vector<std::pair<int, int>> dimTagsContainers;
        std::vector<std::pair<int, int>> dimTagsBridges;
        std::vector<std::pair<int, int>> dimTagsInterstitial;

        // Store dimTags and dimTagsMap
        std::vector<std::pair<int, int> > dimTagsBeadsInside;
        std::vector<std::pair<int, int> > dimTagsBeadsInPeriodicInlet;
        std::vector<std::pair<int, int> > dimTagsBeadsInPeriodicOutlet;
        std::vector<std::pair<int, int> > dimTagsFused;
        std::vector<std::pair<int, int> > dimTagsFragmented;
        std::vector<std::pair<int, int> > dimTagsFragmentedPeriodicInlet;
        std::vector<std::pair<int, int> > dimTagsFragmentedPeriodicOutlet;
        std::vector<std::pair<int, int> > dimTagsDummy;

        void createContainer(PackedBed * pb, Parameters * prm, std::vector<std::pair<int,int>> &dimTagsContainers);
        void createPackedBed(PackedBed * pb, Parameters * prm, std::vector<std::pair<int, int>> &dimTagsBeads);
        void createBridges(PackedBed * pb, Parameters * prm, std::vector<std::pair<int, int>> &dimTagsBridges);
        void operate(Parameters * prm);

};

#endif /* GEOMETRY_H */
