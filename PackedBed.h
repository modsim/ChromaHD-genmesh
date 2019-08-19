/*
 * File:PackedBed.h
 * Desc:Header for PackedBed class. Should contain all the beads.
 * Created By: Rao, Jayghosh Subodh
 * Created On: Fri 05 Apr 2019 04:32:40 PM CEST
 *
 */

#ifndef PACKEDBED_H
#define PACKEDBED_H

#include "Bead.h"
#include "Parameters.h"

#include <string>
#include <float.h>

class PackedBed {
    public:
        PackedBed();
        PackedBed(Parameters * prm);
        virtual ~PackedBed();

        void printPacking();

        std::vector<Bead *> beads;
        double xCyl = 0, yCyl = 0, rCyl = 0;
        double x=0, y=0, z=0, r=0;
        double xMax = -DBL_MAX, xMin = DBL_MAX;
        double yMax = -DBL_MAX, yMin = DBL_MAX;
        double zMax = -DBL_MAX, zMin = DBL_MAX;
        double zBot=0, zTop=0;
        double radius_avg=0, radius_min = 0, radius_max = 0;

        double bedLength = 0;
        double vol_bed_cyl = 0;
        double por_real_bed  = 0;
        double por_geom_bed  = 0;

    private:
        int tCyl;

        int nBeadsMax = 0;
        int nBeads = 0;
        int uLimit = 0;
        int lLimit = 0;

        void getBeads(Parameters * prm);
        void transformBeads(Parameters * prm);
};



#endif /* PACKEDBED_H */
