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
#include <vector>
#include <string>
#include "Parameters.h"

class PackedBed{
    public:
        PackedBed();
        //PackedBed(const PackedBed& packedBed);
        PackedBed(Parameters * prm);
        virtual ~PackedBed();

        std::vector<Bead *> beads;
        void createGeometry();
        void printPacking();
        void mesh(std::string outfile);
        void createBridge(double relativeBridgeRadius, double eps);
        void createGMSHSphere(double x, double y, double z, double r, double lc_surface, double lc_center, std::vector<int> &shells, std::vector<int> &volumes);
        void GMSHGeom(std::string outfile);

        int tCyl;

        double radius_avg = 0.0;
        double radius_max = -1.0;
        double radius_min = 999999;
        double nBeadsMax = 0;

        std::vector<int> tBeads;
        std::vector<int> tBeadCPs;
        std::vector<std::pair<int, int>> dimTagsBeads;
        std::vector<std::pair<int, int>> dimTagsCyl;
        std::vector<std::pair<int, int>> dimTagsBridges;
        std::vector<std::pair<int, int>> dimTagsInterstitial;

        std::vector<std::pair <int, double>> bridgeTagRadiusPairs;

    private:
        Parameters * prm;

        void getBeads(std::string packingFilename);
        void transformBeads();

};

#endif /* PACKEDBED_H */
