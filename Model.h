/*
 * File:Model.h
 * Desc:Header for PackedBed class. Should contain all the beads.
 * Created By: Rao, Jayghosh Subodh
 * Created On: Fri 05 Apr 2019 04:32:40 PM CEST
 *
 */

#ifndef MODEL_H
#define MODEL_H

#include "Bead.h"
#include "Parameters.h"
#include "PackedBed.h"

#include <vector>
#include <string>

class Model{
    public:
        Model();
        Model(Parameters * prm);
        virtual ~Model();

        /* std::vector<Bead *> beads; */
        void createGeometry(PackedBed * packedBed, Parameters * prm);
        /* void printPacking(); */
        void mesh(std::string outfile, Parameters * prm);
        /* void createBridge(double relativeBridgeRadius, double eps); */
        /* void createGMSHSphere(double x, double y, double z, double r, double lc_surface, double lc_center, std::vector<int> &shells, std::vector<int> &volumes); */
        /* void GMSHGeom(std::string outfile); */

        /* int tCyl; */
        /* double radius_avg = 0.0; */
        /* double radius_max = -1.0; */
        /* double radius_min = 999999; */
        /* int nBeadsMax = 0; */
        /* int nBeads = 0; */

        std::vector<int> tBeads;
        std::vector<int> tBeadCPs;
        std::vector<std::pair<int, int>> dimTagsBeads;
        std::vector<std::pair<int, int>> dimTagsCyl;
        std::vector<std::pair<int, int>> dimTagsBridges;
        std::vector<std::pair<int, int>> dimTagsInterstitial;

        std::vector<std::pair <int, double>> bridgeTagRadiusPairs;

    /* private: */
    /*     Parameters * prm; */
    /*     void getBeads(std::string packingFilename); */
    /*     void transformBeads(); */

};

#endif /* MODEL_H */
