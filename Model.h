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
        Model(Parameters * prm, std::string geometryFile);
        virtual ~Model();

        void createGeometry(PackedBed * packedBed, Parameters * prm);
        void mesh(std::string outfile, Parameters * prm);
        void createNamedGroups(std::vector<std::pair<int,int>> ov, int containerShape);

        std::vector<int> tBeads;
        std::vector<int> tBeadCPs;
        std::vector<std::pair<int, int>> dimTagsBeads;
        std::vector<std::pair<int, int>> dimTagsCyl;
        std::vector<std::pair<int, int>> dimTagsBridges;
        std::vector<std::pair<int, int>> dimTagsInterstitial;

        std::vector<int> tVBeads, tVInt, tSBeads, tSWall, tSOutlet, tSInlet;
        std::vector<int> tXLeftWallInt, tYLeftWallInt, tZLeftWallInt, tXRightWallInt, tYRightWallInt, tZRightWallInt;
        /* int tXLeftWallInt, tYLeftWallInt, tZLeftWallInt, tXRightWallInt, tYRightWallInt, tZRightWallInt; */
        std::vector<std::pair <int, double>> bridgeTagRadiusPairs;

};

#endif /* MODEL_H */
