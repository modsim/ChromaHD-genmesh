/*
 * File: Parameters.h
 * Desc: Parameter class handles inputs from mesher.in
 * Created by: guowei
 *
 * Created on September 24, 2012, 10:34 AM
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include "mixd.hpp"
#include <vector>

class Parameters {
public:
    Parameters(const std::string& fname);
    Parameters(const Parameters& orig);
    virtual ~Parameters();

    void print();
    void write(std::string filename);

    double zBot, zTop, rCyl, xCyl, yCyl;
    double zCylMin, zCylMax;
    double inlet, outlet;
    double radius, rFactor;
    double lc, lc_beads;
    double bridgeTol, db_dp;

    double bridgeOffsetFactor;

    int GeometryOCCParallel;
    int MeshAlgorithm, MeshAlgorithm3D;
    int MeshRefineSteps;
    int MeshOptimize;
    int MeshOptimizeNetgen;
    int MeshGenerate;

    int fuseBeadsAndBridges;
    int cutBeadsAndBridges;
    int fragment;

    int NamedBeadVolume;
    int NamedInterstitialVolume;
    int NamedBeadSurface;
    int NamedInlet;
    int NamedOutlet;
    int NamedWall;
    int NamedOuterSurface;


    int nBeadsInPack;

    bool dryRun;

    bool copyBeads, periodic;

    std::string packfile, outpath;


private:
void decide(const std::string & key, const std::vector<std::string> & val)
throw (mixd::MixdException);

};

#endif	/* PARAMETERS_H */

