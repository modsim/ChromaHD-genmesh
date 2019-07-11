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
    void update();

    double zBot=0.0, zTop=0.0, rCyl=0.0, xCyl=0.0, yCyl=0.0;
    /* double zCylMin=1.0, zCylMax=-1.0; */
    double inlet=-1.0, outlet=-1.0;
    double radius=0.0, rFactor=1.0;
    double lc=0.0, lc_beads=0.0;
    double lc_max = 0.0, lc_min = 0.0;
    double bridgeTol=-999, relativeBridgeRadius=0.0;
    double bridgeOffsetRatio=0.0;

    double GeometryScalingFactor=1.0;
    double MeshScalingFactor=1.0;
    double preScalingFactor=1.0;

    int MeshSmoothing=1;
    double MeshSmoothRatio=1.8;

    double GeometryTolerance = 1e-8;
    double GeometryToleranceBoolean = 0.0;

    int GeometryOCCParallel=1;
    int MeshAlgorithm=2, MeshAlgorithm3D=1;
    int MeshRefineSteps=10;
    int MeshOptimize=1;
    int MeshOptimizeNetgen=0;
    int MeshGenerate=3;

    double fieldExtensionFactor = 1.00;


    int GeneralNumThreads=8;

    int MeshCharacteristicLengthExtendFromBoundary=1;
    int MeshCharacteristicLengthFromCurvature=1;
    int MeshCharacteristicLengthFromPoints=1;

    int MeshMaxNumThreads = 0;
    double MeshCharacteristicLengthFactor = 1;
    int MeshMinimumCirclePoints = 7;
    double MeshOptimizeThreshold = 0.3;

    double MeshCharacteristicLengthMin=0;
    double  MeshCharacteristicLengthMax=1e22;

    double bridged = -1;
    double reduced = -1;
    double capped = -1;

    int beadType = 0; //

    int booleanOperation    = 0;
    /* int fuseBeadsAndBridges = 0; */
    /* int cutBeadsAndBridges  = 0; */
    int fragment            = 1;

    int NamedBeadVolume=1;
    int NamedInterstitialVolume=1;
    int NamedBeadSurface=1;
    int NamedOuterSurface=1;

    int NamedInlet;
    int NamedOutlet;
    int NamedWall;

    int nBeads=0;
    int dryRun=0;

    bool copyBeads, periodic;

    std::string packfile, outpath="output/";


private:
void decide(const std::string & key, const std::vector<std::string> & val);

};

#endif	/* PARAMETERS_H */

