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
#include <float.h>

class Parameters {
public:
    Parameters(const std::string& fname);
    Parameters(const Parameters& orig);
    virtual ~Parameters();

    void print();
    void write(std::string filename);
    void update();
    void setGMSHOptions();

    int meshSizeMethod                              = 1;
    int outputFragments                             = 1;
    int MeshSmoothing                               = 1;
    int GeneralNumThreads                           = 8;
    int MeshCharacteristicLengthExtendFromBoundary  = 1;
    int MeshCharacteristicLengthFromCurvature       = 1;
    int MeshCharacteristicLengthFromPoints          = 1;
    int GeometryOCCParallel                         = 1;
    int MeshAlgorithm                               = 2;
    int MeshAlgorithm3D                             = 1;
    int MeshRefineSteps                             = 10;
    int MeshOptimize                                = 1;
    int MeshHighOrderOptimize                       = 0;
    int MeshElementOrder                            = 1;
    int MeshOptimizeNetgen                          = 0;
    int MeshGenerate                                = 3;
    int beadType                                    = 0;
    int booleanOperation                            = 0;
    int fragment                                    = 1;
    int NamedBeadVolume                             = 1;
    int NamedInterstitialVolume                     = 1;
    int NamedBeadSurface                            = 1;
    int NamedOuterSurface                           = 1;
    int nBeads                                      = 0;
    int dryRun                                      = 0;
    int MeshMaxNumThreads                           = 0;
    int MeshMinimumCirclePoints                     = 7;
    int packingPrecision                            = 4;

    double zBot                                     = 0.0;
    double zTop                                     = 0.0;
    double rCyl                                     = 0.0;
    double xCyl                                     = 0.0;
    double yCyl                                     = 0.0;
    double inlet                                    = -1.0;
    double outlet                                   = -1.0;
    double rFactor                                  = 1.0;
    double lc                                       = 0.0;
    double lc_beads                                 = 0.0;
    double lc_bridge                                = 0.0;
    double lc_out                                   = 0.0;
    double refBeadRadius                            = 0.0;
    double bridgeTol                                = -999;
    double relativeBridgeRadius                     = 0.0;
    double bridgeOffsetRatio                        = 0.0;
    double GeometryScalingFactor                    = 1.0;
    double MeshScalingFactor                        = 1.0;
    double preScalingFactor                         = 1.0;
    double MeshSmoothRatio                          = 1.8;
    double GeometryTolerance                        = 1e-8;
    double GeometryToleranceBoolean                 = 0.0;
    double fieldThresholdMinFactor                  = 1.00;
    double fieldThresholdMaxFactor                  = 1.00;
    double MeshCharacteristicLengthFactor           = 1;
    double MeshOptimizeThreshold                    = 0.3;
    double MeshCharacteristicLengthMin              = 0;
    double MeshCharacteristicLengthMax              = 1e22;
    double bridged                                  = -1;
    double reduced                                  = -1;
    double capped                                   = -1;
    double rCylDelta                                = 0.0;
    double por_target                               = 0.0;
    double por_eps                                  = DBL_MAX;


    std::string packfile;
    std::string geomOutfile;
    std::string geomInfile;
    std::string outpath="output/";
    std::string refBeadSize = "avg";
    std::string fragmentFormat = "vtk";


private:
void decide(const std::string & key, const std::vector<std::string> & val);

};

#endif	/* PARAMETERS_H */

