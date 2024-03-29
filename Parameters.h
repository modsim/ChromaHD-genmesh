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

    int meshSizeMethod                             = 1;
    int outputFragments                            = 1;
    int MeshSmoothing                              = 1;
    int GeneralNumThreads                          = 0;
    int GeneralVerbosity                           = 5;
    int MeshCharacteristicLengthExtendFromBoundary = 1;
    int MeshCharacteristicLengthFromCurvature      = 1;
    int MeshCharacteristicLengthFromPoints         = 1;
    int GeometryOCCParallel                        = 1;
    int MeshAlgorithm                              = 5;
    int MeshAlgorithm3D                            = 10;
    int MeshRefineSteps                            = 10;
    int MeshOptimize                               = 1;
    int MeshHighOrderOptimize                      = 0;
    int MeshElementOrder                           = 1;
    int MeshOptimizeNetgen                         = 1;
    int MeshGenerate                               = 3;
    int booleanOperation                           = 0;
    int fragment                                   = 1;
    int NamedBeadVolume                            = 1;
    int NamedInterstitialVolume                    = 1;
    int NamedBeadSurface                           = 1;
    int NamedOuterSurface                          = 1;
    int nBeads                                     = 0;
    int dryRun                                     = 0;
    int MeshMaxNumThreads                          = 0;
    int MeshMinimumCirclePoints                    = 7;
    int packingPrecision                           = 4;
    int fixPorosityMethod                          = 2;
    int containerShape                             = 0;
    int autoContainment                            = 1;

    double zBot                           = 0.0;
    double zTop                           = 0.0;
    double rCyl                           = 0.0;
    double xCyl                           = 0.0;
    double yCyl                           = 0.0;
    double inlet                          = -1.0;
    double outlet                         = -1.0;
    double rFactor                        = 1.0;
    double lc                             = 0.0;        // unused
    double lc_beads                       = 0.0;
    double lc_bridge                      = 0.0;
    double lc_out                         = 0.0;
    double refBeadRadius                  = 0.0;
    double bridgeTol                      = -999;
    double relativeBridgeRadius           = 0.0;
    double bridgeOffsetRatio              = 0.0;
    double GeometryScalingFactor          = 1.0;
    double MeshScalingFactor              = 1.0;
    double preScalingFactor               = 1.0;
    double MeshSmoothRatio                = 1.8;            // TODO: Why?
    double GeometryTolerance              = 1e-8;
    double GeometryToleranceBoolean       = 0.0;            // TODO: Why?
    double fieldThresholdMinFactor        = 1.00;
    double fieldThresholdMaxFactor        = 1.00;
    double MeshCharacteristicLengthFactor = 1;
    double MeshOptimizeThreshold          = 0.3;            // TODO: Why?
    double MeshCharacteristicLengthMin    = 0;
    double MeshCharacteristicLengthMax    = 1e22;
    double bridged                        = -1;
    double reduced                        = -1;
    double capped                         = -1;
    double rCylDelta                      = 0.0;
    double por_target                     = 0.0;
    double por_eps                        = DBL_MAX;
    double x0                             = 0.0;
    double y0                             = 0.0;
    double z0                             = 0.0;
    double dx                             = 0.0;
    double dy                             = 0.0;
    double dz                             = 0.0;
    double tOffX                          = 0.0;    // translateOffsets
    double tOffY                          = 0.0;
    double tOffZ                          = 0.0;
    double pOffX                          = 0.0;    // periodicOffsets
    double pOffY                          = 0.0;
    double pOffZ                          = 0.0;
    double periodicInlet = 0.0;
    double periodicOutlet = 0.0;


    std::string periodic         = "off";
    std::string translateOffsets = "auto"; // to maintain some backwards compatibility
    std::string periodicOffsets  = "auto";
    std::string outpath          =".";
    std::string refBeadSize      = "avg";
    std::string fragmentFormat   = "vtk";

    std::string packfile;
    std::string geomOutfile;
    std::string geomInfile;


private:
void decide(const std::string & key, const std::vector<std::string> & val);

};

#endif	/* PARAMETERS_H */

