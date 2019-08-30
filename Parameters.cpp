/*
 * File:   Parameters.cpp
 * Author: Jayghosh Rao
 *
 * Based on the file created by guowei on September 24, 2012, 10:34 AM
 */

#include "Parameters.h"
#include "mixd.hpp"
#include <sstream>
#include <cmath>
#include<gmsh.h>

Parameters::Parameters(const std::string& fname)
{
    using namespace std;
    using namespace mixd;

    ifstream file(fname.c_str(), ifstream::in);

    if(!file.is_open())
        throw MixdException("could not open file " + fname);

    string line;

    while(getline(file, line))
    {
        line = ltrim(line);

        // skip empty and comment lines
        if(valid(line))
        {
            pair< string,vector<string> > keyval = getKeyVal(line);

            string         key = keyval.first;
            vector<string> val = keyval.second;

            decide(key, val);
        }
    }

    /*
     * If bridgeOffsetRatio is not defined, calculate it ourselves.
     * If booleanOperation == 2 (Cut), then round UP to two decimal places.
     *      Else round DOWN to two decimal places.
     */
    if (bridgeOffsetRatio == 0)
    {
        bridgeOffsetRatio = ( ( sqrt(1 - pow(relativeBridgeRadius,2)) * 100 ) + (booleanOperation == 2)*1 - (booleanOperation != 2)*2 ) / 100;
        std::cout << "#bridgeOffsetRatio automatically calculated!" << std::endl;

    }

    print();
    update();
    setGMSHOptions();




    file.close();
}


Parameters::Parameters(const Parameters& orig) {
}

Parameters::~Parameters() {
}

void Parameters::update()
{
    //using a separate method to update instead of putting the code
    //in the constructor allows us to print the parameters in the stdout
    //and re-use the parameters directly in the next simulation.

    //remnants of a previous way to handle packings with different scales
    /* bridgeTol = preScalingFactor * bridgeTol; */
    /* lc        = preScalingFactor * lc; */
    /* lc_beads  = preScalingFactor * lc_beads; */
    /* lc_out    = preScalingFactor * lc_out; */


    // zBot and zTop in the inputs are based on the column with bead size = 1
    zBot      = zBot / preScalingFactor;
    zTop      = zTop / preScalingFactor;


}

void Parameters::decide (const std::string & key, const std::vector<std::string> & val)
{

    if(key == "zBot")                                             zBot                                       = atof(val.at(0).c_str());

    else if(key == "zTop")                                        zTop                                       = atof(val.at(0).c_str());
    else if(key == "rCyl")                                        rCyl                                       = atof(val.at(0).c_str());
    else if(key == "rCylDelta")                                   rCylDelta                                  = atof(val.at(0).c_str());
    else if(key == "xCyl")                                        xCyl                                       = atof(val.at(0).c_str());
    else if(key == "yCyl")                                        yCyl                                       = atof(val.at(0).c_str());
    else if(key == "inlet")                                       inlet                                      = atof(val.at(0).c_str());
    else if(key == "outlet")                                      outlet                                     = atof(val.at(0).c_str());
    else if(key == "preScalingFactor")                            preScalingFactor                           = atof(val.at(0).c_str());
    else if(key == "rFactor")                                     rFactor                                    = atof(val.at(0).c_str());
    else if(key == "bridgeTol")                                   bridgeTol                                  = atof(val.at(0).c_str());
    else if(key == "relativeBridgeRadius")                        relativeBridgeRadius                       = atof(val.at(0).c_str());
    else if(key == "lc")                                          lc                                         = atof(val.at(0).c_str());
    else if(key == "lc_beads")                                    lc_beads                                   = atof(val.at(0).c_str());
    else if(key == "lc_out")                                      lc_out                                     = atof(val.at(0).c_str());
    else if(key == "fieldExtensionFactor")                        fieldExtensionFactor                       = atof(val.at(0).c_str());
    else if(key == "nBeads")                                      nBeads                                     = atoi(val.at(0).c_str());
    else if(key == "meshSizeMethod")                              meshSizeMethod                             = atoi(val.at(0).c_str());
    else if(key == "outputFragments")                             outputFragments                            = atoi(val.at(0).c_str());
    else if(key == "bridgeOffsetRatio")                           bridgeOffsetRatio                          = atof(val.at(0).c_str());
    else if(key == "booleanOperation")                            booleanOperation                           = atoi(val.at(0).c_str());
    else if(key == "fragment")                                    fragment                                   = atoi(val.at(0).c_str());
    else if(key == "Named.beadVolume")                            NamedBeadVolume                            = atoi(val.at(0).c_str());
    else if(key == "Named.interstitialVolume")                    NamedInterstitialVolume                    = atoi(val.at(0).c_str());
    else if(key == "Named.beadSurface")                           NamedBeadSurface                           = atoi(val.at(0).c_str());
    else if(key == "Named.outerSurface")                          NamedOuterSurface                          = atoi(val.at(0).c_str());
    else if(key == "General.NumThreads")                          GeneralNumThreads                          = atoi(val.at(0).c_str());
    else if(key == "Geometry.OCCParallel")                        GeometryOCCParallel                        = atoi(val.at(0).c_str());
    else if(key == "Geometry.ScalingFactor")                      GeometryScalingFactor                      = atof(val.at(0).c_str());
    else if(key == "Geometry.Tolerance")                          GeometryTolerance                          = atof(val.at(0).c_str());
    else if(key == "Geometry.ToleranceBoolean")                   GeometryToleranceBoolean                   = atof(val.at(0).c_str());
    else if(key == "Mesh.ScalingFactor")                          MeshScalingFactor                          = atof(val.at(0).c_str());
    else if(key == "Mesh.Smoothing")                              MeshSmoothing                              = atof(val.at(0).c_str());
    else if(key == "Mesh.SmoothRatio")                            MeshSmoothRatio                            = atof(val.at(0).c_str());
    else if(key == "Mesh.CharacteristicLengthExtendFromBoundary") MeshCharacteristicLengthExtendFromBoundary = atoi(val.at(0).c_str());
    else if(key == "Mesh.CharacteristicLengthMin")                MeshCharacteristicLengthMin                = atof(val.at(0).c_str());
    else if(key == "Mesh.CharacteristicLengthMax")                MeshCharacteristicLengthMax                = atof(val.at(0).c_str());
    else if(key == "Mesh.CharacteristicLengthFromCurvature")      MeshCharacteristicLengthFromCurvature      = atoi(val.at(0).c_str());
    else if(key == "Mesh.CharacteristicLengthFromPoints")         MeshCharacteristicLengthFromPoints         = atoi(val.at(0).c_str());
    else if(key == "Mesh.CharacteristicLengthFactor")             MeshCharacteristicLengthFactor             = atof(val.at(0).c_str());
    else if(key == "Mesh.MinimumCirclePoints")                    MeshMinimumCirclePoints                    = atoi(val.at(0).c_str());
    else if(key == "Mesh.MaxNumThreads")                          MeshMaxNumThreads                          = atof(val.at(0).c_str());
    else if(key == "Mesh.OptimizeThreshold")                      MeshOptimizeThreshold                      = atof(val.at(0).c_str());
    else if(key == "Mesh.Algorithm")                              MeshAlgorithm                              = atoi(val.at(0).c_str());
    else if(key == "Mesh.Algorithm3D")                            MeshAlgorithm3D                            = atoi(val.at(0).c_str());
    else if(key == "Mesh.Optimize")                               MeshOptimize                               = atoi(val.at(0).c_str());
    else if(key == "Mesh.OptimizeNetgen")                         MeshOptimizeNetgen                         = atoi(val.at(0).c_str());
    else if(key == "Mesh.RefineSteps")                            MeshRefineSteps                            = atoi(val.at(0).c_str());
    else if(key == "Mesh.Generate")                               MeshGenerate                               = atoi(val.at(0).c_str());
    else if(key == "dryRun")                                      dryRun                                     = atoi(val.at(0).c_str());

    /* else if(key == "Named.inlet")                                 NamedInlet                                 = atoi(val.at(0).c_str()); */
    /* else if(key == "Named.outlet")                                NamedOutlet                                = atoi(val.at(0).c_str()); */
    /* else if(key == "Named.wall")                                  NamedWall                                  = atoi(val.at(0).c_str()); */


    else if(key == "reduced")
    {
        rFactor = atof(val.at(0).c_str());
        bridgeTol = -999;
        booleanOperation = 0;
    }
    else if(key == "bridged")
    {
        relativeBridgeRadius = atof(val.at(0).c_str());
        bridgeTol = lc/10;
        booleanOperation = 1;
    }
    else if(key == "capped")
    {
        relativeBridgeRadius = atof(val.at(0).c_str());
        bridgeTol = lc/10;
        booleanOperation = 2;
    }

    else if(key == "packing") packfile = val.at(0);
    else if(key == "geomInfile") geomInfile = val.at(0);
    else if(key == "geomOutfile") geomOutfile = val.at(0);
    else if(key == "outpath") outpath  = val.at(0);
    else throw mixd::MixdException("Unknown keyword: " + key);



}

void Parameters::print()
{

    std::cout << std::setprecision(10);
    std::cout << "zBot                                        "<< this->zBot                                       << std::endl;
    std::cout << "zTop                                        "<< this->zTop                                       << std::endl;
    std::cout << "rCyl                                        "<< this->rCyl                                       << std::endl;
    std::cout << "rCylDelta                                   "<< this->rCylDelta                                  << std::endl;
    std::cout << "xCyl                                        "<< this->xCyl                                       << std::endl;
    std::cout << "yCyl                                        "<< this->yCyl                                       << std::endl;
    std::cout << "inlet                                       "<< this->inlet                                      << std::endl;
    std::cout << "outlet                                      "<< this->outlet                                     << std::endl;
    std::cout << "nBeads                                      "<< this->nBeads                                     << std::endl;
    std::cout << "meshSizeMethod                              "<< this->meshSizeMethod                             << std::endl;
    std::cout << "outputFragments                             "<< this->outputFragments                            << std::endl;
    std::cout << "rFactor                                     "<< this->rFactor                                    << std::endl;
    std::cout << "preScalingFactor                            "<< this->preScalingFactor                           << std::endl;
    std::cout << "lc                                          "<< this->lc                                         << std::endl;
    std::cout << "lc_beads                                    "<< this->lc_beads                                   << std::endl;
    std::cout << "lc_out                                      "<< this->lc_out                                     << std::endl;
    std::cout << "fieldExtensionFactor                        "<< this->fieldExtensionFactor                       << std::endl;
    std::cout << "bridgeTol                                   "<< this->bridgeTol                                  << std::endl;
    std::cout << "relativeBridgeRadius                        "<< this->relativeBridgeRadius                       << std::endl;
    std::cout << "bridgeOffsetRatio                           "<< this->bridgeOffsetRatio                          << std::endl;
    std::cout << "booleanOperation                            "<< this->booleanOperation                           << std::endl;
    std::cout << "fragment                                    "<< this->fragment                                   << std::endl;
    std::cout << "Named.beadVolume                            "<< this->NamedBeadVolume                            << std::endl;
    std::cout << "Named.interstitialVolume                    "<< this->NamedInterstitialVolume                    << std::endl;
    std::cout << "Named.beadSurface                           "<< this->NamedBeadSurface                           << std::endl;
    std::cout << "Named.outerSurface                          "<< this->NamedOuterSurface                          << std::endl;
    std::cout << "General.NumThreads                          "<< this->GeneralNumThreads                          << std::endl;
    std::cout << "Geometry.OCCParallel                        "<< this->GeometryOCCParallel                        << std::endl;
    std::cout << "Geometry.ScalingFactor                      "<< this->GeometryScalingFactor                      << std::endl;
    std::cout << "Geometry.Tolerance                          "<< this->GeometryTolerance                          << std::endl;
    std::cout << "Geometry.ToleranceBoolean                   "<< this->GeometryToleranceBoolean                   << std::endl;
    std::cout << "Mesh.ScalingFactor                          "<< this->MeshScalingFactor                          << std::endl;
    std::cout << "Mesh.CharacteristicLengthExtendFromBoundary "<< this->MeshCharacteristicLengthExtendFromBoundary << std::endl;
    std::cout << "Mesh.CharacteristicLengthMin                "<< this->MeshCharacteristicLengthMin                << std::endl;
    std::cout << "Mesh.CharacteristicLengthMax                "<< this->MeshCharacteristicLengthMax                << std::endl;
    std::cout << "Mesh.CharacteristicLengthFromCurvature      "<< this->MeshCharacteristicLengthFromCurvature      << std::endl;
    std::cout << "Mesh.CharacteristicLengthFromPoints         "<< this->MeshCharacteristicLengthFromPoints         << std::endl;
    std::cout << "Mesh.CharacteristicLengthFactor             "<< this->MeshCharacteristicLengthFactor             << std::endl;
    std::cout << "Mesh.MinimumCirclePoints                    "<< this->MeshMinimumCirclePoints                    << std::endl;
    std::cout << "Mesh.OptimizeThreshold                      "<< this->MeshOptimizeThreshold                      << std::endl;
    std::cout << "Mesh.MaxNumThreads                          "<< this->MeshMaxNumThreads                          << std::endl;
    std::cout << "Mesh.Smoothing                              "<< this->MeshSmoothing                              << std::endl;
    std::cout << "Mesh.SmoothRatio                            "<< this->MeshSmoothRatio                            << std::endl;
    std::cout << "Mesh.Algorithm                              "<< this->MeshAlgorithm                              << std::endl;
    std::cout << "Mesh.Algorithm3D                            "<< this->MeshAlgorithm3D                            << std::endl;
    std::cout << "Mesh.Optimize                               "<< this->MeshOptimize                               << std::endl;
    std::cout << "Mesh.OptimizeNetgen                         "<< this->MeshOptimizeNetgen                         << std::endl;
    std::cout << "Mesh.RefineSteps                            "<< this->MeshRefineSteps                            << std::endl;
    std::cout << "Mesh.Generate                               "<< this->MeshGenerate                               << std::endl;
    std::cout << "packing                                     "<< this->packfile                                   << std::endl;
    std::cout << "geomInfile                                  "<< this->geomInfile                                 << std::endl;
    std::cout << "geomOutfile                                 "<< this->geomOutfile                                << std::endl;
    std::cout << "outpath                                     "<< this->outpath                                    << std::endl;
    std::cout << "dryRun                                      "<< this->dryRun                                     << std::endl;

}


void Parameters::setGMSHOptions()
{
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::option::setNumber("General.NumThreads", this->GeneralNumThreads);
    gmsh::option::setNumber("Geometry.OCCParallel", this->GeometryOCCParallel);
    gmsh::option::setNumber("Geometry.ScalingFactor", this->GeometryScalingFactor);
    gmsh::option::setNumber("Geometry.Tolerance", this->GeometryTolerance);
    gmsh::option::setNumber("Geometry.ToleranceBoolean", this->GeometryToleranceBoolean);
    gmsh::option::setNumber("Mesh.ScalingFactor", this->MeshScalingFactor);
    gmsh::option::setNumber("Mesh.Smoothing", this->MeshSmoothing);
    gmsh::option::setNumber("Mesh.SmoothRatio", this->MeshSmoothRatio);

    gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary", this->MeshCharacteristicLengthExtendFromBoundary);

    gmsh::option::setNumber("Mesh.CharacteristicLengthMin", this->MeshCharacteristicLengthMin);
    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", this->MeshCharacteristicLengthMax);

    gmsh::option::setNumber("Mesh.CharacteristicLengthFromCurvature", this->MeshCharacteristicLengthFromCurvature);
    gmsh::option::setNumber("Mesh.CharacteristicLengthFromPoints", this->MeshCharacteristicLengthFromPoints);

    gmsh::option::setNumber("Mesh.MinimumCirclePoints", this->MeshMinimumCirclePoints); //Default = 7

    gmsh::option::setNumber("Mesh.Optimize", this->MeshOptimize); //Default = 1
    gmsh::option::setNumber("Mesh.OptimizeNetgen", this->MeshOptimizeNetgen); //Default = 0
    gmsh::option::setNumber("Mesh.RefineSteps", this->MeshRefineSteps); //Default = 10

    gmsh::option::setNumber("Mesh.MaxNumThreads1D", this->MeshMaxNumThreads);
    gmsh::option::setNumber("Mesh.MaxNumThreads2D", this->MeshMaxNumThreads);
    gmsh::option::setNumber("Mesh.MaxNumThreads3D", this->MeshMaxNumThreads);
    gmsh::option::setNumber("Mesh.CharacteristicLengthFactor", this->MeshCharacteristicLengthFactor);
    gmsh::option::setNumber("Mesh.OptimizeThreshold", this->MeshOptimizeThreshold);

    //1: MeshAdapt | 2: Auto | 5: Delaunay | 6: Frontal | 7: BAMG | 8: DelQuad
    gmsh::option::setNumber("Mesh.Algorithm", this->MeshAlgorithm); //Default = 2

    //1: Delaunay | 4: Frontal | 5: Frontal Delaunay | 6: Frontal Hex | 7: MMG3D | 9: RTree | 10: HXT
    gmsh::option::setNumber("Mesh.Algorithm3D", this->MeshAlgorithm3D); //Default = 1

    gmsh::option::setNumber("Print.GeoLabels", 1); //Default = 1
    gmsh::option::setNumber("Print.GeoOnlyPhysicals", 1); //Default = 1

}
