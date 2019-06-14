/*
 * File:   Parameters.cpp
 * Author: guowei
 *
 * Created on September 24, 2012, 10:34 AM
 */

#include "Parameters.h"
#include "mixd.hpp"
#include <sstream>

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

    file.close();
}


Parameters::Parameters(const Parameters& orig) {
}

Parameters::~Parameters() {
}


void Parameters::decide (const std::string & key, const std::vector<std::string> & val)
{

    if(key == "zBot")                                        zBot                                       = atof(val.at(0).c_str());

    else if(key == "zTop")                                        zTop                                       = atof(val.at(0).c_str());
    else if(key == "rCyl")                                        rCyl                                       = atof(val.at(0).c_str());
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
    else if(key == "nbeads")                                      nBeadsInPack                               = atoi(val.at(0).c_str());
    else if(key == "bridgeOffsetRatio")                           bridgeOffsetRatio                          = atof(val.at(0).c_str());
    else if(key == "booleanOperation")                            booleanOperation                           = atoi(val.at(0).c_str());
    else if(key == "fragment")                                    fragment                                   = atoi(val.at(0).c_str());
    else if(key == "Named.beadVolume")                            NamedBeadVolume                            = atoi(val.at(0).c_str());
    else if(key == "Named.interstitialVolume")                    NamedInterstitialVolume                    = atoi(val.at(0).c_str());
    else if(key == "Named.beadSurface")                           NamedBeadSurface                           = atoi(val.at(0).c_str());
    else if(key == "Named.inlet")                                 NamedInlet                                 = atoi(val.at(0).c_str());
    else if(key == "Named.outlet")                                NamedOutlet                                = atoi(val.at(0).c_str());
    else if(key == "Named.wall")                                  NamedWall                                  = atoi(val.at(0).c_str());
    else if(key == "Named.outerSurface")                          NamedOuterSurface                          = atoi(val.at(0).c_str());
    else if(key == "General.NumThreads")                          GeneralNumThreads                          = atoi(val.at(0).c_str());
    else if(key == "Geometry.OCCParallel")                        GeometryOCCParallel                        = atoi(val.at(0).c_str());
    else if(key == "Geometry.ScalingFactor")                      GeometryScalingFactor                      = atof(val.at(0).c_str());
    else if(key == "Mesh.ScalingFactor")                          MeshScalingFactor                          = atof(val.at(0).c_str());
    else if(key == "Mesh.Smoothing")                              MeshSmoothing                              = atof(val.at(0).c_str());
    else if(key == "Mesh.SmoothRatio")                            MeshSmoothRatio                            = atof(val.at(0).c_str());
    else if(key == "Mesh.CharacteristicLengthExtendFromBoundary") MeshCharacteristicLengthExtendFromBoundary = atoi(val.at(0).c_str());
    else if(key == "Mesh.CharacteristicLengthMin")                MeshCharacteristicLengthMin                = atof(val.at(0).c_str());
    else if(key == "Mesh.CharacteristicLengthFromCurvature")      MeshCharacteristicLengthFromCurvature      = atoi(val.at(0).c_str());
    else if(key == "Mesh.CharacteristicLengthFromPoints")         MeshCharacteristicLengthFromPoints         = atoi(val.at(0).c_str());
    else if(key == "Mesh.Algorithm")                              MeshAlgorithm                              = atoi(val.at(0).c_str());
    else if(key == "Mesh.Algorithm3D")                            MeshAlgorithm3D                            = atoi(val.at(0).c_str());
    else if(key == "Mesh.Optimize")                               MeshOptimize                               = atoi(val.at(0).c_str());
    else if(key == "Mesh.OptimizeNetgen")                         MeshOptimizeNetgen                         = atoi(val.at(0).c_str());
    else if(key == "Mesh.RefineSteps")                            MeshRefineSteps                            = atoi(val.at(0).c_str());
    else if(key == "Mesh.Generate")                               MeshGenerate                               = atoi(val.at(0).c_str());
    else if(key == "dryRun")                                      dryRun                                     = atoi(val.at(0).c_str());

    else if(key == "packing") packfile = val.at(0);
    else if(key == "outpath") outpath  = val.at(0);

    /* else if(key == "copybeads") */
    /* { */
    /*     if(val.at(0) == "on") copyBeads = true; */
    /*     else if(val.at(0) == "off") copyBeads = false; */
    /*     else throw mixd::MixdException("Unknown option for copybeads: " + val.at(0)); */
    /* } */
    /* else if(key == "periodic") */
    /* { */
    /*     if(val.at(0) == "on") periodic = true; */
    /*     else if(val.at(0) == "off") periodic = false; */
    /*     else throw mixd::MixdException("Unknown option for periodic: " + val.at(0)); */
    /* } */
    /* else if(key == "dryRun") */
    /* { */
    /*     if(val.at(0) == "on") dryRun = true; */
    /*     else if(val.at(0) == "off") dryRun = false; */
    /*     else throw mixd::MixdException("Unknown option for dryRun: " + val.at(0)); */
    /* } */
    else
        throw mixd::MixdException("Unknown keyword: " + key);
}

void Parameters::print()
{

    std::cout << "zBot                                        "<< this->zBot                                       << std::endl;
    std::cout << "zTop                                        "<< this->zTop                                       << std::endl;
    std::cout << "rCyl                                        "<< this->rCyl                                       << std::endl;
    std::cout << "xCyl                                        "<< this->xCyl                                       << std::endl;
    std::cout << "yCyl                                        "<< this->yCyl                                       << std::endl;
    std::cout << "inlet                                       "<< this->inlet                                      << std::endl;
    std::cout << "outlet                                      "<< this->outlet                                     << std::endl;
    std::cout << "rFactor                                     "<< this->rFactor                                    << std::endl;
    std::cout << "preScalingFactor                            "<< this->preScalingFactor                           << std::endl;
    std::cout << "bridgeTol                                   "<< this->bridgeTol                                  << std::endl;
    std::cout << "relativeBridgeRadius                        "<< this->relativeBridgeRadius                       << std::endl;
    std::cout << "lc                                          "<< this->lc                                         << std::endl;
    std::cout << "lc_beads                                    "<< this->lc_beads                                   << std::endl;
    std::cout << "nbeads                                      "<< this->nBeadsInPack                               << std::endl;
    std::cout << "bridgeOffsetRatio                           "<< this->bridgeOffsetRatio                          << std::endl;
    std::cout << "booleanOperation                            "<< this->booleanOperation                           << std::endl;
    std::cout << "fragment                                    "<< this->fragment                                   << std::endl;
    std::cout << "Named.beadVolume                            "<< this->NamedBeadVolume                            << std::endl;
    std::cout << "Named.interstitialVolume                    "<< this->NamedInterstitialVolume                    << std::endl;
    std::cout << "Named.beadSurface                           "<< this->NamedBeadSurface                           << std::endl;
    std::cout << "Named.outerSurface                          "<< this->NamedOuterSurface                          << std::endl;
    /* std::cout << "Named.inlet                                 "<< this->NamedInlet                                 << std::endl; */
    /* std::cout << "Named.outlet                                "<< this->NamedOutlet                                << std::endl; */
    /* std::cout << "Named.wall                                  "<< this->NamedWall                                  << std::endl; */
    std::cout << "General.NumThreads                          "<< this->GeneralNumThreads                          << std::endl;
    std::cout << "Geometry.OCCParallel                        "<< this->GeometryOCCParallel                        << std::endl;
    std::cout << "Geometry.ScalingFactor                      "<< this->GeometryScalingFactor                      << std::endl;
    std::cout << "Mesh.ScalingFactor                          "<< this->MeshScalingFactor                          << std::endl;
    std::cout << "Mesh.CharacteristicLengthExtendFromBoundary "<< this->MeshCharacteristicLengthExtendFromBoundary << std::endl;
    std::cout << "Mesh.CharacteristicLengthMin                "<< this->MeshCharacteristicLengthMin                << std::endl;
    std::cout << "Mesh.CharacteristicLengthFromCurvature      "<< this->MeshCharacteristicLengthFromCurvature      << std::endl;
    std::cout << "Mesh.CharacteristicLengthFromPoints         "<< this->MeshCharacteristicLengthFromPoints         << std::endl;
    std::cout << "Mesh.Smoothing                              "<< this->MeshSmoothing                              << std::endl;
    std::cout << "Mesh.SmoothRatio                            "<< this->MeshSmoothRatio                            << std::endl;
    std::cout << "Mesh.Algorithm                              "<< this->MeshAlgorithm                              << std::endl;
    std::cout << "Mesh.Algorithm3D                            "<< this->MeshAlgorithm3D                            << std::endl;
    std::cout << "Mesh.Optimize                               "<< this->MeshOptimize                               << std::endl;
    std::cout << "Mesh.OptimizeNetgen                         "<< this->MeshOptimizeNetgen                         << std::endl;
    std::cout << "Mesh.RefineSteps                            "<< this->MeshRefineSteps                            << std::endl;
    std::cout << "Mesh.Generate                               "<< this->MeshGenerate                               << std::endl;
    std::cout << "packing                                     "<< this->packfile                                   << std::endl;
    std::cout << "outpath                                     "<< this->outpath                                    << std::endl;
    std::cout << "dryRun                                      "<< this->dryRun                                     << std::endl;

}

