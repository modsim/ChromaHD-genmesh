/*
 * File:   Parameters.cpp
 * Author: guowei
 *
 * Created on September 24, 2012, 10:34 AM
 */

#include "Parameters.h"
#include "mixd.hpp"
#include <sstream>

Parameters::Parameters(const std::string& fname) :
  zBot(0.0), zTop(0.0), rCyl(0.0), xCyl(0.0), yCyl(0.0),
  zCylMin(1.0), zCylMax(-1.0), inlet(-1.0), outlet(-1.0), radius(0.0),
  rFactor(0.0), lc(0.0), lc_beads(0.0), nBeadsInPack(0),
  copyBeads(true), periodic(false), packfile(""), outpath(""), bridgeTol(0.0), db_dp(0.0), dryRun(true), bridgeOffsetFactor(0.0),
  GeometryOCCParallel(1), MeshAlgorithm(2), MeshAlgorithm3D(1), MeshRefineSteps(10), MeshOptimizeNetgen(1), MeshGenerate(3)
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


void Parameters::decide
(const std::string & key, const std::vector<std::string> & val)
throw (mixd::MixdException)
{

         if(key == "zBot")                             zBot                 = atof(val.at(0).c_str());
    else if(key == "zTop")                             zTop                 = atof(val.at(0).c_str());
    else if(key == "rCyl")                             rCyl                 = atof(val.at(0).c_str());
    else if(key == "xCyl")                             xCyl                 = atof(val.at(0).c_str());
    else if(key == "yCyl")                             yCyl                 = atof(val.at(0).c_str());
    else if(key == "inlet")                            inlet                = atof(val.at(0).c_str());
    else if(key == "outlet")                           outlet               = atof(val.at(0).c_str());
    else if(key == "zCylMin")                          zCylMin              = atof(val.at(0).c_str());
    else if(key == "zCylMax")                          zCylMax              = atof(val.at(0).c_str());
    else if(key == "rFactor")                          rFactor              = atof(val.at(0).c_str());
    else if(key == "bridgeTol")                        bridgeTol            = atof(val.at(0).c_str());
    else if(key == "db_dp")                            db_dp                = atof(val.at(0).c_str());
    else if(key == "lc")                               lc                   = atof(val.at(0).c_str());
    else if(key == "lc_beads")                         lc_beads             = atof(val.at(0).c_str());
    else if(key == "nbeads")                           nBeadsInPack         = atoi(val.at(0).c_str());
    else if(key == "bridgeOffsetFactor")               bridgeOffsetFactor   = atof(val.at(0).c_str());

    else if(key == "fuseBeadsAndBridges")               fuseBeadsAndBridges            =  atoi(val.at(0).c_str());
    else if(key == "cutBeadsAndBridges")                cutBeadsAndBridges             =  atoi(val.at(0).c_str());
    else if(key == "fragment")                          fragment                       =  atoi(val.at(0).c_str());

    else if(key == "Named.beadVolume")           NamedBeadVolume                               =  atoi(val.at(0).c_str());
    else if(key == "Named.interstitialVolume")   NamedInterstitialVolume                       =  atoi(val.at(0).c_str());
    else if(key == "Named.beadSurface")          NamedBeadSurface                              =  atoi(val.at(0).c_str());
    else if(key == "Named.inlet")                NamedInlet                                    =  atoi(val.at(0).c_str());
    else if(key == "Named.outlet")               NamedOutlet                                   =  atoi(val.at(0).c_str());
    else if(key == "Named.wall")                 NamedWall                                     =  atoi(val.at(0).c_str());
    else if(key == "Named.outerSurface")         NamedOuterSurface                                     =  atoi(val.at(0).c_str());


    else if(key == "Geometry.OCCParallel")             GeometryOCCParallel  = atof(val.at(0).c_str());
    else if(key == "Mesh.Algorithm")                   MeshAlgorithm        = atof(val.at(0).c_str());
    else if(key == "Mesh.Algorithm3D")                 MeshAlgorithm3D      = atof(val.at(0).c_str());
    else if(key == "Mesh.OptimizeNetgen")              MeshOptimizeNetgen   = atof(val.at(0).c_str());
    else if(key == "Mesh.RefineSteps")                 MeshRefineSteps      = atof(val.at(0).c_str());
    else if(key == "Mesh.Generate")                    MeshGenerate         = atof(val.at(0).c_str());

    else if(key == "packing")                          packfile             = val.at(0);
    else if(key == "outpath")                          outpath              = val.at(0);
    else if(key == "copybeads")

    {
        if(val.at(0) == "on") copyBeads = true;
        else if(val.at(0) == "off") copyBeads = false;
        else throw mixd::MixdException("Unknown option for copybeads: " + val.at(0));
    }
    else if(key == "periodic")
    {
        if(val.at(0) == "on") periodic = true;
        else if(val.at(0) == "off") periodic = false;
        else throw mixd::MixdException("Unknown option for periodic: " + val.at(0));
    }
    else if(key == "dryRun")
    {
        if(val.at(0) == "on") dryRun = true;
        else if(val.at(0) == "off") dryRun = false;
        else throw mixd::MixdException("Unknown option for dryRun: " + val.at(0));
    }
    else
        throw mixd::MixdException("Unknown keyword: " + key);
}

void Parameters::print()
{

    std::cout << "zBot:                    " << this->zBot << std::endl;
    std::cout << "zTop:                    " << this->zTop << std::endl;
    std::cout << "rCyl:                    " << this->rCyl << std::endl;
    std::cout << "xCyl:                    " << this->xCyl << std::endl;
    std::cout << "yCyl:                    " << this->yCyl << std::endl;
    std::cout << "zCylMin:                 " << this->zCylMin << std::endl;
    std::cout << "zCylMax:                 " << this->zCylMax << std::endl;
    std::cout << "inlet:                   " << this->inlet << std::endl;
    std::cout << "outlet:                  " << this->outlet << std::endl;
    std::cout << "radius:                  " << this->radius << std::endl;
    std::cout << "rFactor:                 " << this->rFactor << std::endl;
    std::cout << "bridgeTol:               " << this->bridgeTol << std::endl;
    std::cout << "db_dp:                   " << this->db_dp << std::endl;
    std::cout << "bridgeOffsetFactor:      " << this->bridgeOffsetFactor << std::endl;
    std::cout << "lc:                      " << this->lc << std::endl;
    std::cout << "lc_beads:                " << this->lc_beads << std::endl;
    std::cout << "nBeads:                  " << this->nBeadsInPack << std::endl;
    std::cout << "copyBeads:               " << this->copyBeads << std::endl;
    std::cout << "periodic:                " << this->periodic << std::endl;
    std::cout << "packfile:                " << this->packfile << std::endl;
    std::cout << "outpath:                 " << this->outpath << std::endl;
    std::cout << "dryRun:                  " << this->dryRun << std::endl;

}

