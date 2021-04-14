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
#include "Geometry.h"
#include "Data.h"
#include "Column.h"

#include <vector>
#include <string>

class Model{
    public:
        Model();
        Model(Parameters * prm, Geometry * geom);
        Model(Parameters * prm, std::string geometryFile);
        virtual ~Model();

        Column column, columnInlet, columnOutlet;

        void mesh(std::string outfile, Parameters * prm);
        void write(std::string outfile, Parameters * prm);
};

#endif /* MODEL_H */
