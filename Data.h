#ifndef DATA_H
#define DATA_H

#include <vector>
#include <string>

struct tagsWalls{
    std::vector<int> tXLeft, tYLeft, tZLeft, tXRight, tYRight, tZRight;
};

struct Fragment{
    double dx, dy, dz;
    tagsWalls Outer;
    tagsWalls Beads;
};

struct Volumes {
    std::vector<int> beads, interstitial, all;
};

struct Surfaces {
    std::vector<int> beads, inlet, outlet, walls;
} ;

struct Walls {
    std::vector<int> xleft, yleft, zleft, xright, yright, zright;
};


// should be part of Geometry?
struct ColumnDataStruct {

    std::string periodic;
    Volumes volumes;
    Surfaces surfaces;
    Walls outerWalls, beadWalls;
};

#endif
