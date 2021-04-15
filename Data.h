#ifndef DATA_H
#define DATA_H

#include <vector>

struct Volumes {
    std::vector<int> beads, interstitial, all;
};

struct Surfaces {
    std::vector<int> beads, inlet, outlet, walls;
} ;

struct Walls {
    std::vector<int> xleft, yleft, zleft, xright, yright, zright;
};

#endif
