#include "Tools.h"

/* #include <iterator> */
#include <algorithm>

void subtractTags(std::vector<int>& mainVec, std::vector<int>& subVec)
{
    /*
     * Remove bead surfaces that are also part of the wall
     * from being categorized as bead surfaces.
     * These 'tSBeads' will go on to be doubled in `gmsh2mixdv2 -d 4`
     */
    for (auto it : subVec)
    {
        auto it2 = std::find(mainVec.begin(), mainVec.end(), it);
        if (it2 != mainVec.end()) mainVec.erase(it2);
    }

}
