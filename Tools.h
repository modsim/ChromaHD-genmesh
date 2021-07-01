#ifndef TOOLS_H
#define TOOLS_H

#include <vector>

void subtractTags(std::vector<int>& mainVec, std::vector<int>& subVec);
void addTags(std::vector<int>& mainVec, const std::vector<int> subVec);

void subtractDimTags(std::vector<std::pair<int,int>>& mainVec, std::vector<std::pair<int,int>> subVec);
void addDimTags(std::vector<std::pair<int,int>>& mainVec, const std::vector<std::pair<int,int>> subVec);
void printDimTags(std::vector<std::pair<int,int>> dimTags);

void extractSurfacesWithNormal(std::vector<std::pair<int,int>> dt_input, std::vector<double> _normals, std::vector<std::pair<int,int>>& dt_output);

void findCutBeads(std::vector<std::pair<int,int>> dt_beads, std::vector<std::pair<int,int>> dt_output);
void findUncutBeads(std::vector<std::pair<int,int>> dt_beads, std::vector<std::pair<int,int>> dt_output);

void filterIfSurfaceWithNormal(std::vector<std::pair<int,int>> dt_input, std::vector<double> _normals, std::vector<std::pair<int,int>>& dt_output);

void printDimTagsMap(const std::vector<std::vector<std::pair<int,int>>>& ovv);

void extractTags(const std::vector<std::pair<int,int>>& dt, std::vector<int>& tags);

#endif
