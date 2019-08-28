/*
 *
 *
 *
 */

#ifndef FILES_H
#define FILES_H

#include <vector>
#include <fstream>
#include <string.h>

bool isBigEndian();
void swapbytes(char *, int, int);
void writeIntVecToBin(std::vector<int> vec, std::ofstream& outfile);
void readBinToIntVec(std::vector<int>& vec, std::ifstream& infile);
std::string remove_extension(const std::string& filename);
void create_directory(const std::string& path);

#endif /* FILES_H */
