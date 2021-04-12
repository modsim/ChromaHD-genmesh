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
template<typename T> void writeVecToBin(std::vector<T> vec, std::ofstream& outfile);
void writeIntVecToBin(std::vector<int> vec, std::ofstream& outfile);
void readBinToIntVec(std::vector<int>& vec, std::ifstream& infile);
std::string remove_extension(const std::string& filename);
std::string get_extension(const std::string& filename);
void create_directory(const std::string& path);

#endif /* FILES_H */
