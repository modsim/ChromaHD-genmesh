/*
 *
 *
 *
 */

#ifndef FILES_H
#define FILES_H

#include <vector>
#include <fstream>

bool isBigEndian();
void swapbytes(char *, int, int);
void writeIntVecToBin(std::vector<int> vec, std::ofstream& outfile);
void readBinToIntVec(std::vector<int>& vec, std::ifstream& infile);


#endif /* FILES_H */
