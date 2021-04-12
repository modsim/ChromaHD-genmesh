#include "Files.h"
#include <stdlib.h>

#include <iostream>
#include <sys/stat.h>


bool isBigEndian() {
    short word = 0x4321;
    if ((* (char*) & word) != 0x21)
        return true;
    else
        return false;
}

// swap bytes according to machine's endianness

void swapbytes(char *array, int nelem, int elsize) {
    register int sizet, sizem, i, j;
    char *bytea, *byteb;
    sizet = elsize;
    sizem = sizet - 1;
    bytea = (char*) malloc(sizet);
    byteb = (char*) malloc(sizet);

    for (i = 0; i < nelem; i++) {
        memcpy((void *) bytea, (void *) (array + i * sizet), sizet);
        for (j = 0; j < sizet; j++) byteb[j] = bytea[sizem - j];
        memcpy((void *) (array + i * sizet), (void *) byteb, sizet);
    }

    delete bytea;
    delete byteb;
}

/* void writeIntVecToBin(std::vector<int> vec, std::ofstream& outfile) */
/* { */
/*     size_t size = vec.size(); */
/*     outfile.write(reinterpret_cast<char *>(&size), sizeof(size)); */
/*     for (std::vector<int>::iterator it = vec.begin(); it!=vec.end(); it++) outfile.write( reinterpret_cast<char *>(&*it), sizeof(int)); */
/* } */

template<typename T>
void writeVecToBin(std::vector<T> vec, std::ofstream& outfile)
{
    size_t size = vec.size();
    outfile.write(reinterpret_cast<char *>(&size), sizeof(size));
    for (typename std::vector<T>::iterator it = vec.begin(); it!=vec.end(); it++) outfile.write( reinterpret_cast<char *>(&*it), sizeof(T));
}

void readBinToIntVec(std::vector<int>& vec, std::ifstream& infile)
{

    size_t size;
    infile.read(reinterpret_cast<char *>(&size), sizeof(size));
    int dummy;
    for(int i=0; i<size; i++)
    {
        infile.read( reinterpret_cast<char *>(&dummy), sizeof(int));
        vec.push_back(dummy);
    }

}

std::string remove_extension(const std::string& filename)
{
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
}

std::string get_extension(const std::string& filename)
{
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(lastdot, std::string::npos);
}


void create_directory(const std::string& path)
{
        if (mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
        {
            if( errno == EEXIST )
            {
                // alredy exists
                std::cout << "Output directory exists!" << std::endl;
            }
            else
            {
                // something else
                std::cout << "Error creating output directory! " << strerror(errno) << std::endl;
                /* throw std::runtime_exception( strerror(errno) ); */
            }
        }

}
