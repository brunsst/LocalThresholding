#ifndef AUXILIARY_H
#define AUXILIARY_H

#include <iostream>
#include <string.h>
#include <cstdint>
#include <vector>

namespace aux
{
    /*String-Manipulation
    *********************************************************/
    std::string zfill_int2string(int inint, const unsigned int &zfill);

    /*Numpy-like
    *********************************************************/
    std::vector<double> linspace(double startval, double endval, uint64_t bins);

    double get_otsuthreshold(std::vector<float> &greyscales, float minval, float maxval);
}

#endif // AUXILIARY_H

