#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include "histogram.h"

namespace aux
{
    /*String-Manipulation
    *********************************************************/
    std::string zfill_int2string(int inint, const unsigned int &zfill)
    {
        std::string outstring = std::to_string(inint);
        while(outstring.length() < zfill)
            outstring = "0" + outstring;
        return outstring;
    }

    /*Numpy-like
    *********************************************************/
    std::vector<double> linspace(double startval, double endval, uint64_t bins)
    {
        std::vector<double> linspaced(bins);
        double delta = (endval-startval)/(bins-1);
        for(uint64_t i = 0; i < (bins-1); i++)
        {
            linspaced[i] = startval + delta * i;
        }
        linspaced[bins-1] = endval;
        return linspaced;
    }

    double get_otsuthreshold(std::vector<float> &greyscales, float minval, float maxval)
    {
        int n_bins = 1000;

        //get a detailed histogram
        histo::Histogram histo;
        std::vector<double> binedges, bincenters;
        std::vector<uint64_t> bincounts;

        histo.calculatehistogram(greyscales, n_bins, minval, maxval, bincounts, binedges);
        bincenters = histo.binedges2bincenter(binedges);

        uint64_t n_voxels = std::accumulate(bincounts.begin(), bincounts.end(), 0);
        double inputsum = 0.;
        for (uint64_t idx = 0; idx < bincounts.size(); idx++)
            inputsum += bincounts[idx]*bincenters[idx];

        //get otsu threshold
        double sumB = 0., wB = 0., wF = 0., mB = 0., mF = 0., best = 0., bw = 0., th1 = 0., th2 = 0.;
        for (uint16_t idx = 0; idx < bincounts.size(); idx++)
        {
            wB += bincounts[idx];
            double center = bincenters[idx];

            if (wB == 0) continue;

            wF = n_voxels-wB;
            if (wF == 0) break;

            sumB += center * bincounts[idx];
            mB = sumB/wB;
            mF = (inputsum-sumB)/wF;
            bw = wB*wF*(mB-mF)*(mB-mF);

            if (bw >= best)
            {
                th1 = center;
                if (bw > best)
                    th2 = center;
                best = bw;
            }
        }

        return (th1+th2)/2.;
    }
}

