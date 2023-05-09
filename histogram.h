#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <vector>
#include <string.h>
#include <iostream>
#include <cstdint>

namespace histo
{
    typedef float dtype;

    class Histogram
    {
    public:
        double histocutoff = 0.001;

        std::pair<dtype,dtype> get_effectivehistobounds(std::vector<dtype> &data, int n_bins);
        void calculateeffectivehistogram(std::vector<dtype> &data, int n_bins, std::vector<uint64_t> &out_histobins, std::vector<double> &out_histoedges);
        void calculateeffectivehistogram(std::vector<dtype> &data, int n_bins, std::vector<double> &out_histobins, std::vector<double> &out_histoedges);
        void calculatehistogram(std::vector<dtype> &data, int n_bins, dtype lowerbound, dtype upperbound, std::vector<uint64_t> &out_histobins, std::vector<double> &out_histoedges);
        void calculatehistogram(std::vector<dtype> &data, int n_bins, dtype lowerbound, dtype upperbound, std::vector<double> &out_histobins, std::vector<double> &out_histoedges);
        std::vector<double> binedges2bincenter(std::vector<double> &binedges);
        std::vector<double> normalize(std::vector<uint64_t> histobins, std::string normalization);
        std::vector<double> normalize(std::vector<double> histobins, std::string normalization);

        void set_radius(int radius_, int shape_[3]){
            radius = radius_;
            shape[0] = shape_[0];
            shape[1] = shape_[1];
            shape[2] = shape_[2];
        }

    private:
        int shape[3] = {0,0,0}; //set the shape if you set a radius to exclude non-reconstructed
        int radius = 0; //exclude out of this radius (also need to set shape)

        std::vector<double> linspace(double startval, double endval, uint64_t bins);
    };
}

#endif // HISTOGRAM_H

