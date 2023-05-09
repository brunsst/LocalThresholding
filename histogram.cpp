#include "histogram.h"
#include <numeric>
#include <algorithm>
#include <omp.h>

namespace histo
{
    std::pair<dtype,dtype> Histogram::get_effectivehistobounds(std::vector<dtype> &data, int n_bins)
    {
        std::vector<uint64_t> out_histobins;
        std::vector<double> out_histoedges;

        dtype lowerbound = *std::min_element(data.begin(),data.end());
        dtype upperbound = *std::max_element(data.begin(),data.end());

        //First determine the appropriate lower and upper bound from a histogram covering all data
        calculatehistogram(data, n_bins, lowerbound, upperbound, out_histobins, out_histoedges);

        uint64_t running_sum = 0;
        for (int i = 0; i < n_bins; i++)
        {
            running_sum += out_histobins[i];
            if (running_sum/((double) data.size()) >= histocutoff)
            {
                //lower bound found
                lowerbound = out_histoedges[i];
                break;
            }
        }

        running_sum = 0;
        for (int i = n_bins-1; i >= 0; i--)
        {
            running_sum += out_histobins[i];
            if (running_sum/((double) data.size()) >= histocutoff)
            {
                //upper bound found
                upperbound = out_histoedges[i+1];
                break;
            }
        }

        return {lowerbound, upperbound};
    }
    void Histogram::calculateeffectivehistogram(std::vector<dtype> &data, int n_bins, std::vector<uint64_t> &out_histobins, std::vector<double> &out_histoedges)
    {
        dtype lowerbound = *std::min_element(data.begin(),data.end());
        dtype upperbound = *std::max_element(data.begin(),data.end());

        //First determine the appropriate lower and upper bound from a histogram covering all data
        calculatehistogram(data, n_bins, lowerbound, upperbound, out_histobins, out_histoedges);

        uint64_t running_sum = 0;
        for (int i = 0; i < n_bins; i++)
        {
            running_sum += out_histobins[i];
            if (running_sum/((double) data.size()) >= histocutoff)
            {
                //lower bound found
                lowerbound = out_histoedges[i];
                break;
            }
        }

        running_sum = 0;
        for (int i = n_bins-1; i >= 0; i--)
        {
            running_sum += out_histobins[i];
            if (running_sum/((double) data.size()) >= histocutoff)
            {
                //upper bound found
                upperbound = out_histoedges[i+1];
                break;
            }
        }
        //Now recalculate with effctive data range:

        uint64_t splitsize = 2e9;
        if (data.size() > splitsize)
        {
            std::vector<float> block1(data.begin(), data.begin()+splitsize);
            calculatehistogram(block1, n_bins, lowerbound, upperbound, out_histobins, out_histoedges);
            std::vector<float> block2(data.begin()+splitsize, data.end());
            std::vector<uint64_t> tmpbins(out_histobins.size(), 0);
            calculatehistogram(block1, n_bins, lowerbound, upperbound, tmpbins, out_histoedges);
            for (uint64_t idx = 0; idx < out_histobins.size(); idx++)
                out_histobins[idx] += tmpbins[idx];
        }
        else
            calculatehistogram(data, n_bins, lowerbound, upperbound, out_histobins, out_histoedges);

        return;
    }
    void Histogram::calculateeffectivehistogram(std::vector<dtype> &data, int n_bins, std::vector<double> &out_histobins, std::vector<double> &out_histoedges)
    {
        std::vector<uint64_t> out_histobins_tmp;
        calculateeffectivehistogram(data, n_bins, out_histobins_tmp, out_histoedges);
        out_histobins.assign(out_histobins_tmp.size(),0);
        for(uint64_t idx = 0; idx < out_histobins.size(); idx++)
            out_histobins[idx] = (double) out_histobins_tmp[idx];

        return;
    }

    void Histogram::calculatehistogram(std::vector<dtype> &data, int n_bins, dtype lowerbound, dtype upperbound, std::vector<uint64_t> &out_histobins, std::vector<double> &out_histoedges)
    {
        double lb = (double) lowerbound;
        double ub = (double) upperbound;
        double range = ub-lb;

        long long int nslice = shape[0]*shape[1];
        double centerx = shape[0]/2.;
        double centery = shape[1]/2.;

        out_histoedges = linspace(lowerbound,upperbound,n_bins+1);
        out_histobins.clear();
        out_histobins.resize(n_bins,0);

        uint64_t this_bin;
        dtype datavalue;
        //Iterate one time to determine the appropriate lower and upper bound
        for (uint64_t idx = 0; idx < data.size(); idx++)
        {
            datavalue = data[idx];

            if((datavalue < lowerbound) || (datavalue > upperbound))
                continue;

            if (radius > 0)
            {
                int z = idx/nslice;
                int y = (idx-z*nslice)/shape[0];
                int x = (idx-z*nslice-y*shape[0]);

                double sqdist = (x-centerx)*(x-centerx)+(y-centery)*(y-centery);
                if (sqdist >= radius*radius) continue;
            }

            this_bin = (datavalue-lowerbound)/range*n_bins;
            out_histobins[this_bin]++;
        }
        return;
    }
    void Histogram::calculatehistogram(std::vector<dtype> &data, int n_bins, dtype lowerbound, dtype upperbound, std::vector<double> &out_histobins, std::vector<double> &out_histoedges)
    {
        double lb = (double) lowerbound;
        double ub = (double) upperbound;
        double range = ub-lb;

        out_histoedges = linspace(lowerbound,upperbound,n_bins+1);
        out_histobins.clear();
        out_histobins.resize(n_bins,0);

        uint64_t this_bin;
        dtype datavalue;
        //Iterate one time to determine the appropriate lower and upper bound
        for (uint64_t idx = 0; idx < data.size(); idx++)
        {
            datavalue = data[idx];

            if((datavalue < lowerbound) || (datavalue > upperbound))
                continue;

            this_bin = (datavalue-lowerbound)/range*n_bins;
            out_histobins[this_bin]++;
        }
        return;
    }

    std::vector<double> Histogram::binedges2bincenter(std::vector<double> &binedges)
    {
        std::vector<double> bincenters(binedges.size()-1,0.);
        for (int i = 0; i < bincenters.size(); i++)
            bincenters[i] = (binedges[i]+binedges[i+1])/2.;
        return bincenters;
    }

    std::vector<double> Histogram::normalize(std::vector<uint64_t> histobins, std::string normalization)
    {
        std::vector<double> output(histobins.size(), 0.);
        double weight = 1.;

        if (normalization == "area")
            weight = std::accumulate(histobins.begin(), histobins.end(), 0.);
        else if (normalization == "height")
            weight = (double) *std::max_element(histobins.begin(), histobins.end());

        for(uint64_t idx; idx < histobins.size(); idx++)
            output[idx] = ((double) histobins[idx])/weight;

        return output;
    }
    std::vector<double> Histogram::normalize(std::vector<double> histobins, std::string normalization)
    {
        std::vector<double> output(histobins.size(), 0.);
        double weight = 1.;

        if (normalization == "area")
            weight = std::accumulate(histobins.begin(), histobins.end(), 0.);
        else if (normalization == "height")
            weight = *std::max_element(histobins.begin(), histobins.end());

        for(uint64_t idx; idx < histobins.size(); idx++)
            output[idx] = (histobins[idx])/weight;

        return output;
    }

    /*****************************************************************************************/

    std::vector<double> Histogram::linspace(double startval, double endval, uint64_t bins)
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
}

