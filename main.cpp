#include <iostream>
#include <omp.h>
#include <math.h>

#include <libgen.h>
#include <unistd.h>
#include <linux/limits.h>

#include "hdcommunication.h"
#include "auxiliary.h"
#include "histogram.h"

using namespace std;

int main(int argc, char* argv[])
{
    string inpath = "";
    string outpath = "";

    string threshhold_mode = "otsu"; //only mode right now. kmeans would be reasonable addition
    int stepsize = 50;
    int radius = 200;

    bool full_reconstruction = true; //assuming 0.0 outside the reconstructed area
    float rel_inside = 0.9; //if closer to center 0.0 is still a valid value
    uint8_t oob_value = 0; //cvalue outside reconstructed area

    bool save_thresholds = false;

    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    string rootpath;

    if ("extract comman line arguments)")
    {
        for (uint16_t i = 1; i < argc; i++)
        {
            if ((string(argv[i]) == "-i") || (string(argv[i]) == "-input"))
            {
                i++;
                inpath = string(argv[i]);
            }
            else if((string(argv[i]) == "-o") || (string(argv[i]) == "-output"))
            {
                i++;
                outpath = string(argv[i]);
            }
            else if((string(argv[i]) == "-step") || (string(argv[i]) == "-stepsize"))
            {
                i++;
                stepsize = atoi(argv[i]);
            }
            else if((string(argv[i]) == "-r") || (string(argv[i]) == "-radius"))
            {
                i++;
                radius = atoi(argv[i]);
            }
            else if (string(argv[i]) == "--subregion")
            {
                full_reconstruction = false;
            }
            else if (string(argv[i]) == "--demo")
            {
                char result[PATH_MAX];
                ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
                const char *tmppath;
                if (count != -1) {
                    tmppath = dirname(result);
                }
                inpath = string(tmppath) +"/Demo/";
            }
        }
        if (outpath.length() == 0){
            rootpath = inpath.substr(0, inpath.rfind("/", inpath.length()-2)+1);
            outpath = rootpath + "/localseg/";}
        else
            rootpath = outpath.substr(0, outpath.rfind("/", outpath.length()-2)+1);

        outpath += "/";
    }

    //Grab data
    ////////////////////////////////////////////////////////////
    hdcom::HdCommunication hdcom;

    int shape[3];
    float *input = hdcom.GetTif_unknowndim_32bit(inpath, shape, true);

    long long int nslice = shape[0]*shape[1];
    long long int nstack = shape[2]*nslice;

    float center[3] = {shape[0]/2., shape[1]/2., shape[2]/2.};

    int n_steps[3] = {max(2,shape[0]/stepsize+1), max(2,shape[1]/stepsize+1), max(2,shape[2]/stepsize+1)};
    long long int nslice_th = n_steps[0]*n_steps[1];
    long long int nstack_th = n_steps[2]*nslice_th;

    float *th_map = (float*) malloc(nstack_th*sizeof(*th_map));
    ////////////////////////////////////////////////////////////

    //Global histobounds
    ////////////////////////////////////////////////////////////
    cout << "global histobounds: ";
    cout.flush();

    std::vector<float> temp;
    for (long long int idx = 0; idx < nstack; idx++)
    {
        int z = idx/nslice;
        int y = (idx-nslice*z)/shape[0];
        int x = (idx-nslice*z-y*shape[0]);

        if (full_reconstruction)
        {
            float sqdist = (x-center[0])*(x-center[0])+(y-center[1])*(y-center[1]);

            if (sqdist > rel_inside/2.f*shape[0] && input[idx] == 0.0f)
                continue;
        }

        temp.push_back(input[idx]);
    }
    histo::Histogram histo;
    std::pair<float, float> histobounds = histo.get_effectivehistobounds(temp, 1000);
    cout << histobounds.first << " " << histobounds.second << endl;

    temp.clear();
    ////////////////////////////////////////////////////////////

    //create sparse thresholds
    ////////////////////////////////////////////////////////////
    cout << "creating sparse thresholds...";
    cout.flush();

    //omp_set_num_threads(4);

    #pragma omp parallel for
    for (long long int thidx = 0; thidx < nstack_th; thidx++)
    {
        int th_z = thidx/nslice_th;
        int th_y = (thidx-th_z*nslice_th)/n_steps[0];
        int th_x = (thidx-th_z*nslice_th-th_y*n_steps[0]);

        int z0 = th_z*(shape[2]-1)/((float) (n_steps[2]-1));
        int y0 = th_y*(shape[1]-1)/((float) (n_steps[1]-1));
        int x0 = th_x*(shape[0]-1)/((float) (n_steps[0]-1));

        //Grab valid values
        std::vector<float> values;
        uint64_t invalid = 0;
        uint64_t valid = 0;

        for (int z1 = max(0, z0-radius); z1 < min(z0+radius, shape[2]); z1++)
        {
            for (int y1 = max(0, y0-radius); y1 < min(y0+radius, shape[1]); y1++)
            {
                for (int x1 = max(0, x0-radius); x1 < min(x0+radius, shape[0]); x1++)
                {
                    //not spherical
                    float sqdist = (x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0);
                    if (sqdist > radius*radius) continue;

                    float val = input[z1*nslice+y1*shape[0]+x1];

                    //out of bounds
                    if (full_reconstruction)
                    {
                        sqdist = (x1-center[0])*(x1-center[0])+(y1-center[1])*(y1-center[1]);
                        if (sqdist > rel_inside/2.f*shape[0] && val == 0.0f)
                        {
                            invalid++;
                            continue;
                        }
                    }

                    valid++;
                    values.push_back(val);
                }
            }
        }

        if (values.size() > 1 && valid > invalid)
        {
            double this_th = aux::get_otsuthreshold(values, histobounds.first, histobounds.second);
            th_map[thidx] = this_th;
        }
        else
            th_map[thidx] = histobounds.second+1e-6f;
    }
    cout << endl;
    ////////////////////////////////////////////////////////////

    //hdcom.SaveTifSequence_32bit(th_map, n_steps, outpath, "thresholds", true);
    //return 0;

    //Linear interpolation
    ////////////////////////////////////////////////////////////
    cout << "interpolating..." << endl;

    std::vector<uint8_t> output(nstack, 0.0);

    #pragma omp parallel for
    for (long long int idx = 0; idx < nstack; idx++)
    {
        int z = idx/nslice;
        int y = (idx-nslice*z)/shape[0];
        int x = (idx-nslice*z-y*shape[0]);

        if (full_reconstruction)
        {
            float sqdist = (x-center[0])*(x-center[0])+(y-center[1])*(y-center[1]);

            if (sqdist > rel_inside/2.f*shape[0] && input[idx] == 0.0f)
            {
                input[idx] = histobounds.second;
                output[idx] = oob_value;
                continue;
            }
        }

        float th_x = ((float) x)/(shape[0]-1)*(n_steps[0]-1);
        float th_y = ((float) y)/(shape[1]-1)*(n_steps[1]-1);
        float th_z = ((float) z)/(shape[2]-1)*(n_steps[2]-1);

        int xf = std::floor(th_x);
        int yf = std::floor(th_y);
        int zf = std::floor(th_z);

        int xc = min(n_steps[0]-1, xf+1);
        int yc = min(n_steps[1]-1, yf+1);
        int zc = min(n_steps[2]-1, zf+1);

        float wx = th_x-xf;
        float wy = th_y-yf;
        float wz = th_z-zf;

        float val0 = th_map[zf*nslice_th + yf*n_steps[0] + xf];
        float val1 = th_map[zf*nslice_th + yf*n_steps[0] + xc];
        float val2 = th_map[zf*nslice_th + yc*n_steps[0] + xf];
        float val3 = th_map[zf*nslice_th + yc*n_steps[0] + xc];

        float val4 = th_map[zc*nslice_th + yf*n_steps[0] + xf];
        float val5 = th_map[zc*nslice_th + yf*n_steps[0] + xc];
        float val6 = th_map[zc*nslice_th + yc*n_steps[0] + xf];
        float val7 = th_map[zc*nslice_th + yc*n_steps[0] + xc];

        if (val0 > histobounds.second) val0 = val1;
        if (val1 > histobounds.second) val1 = val0;
        if (val2 > histobounds.second) val2 = val3;
        if (val3 > histobounds.second) val3 = val2;

        if (val4 > histobounds.second) val4 = val5;
        if (val5 > histobounds.second) val5 = val4;
        if (val6 > histobounds.second) val6 = val7;
        if (val7 > histobounds.second) val7 = val6;

        val0 = (1.f-wx)*val0 + wx*val1;
        val1 = (1.f-wx)*val2 + wx*val3;
        val2 = (1.f-wx)*val4 + wx*val5;
        val3 = (1.f-wx)*val6 + wx*val7;

        if (val0 > histobounds.second) val0 = val1;
        if (val1 > histobounds.second) val1 = val0;
        if (val2 > histobounds.second) val2 = val3;
        if (val3 > histobounds.second) val3 = val2;

        val0 = (1.f-wy)*val0 + wy*val1;
        val1 = (1.f-wy)*val2 + wy*val3;

        if (val0 > histobounds.second) val0 = val1;
        if (val1 > histobounds.second) val1 = val0;

        val0 = (1.f-wz)*val0 + wz*val1;

        float inval = input[idx];

        if (save_thresholds)
            input[idx] = val0;

        if (val0 > histobounds.second) output[idx] = oob_value;
        else if (inval <= val0) output[idx] = 0;
        else output[idx] = 128;
    }
    ////////////////////////////////////////////////////////////

    //Save
    if (save_thresholds)
        hdcom.SaveTifSequence_32bit(input, shape, outpath +"/localthreshholds/", "thresholds", true);

    hdcom.SaveTifSequence_8bit(output, shape, outpath, "segmented", 0);

    if ("append logfile"){
        time_t now = time(0);
        ofstream logfile;
        logfile.open(rootpath + "/logfile.txt", fstream::in | fstream::out | fstream::app);
        logfile << ctime(&now);
        logfile << "ran LocalThresholding:\n";
        logfile << "-------------------------------------------------------------------------\n";
        logfile << "    - threshold_mode: " << threshhold_mode << "\n";
        logfile << "    - stepsize: " << stepsize << "\n";
        logfile << "    - radius: " << radius << "\n";
        if (full_reconstruction)
            logfile << "    - oob_excluded: true \n";
        else
            logfile << "    - oob_excluded: false \n";
        logfile << "-------------------------------------------------------------------------\n\n";
        logfile.close();
    }

    return 0;
}
