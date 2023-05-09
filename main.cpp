#include <iostream>
#include <unistd.h>
#include <chrono>
#include <omp.h>

#include "Geometry/hdcommunication.h"
#include "Geometry/auxiliary.h"

#include "anisodiff_cpufilters.h"
#include "anisodiff_gpufilters.h"

using namespace std;

/*********************************************************************************************************************************************************
 *
 * Location: Helmholtz-Zentrum fuer Material und Kuestenforschung, Max-Planck-Strasse 1, 21502 Geesthacht
 * Author: Stefan Bruns
 * Contact: bruns@nano.ku.dk
 *
 * License: TBA
 *
 *********************************************************************************************************************************************************/

int main(int argc, char* argv[])
{
	/////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////

	string inpath = "";
	string outpath = "";

	bool cpu_mode = false;
	int gpu0 = 0;
	int slice_overlap = 50; //when processing chunks create an overlap region
	int n_threads = 128; //number of CPU_threads (-1 to use all)
	bool output_as_16bit = false;
	bool twopoint5D = false;

	int n_iters = 50;
	int n_smoothings = 2;
	bool threeDfilter = false;

	float a1 = 0.5;
	float a2 = 0.9;
	float a3 = 1.3;
	float dt = 20;
	float edge_th = 5.0;

	int firstslice = 0;
	int lastslice = -1;

	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////

	if ("extract command line arguments"){
		for (uint16_t i = 1; i < argc; i++)
		{
			if ((string(argv[i]) == "-i") || (string(argv[i]) == "-input"))
			{
				i++;
				inpath = string(argv[i]);
			}
			else if ((string(argv[i]) == "-o") || (string(argv[i]) == "-outpath"))
			{
				i++;
				outpath = string(argv[i]);
			}
			else if (string(argv[i]) == "-n_cpu" || string(argv[i]) == "-n_threads")
			{
				i++;
				n_threads = atoi(argv[i]);
			}
			else if (string(argv[i]) == "--16bit")
				output_as_16bit = true;
			else if (string(argv[i]) == "--cpu")
				cpu_mode = true;
			else if (string(argv[i]) == "--3D")
				threeDfilter = true;
			else if (string(argv[i]) == "--2.5D")
				twopoint5D = true;
			else if (string(argv[i]) == "-gpu0")
			{
				i++;
				gpu0 = atoi(argv[i]);
			}
			else if (string(argv[i]) == "-iter" || string(argv[i]) == "-n_iter" || string(argv[i]) == "-n_iters"  || string(argv[i]) == "-iters")
			{
				i++;
				n_iters = atoi(argv[i]);
			}
			else if (string(argv[i]) == "-smoothings" || string(argv[i]) == "-n_smoothings" || string(argv[i]) == "-smoothing" || string(argv[i]) == "-n_smoothing")
			{
				i++;
				n_smoothings = atoi(argv[i]);
			}
			else if (string(argv[i]) == "-a1")
			{
				i++;
				a1 = atof(argv[i]);
			}
			else if (string(argv[i]) == "-a2")
			{
				i++;
				a2 = atof(argv[i]);
			}
			else if (string(argv[i]) == "-a3")
			{
				i++;
				a3 = atof(argv[i]);
			}
			else if (string(argv[i]) == "-dt")
			{
				i++;
				dt = atof(argv[i]);
			}
			else if (string(argv[i]) == "-edge_th")
			{
				i++;
				edge_th = atof(argv[i]);
			}
			else if (string(argv[i]) == "-overlap")
			{
				i++;
				slice_overlap = atof(argv[i]);
			}
			else if (string(argv[i]) == "-zrange")
			{
				i++;
				firstslice = atoi(argv[i]);
				i++;
				lastslice = atoi(argv[i]);
			}
		}
	}

	string rootpath = inpath.substr(0, inpath.rfind("/", inpath.length()-2)+1);
	if (outpath.length() == 0){
		rootpath = inpath.substr(0, inpath.rfind("/", inpath.length()-2)+1);
		outpath = rootpath + "/denoised/";}
	else
		rootpath = outpath.substr(0, outpath.rfind("/", outpath.length()-2)+1);

	if (n_threads > 0) omp_set_num_threads(min(n_threads, omp_get_max_threads()));

	////////////////////////////////////////////////////////////

	cout << "------------------------------" << endl;
	cout << "in: " << inpath << endl;
	cout << "out: " << outpath << endl;

	auto time0 = chrono::high_resolution_clock::now();
	hdcom::HdCommunication hdcom;
	int shape[3];
	float* imgstack;

	if (firstslice <= 0 && lastslice < 0) imgstack = hdcom.GetTif_unknowndim_32bit(inpath, shape, true);
	else
	{
		std::vector<string> filelist = hdcom.GetFilelist(inpath, shape);
		imgstack = hdcom.Get3DTifSequence_32bitPointer(filelist,shape,firstslice,lastslice);
	}

	if (cpu_mode || threeDfilter)
	{
		filter::AnisotropicDiffusion anisofilter;

		anisofilter.n_iters = n_iters;
		anisofilter.n_smoothings = n_smoothings;
		anisofilter.a1 = a1;
		anisofilter.a2 = a2;
		anisofilter.a3 = a3;
		anisofilter.dt = dt;
		anisofilter.edge_threshold = edge_th;

		if (threeDfilter)
		{
			std::cout << "Using 3D-Filter" << std::endl;
			anisofilter.TschumperleDeriche_3D(imgstack, shape);
		}
		else if (twopoint5D)
		{
			std::cout << "Filtering with alternating directions" << std::endl;
			anisofilter.TschumperleDeriche_2D_AlternatingDirections(imgstack, shape);
		}
		else
		{
			std::cout << "Filtering in xy" << std::endl;
			anisofilter.TschumperleDeriche_2D(imgstack, shape);
		}
	}
	else
	{
		gpufilter::AnisotropicDiffusion anisofilter;

		anisofilter.n_iters = n_iters;
		anisofilter.n_smoothings = n_smoothings;
		anisofilter.a1 = a1;
		anisofilter.a2 = a2;
		anisofilter.dt = dt;
		anisofilter.edge_threshold = edge_th;
		anisofilter.overlap = slice_overlap;

		anisofilter.configure_device(shape, gpu0);

		if (twopoint5D)
		{
			std::cout << "Filtering with alternating directions" << std::endl;
			anisofilter.run_TschumperleDeriche_2D(imgstack, shape, true);
		}
		else
		{
			std::cout << "Filtering in xy" << std::endl;
			anisofilter.run_TschumperleDeriche_2D(imgstack, shape);
		}


		anisofilter.free_device();
	}

	if(!output_as_16bit)
		hdcom.SaveTifSequence_32bit(std::max(0, firstslice), imgstack, shape, outpath, "anisofiltered", true);
	else
		hdcom.SaveTifSequence_as16bit(std::max(0, firstslice), imgstack,shape,outpath,"anisofiltered", true);

	auto time1 = chrono::high_resolution_clock::now();

	chrono::duration<double> elapsed_total = time1-time0;
	std::cout << "execution took " << elapsed_total.count() << " s" << std::endl;
	cout << "------------------------------" << endl;
	return 0;
}
