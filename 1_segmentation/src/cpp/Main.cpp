#include <stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "cxxopts.hpp"
#include "CVTBasedPolycube.h"
#include "atlstr.h" 
using namespace std;


int main(int argc, char* argv[])
{

	string fn_in;
	string fn_out;
	string fn_manual_file;
	bool flag_quality = false;
	string fn_out_k;
	string fn_out_vtk;
	int flag_quality_option = 0;
	int flag_sharp = 0;
	int flag_fit = 0;
	int ngref = 0;
	int nlref = 0;
	int manual_seg = 0;
	double thick = 0.0;

	bool flag_spline = false;
	bool flag_analysis = false;


	try
	{
		cxxopts::Options options(argv[0], "CMU Create Segmentation");
		options
			.positional_help("[optional args]")
			.show_positional_help();

		options.add_options("General")
			("h,help", "Print help")
			("i,input", "File name without type of mesh and extention", cxxopts::value<std::string>(fn_in))
			//("m,manual", "Need the manual information, default is 0", cxxopts::value<int>(manual_seg))
			//("f,manual-file", "The manual file name", cxxopts::value<std::string>(fn_manual_file))
			("o,output", "The output file", cxxopts::value<std::string>(fn_out))
			/*("o,output","Output path", cxxopts::value<std::string>(fn_out))*/
			//("e,extrusion", "Extrusion thickness", cxxopts::value<double>(thick))
#ifdef CXXOPTS_USE_UNICODE
			("unicode", u8"A help option with non-ascii: ��. Here the size of the"
				" string should be correct")
#endif
			;

		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			cout << options.help({ "General" }) << endl;
			exit(0);
		}
		if (!result.count("output"))
		{
			fn_out_k = fn_in + "_initial_read.k";
			fn_out_vtk = fn_in + "_initial_read.vtk";
		}
		else
		{
			fn_out_k = fn_out + ".k";
			fn_out_vtk = fn_out + ".vtk";
		}
		if (result.count("input"))
		{
			char *buffer = strdup(fn_in.c_str());
			char *fn_wo_extension_char = buffer;
			PathRemoveExtensionA(fn_wo_extension_char);
			stringstream ss_in;
			ss_in << fn_wo_extension_char;
			ss_in >> fn_in;

			//stringstream ss_ref;
			//ss_ref << fn_in << "_gref_lev" << ngref << "_lgref_lev" << nlref;
			//ss_ref >> fn_out;

			//if (nlref > 0)
			//{
			//	char *buffer1 = strdup(fn_in.c_str());
			//	char *fn_wo_filespec_char = buffer1;
			//	PathRemoveFileSpecA(fn_wo_filespec_char);
			//	stringstream ss_in;
			//	string fld;
			//	ss_in << fn_wo_filespec_char;
			//	ss_in >> fld;
			//	fn_lref = fld + "\\lev";

			//	/*ss_lref << fn_in << "_lref_lev" << nlref;
			//	ss_lref >> fn_out;*/
			//}

		}

		//if (result.count("output"))
		//{
		//	char *buffer = strdup(fn_in.c_str());
		//	char *fn_wo_extension_char = buffer;
		//	PathRemoveExtensionA(fn_wo_extension_char);
		//	stringstream ss_in;
		//	ss_in << fn_wo_extension_char;
		//	ss_in >> fn_in;

		//	//stringstream ss_ref;
		//	//ss_ref << fn_in << "_gref_lev" << ngref << "_lgref_lev" << nlref;
		//	//ss_ref >> fn_out;

		//	//if (nlref > 0)
		//	//{
		//	//	char *buffer1 = strdup(fn_in.c_str());
		//	//	char *fn_wo_filespec_char = buffer1;
		//	//	PathRemoveFileSpecA(fn_wo_filespec_char);
		//	//	stringstream ss_in;
		//	//	string fld;
		//	//	ss_in << fn_wo_filespec_char;
		//	//	ss_in >> fld;
		//	//	fn_lref = fld + "\\lev";

		//	//	/*ss_lref << fn_in << "_lref_lev" << nlref;
		//	//	ss_lref >> fn_out;*/
		//	//}

		//}


	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}


	string inputmeshName = fn_in;




	string inputFullName;

	string outputmeshName;


	CVTBasedPolycube polycube(inputmeshName);


	double progress = 0.5;

	AddProgress(progress);


	inputFullName = inputmeshName + "_tri.raw";//input file, triangle mesh

	polycube.Initialization(inputFullName.c_str());//generate initial segmentation results

	if (manual_seg > 0)
	{
		progress = 0.7;

		AddProgress(progress);
	}
	else
	{
		progress = 1.0;

		AddProgress(progress);
	}







	//if (manual_seg > 0)
	//{ 
	//	polycube.PostProcessing(fn_manual_file.c_str());//
	//	outputmeshName = fn_out + ".vtk";//improved results with manual adjust
	//	polycube.OutputPatchesVTK(outputmeshName.c_str());
	//	outputmeshName = fn_out + "_segment_initial.k";
	//	polycube.WriteKFileBeforePostProcessing(outputmeshName.c_str());
	//	progress = 1.0;

	//	AddProgress(progress);
	//
	//}
	//else if (manual_seg==0)
	//{
	//outputmeshName = fn_out_vtk;//output initial segmentation results
	//polycube.OutputPatchesVTK(outputmeshName.c_str());

	outputmeshName = fn_out_k;//output initial segmentation results
	polycube.WriteKFileBeforePostProcessing(outputmeshName.c_str());
	/*}*/
		

	return 0;

}