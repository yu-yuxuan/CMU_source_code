#include <stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "cxxopts.hpp"
#include "GeneralizedCVT.h"
#include "atlstr.h" 
using namespace std;

int main(int argc, char* argv[])
{




	string fn_in;
	string fn_out;
	string fn_corner_in;
	string fn_in_normal;
	string fn_in_patch;
	string fn_in_patch_k;
	string fn_in_normal_k;
	string fn_manual_file;
	string fn_in_polycube;
	bool flag_quality = false;
	string fn_out_k;
	string fn_out_vtk;
	int iteration_polycube = 100;
	int octree_level = 1;
	int flag_patch = 0;
	int ngref = 0;
	int nlref = 0;
	int manual_seg = 0;
	double step_polycube = 0.2;

	bool flag_spline = false;
	bool flag_analysis = false;


	try
	{
		cxxopts::Options options(argv[0], "CMU All-hex mesh generation");
		options
			.positional_help("[optional args]")
			.show_positional_help();

		options.add_options("General")
			("h,help", "Print help")
			("i,input", "File name of geometry without type of mesh and extention", cxxopts::value<std::string>(fn_in))
			("P,polycube_structure", "Polycube structure file name without extention", cxxopts::value<std::string>(fn_in_polycube))
			("n,input_normal", "Normal file name without extention", cxxopts::value<std::string>(fn_in_normal))
			("f,index_patch_flag", "Flag whether use index patch", cxxopts::value<int>(flag_patch))
			("p,index_patch", "Patch file name without extention", cxxopts::value<std::string>(fn_in_patch))
			("c,corner_input", "Cornerpoint file name with extention", cxxopts::value<std::string>(fn_corner_in))
			//("f,manual-file", "The manual file name", cxxopts::value<std::string>(fn_manual_file))
			("o,output", "The output file with vtk format (opened by paraview)", cxxopts::value<std::string>(fn_out))
			/*("o,output","Output path", cxxopts::value<std::string>(fn_out))*/
			("s,octree_subdivision", "Octree subdivision to get hex mesh", cxxopts::value<int>(octree_level))
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
		if (!result.count("polycube_structure"))
		{
			//fn_out_k = fn_in + "_initial_read.k";
			//fn_out_vtk = fn_in + "_initial_read.vtk";
			cout << "Need polycube structure information" << endl;
			exit(1);
		}
		if (!result.count("corner_input"))
		{
			//fn_out_k = fn_in + "_initial_read.k";
			//fn_out_vtk = fn_in + "_initial_read.vtk";
			cout << "Need corner point information" << endl;
			exit(1);
		}
		if (!result.count("output"))
		{
			//fn_out_k = fn_in + "_initial_read.k";
			//fn_out_vtk = fn_in + "_initial_read.vtk";
			cout << "Need output information" << endl;
			exit(1);
		}
		else
		{
			fn_out_k = fn_out;
			fn_out_vtk = fn_out;
		}
		if (!result.count("input_normal"))
		{
			cout << "Need normal information" << endl;
			exit(1);
		}
		else
		{
			fn_in_normal_k = fn_in_normal + ".k";
		}


		if (result.count("index_patch_flag"))
		{
			if (flag_patch!=0)
			{
				flag_patch = 1;
			}
			
			if (!result.count("index_patch"))
			{
				cout << "Need patch index information" << endl;
				exit(1);
			}
			else
			{
				fn_in_patch_k = fn_in_patch + ".k";
			}
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



	
	string inputmeshName = fn_in;//"heli";
	//string test =;//"heli";
	//std::cout << "Test " << result["f"].as<std::string>() << " arguments" << std::endl;
	//std::cout << "Test " << fn_in << " arguments" << std::endl;
	string inputFullName;


	double progress = 0.2;
	
	AddProgress(progress);


	string outputmeshName;

	inputFullName = inputmeshName + "_tri.raw";

	GeneralizedCVT gCVT(inputmeshName);

	gCVT.Initialization();
	//outputmeshName = inputmeshName + "_output_initial.vtk";
	//gCVT.OutputPatchesVTK(outputmeshName.c_str());

	//gCVT.EdgeWeightedCVT();
	
	gCVT.PostProcessing(fn_in_normal_k);
	//outputmeshName = inputmeshName + "_output_CleanUp.vtk";
	//gCVT.OutputPatchesVTK(outputmeshName.c_str());

	progress = 0.4;

	AddProgress(progress);

	string cornerpoints_output = fn_corner_in;
	gCVT.InitializePolycube(fn_in_patch_k, flag_patch, cornerpoints_output);
	//outputmeshName = inputmeshName + "_output_ModifiedBifur.vtk";
	//gCVT.OutputPatchesVTKBifurcation(outputmeshName.c_str()); // Be careful, the results are different with OutputPatchesVTK()!
	//outputmeshName = inputmeshName + "_output_paraMapping.vtk";
	//gCVT.OutputPatchesVTKPara(outputmeshName.c_str());

	

	outputmeshName = inputmeshName + "_paraHexNew_hex.raw";
	gCVT.HexMeshParametricDomain(outputmeshName.c_str());

	

	outputmeshName = inputmeshName + "_realHex_hex.vtk";
	gCVT.HexMeshRealDomain(outputmeshName.c_str());

	progress = 0.6;

	AddProgress(progress);

	outputmeshName = fn_out_vtk;
	gCVT.HexMeshOctreeSubdivision(outputmeshName.c_str(), octree_level);

	progress = 1.0;

	AddProgress(progress);

	//cout << "test" << endl;

	return 0;

}