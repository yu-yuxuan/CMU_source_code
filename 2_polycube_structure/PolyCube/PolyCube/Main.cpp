#include <stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "cxxopts.hpp"
#include "atlstr.h" 
#include "PolyCube.h"

using namespace std;

int main(int argc, char* argv[])
{




	string fn_in;
	string fn_out;
	string fn_corner_out;
	string fn_in_normal;
	string fn_in_patch;
	string fn_in_patch_k;
	string fn_in_normal_k;
	string fn_manual_file;
	string fn_input_polycube;
	string fn_input_polycube_k;
	bool fn_input_polycube_flag = false;
	string other_normal_information;
	string normal_projection_edit;
	bool flag_quality = false;
	bool fn_in_patch_flag = false;
	bool normaledit_flag = false;
	string fn_out_k;
	string fn_out_vtk;
	int iteration_polycube = 100;
	int flag_sharp = 0;
	int flag_fit = 0;
	int ngref = 0;
	int nlref = 0;
	int manual_seg = 0;
	double step_polycube = 0.2;
	vector<double> other_normal_vector;
	bool flag_spline = false;
	bool flag_analysis = false;




	try
	{
		cxxopts::Options options(argv[0], "CMU Polycube construction ");
		options
			.positional_help("[optional args]")
			.show_positional_help();

		options.add_options("General")
			("h,help", "Print help")
			("i,input", "File name of geometry without type of mesh and extention", cxxopts::value<std::string>(fn_in))
			("I,iteration", "Iteration number, default is 100", cxxopts::value<int>(iteration_polycube))
			("s,step", "Step, default is 0.2", cxxopts::value<double>(step_polycube))
			("n,input_seg", "Normal Segmentation file name without extention", cxxopts::value<std::string>(fn_in_normal))
			("p,index_patch", "Refine segmentation patch file name without extention", cxxopts::value<std::string>(fn_in_patch))
			//("f,manual-file", "The manual file name", cxxopts::value<std::string>(fn_manual_file))
			("o,output", "The output file with vtk format (opened by paraview)", cxxopts::value<std::string>(fn_out))
			("e,normaledit", "Change normal projection way 0: default 1: center", cxxopts::value<std::string>(normal_projection_edit))
			("f,input_polycube", "The edited polycube", cxxopts::value<std::string>(fn_input_polycube))
			("c,corner_output", "The output file of cornerpoint", cxxopts::value<std::string>(fn_corner_out))


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
		if (!result.count("step"))
		{
			//fn_out_k = fn_in + "_initial_read.k";
			//fn_out_vtk = fn_in + "_initial_read.vtk";
			cout << "Need step information" << endl;
			exit(1);
		}
		if (!result.count("iteration"))
		{
			//fn_out_k = fn_in + "_initial_read.k";
			//fn_out_vtk = fn_in + "_initial_read.vtk";
			cout << "Need iteration information" << endl;
			exit(1);
		}
		if (!result.count("corner_output"))
		{
			//fn_out_k = fn_in + "_initial_read.k";
			//fn_out_vtk = fn_in + "_initial_read.vtk";
			cout << "Need corner point output information" << endl;
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
			fn_out_k = fn_out + ".k";
			fn_out_vtk = fn_out + ".vtk";
		}
		if (result.count("normaledit"))
		{
			//fn_out_k = fn_in + "_initial_read.k";
			//fn_out_vtk = fn_in + "_initial_read.vtk";
			normaledit_flag = true;
		}
		
		if (!result.count("input_seg"))
		{
			cout << "Need normal information" << endl;
			exit(1);
		}
		else
		{
			fn_in_normal_k = fn_in_normal + ".k";
		}

		if (result.count("index_patch"))
		{
			fn_in_patch_flag = true;
			fn_in_patch_k = fn_in_patch + ".k";
		}
		if (result.count("input_polycube"))
		{
			fn_input_polycube_flag = true;
			fn_input_polycube_k = fn_input_polycube + ".k";
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
	//string inputmeshName = "sphere";
	//string inputmeshName = "bust_tilt_spcube";
	//cout << inputmeshName << endl;
	string outputmeshName;
	vector<double> cornerPoints;
	PolyCube AxisAligned;

	AxisAligned.Polycube_initial_simple(inputmeshName);


	//AxisAligned.ReadManualFile(fn_input_polycube_k, fn_corner_out);




	



	/*for (int i = 0; i < 2; i++)
	{
		cout << normal_special[i][0] << normal_special[i][1] << normal_special[i][2] << endl;
	}
*/
	fn_in_patch_flag = true;
	fn_in_patch_k = fn_in_patch + ".k";
	if (normaledit_flag)
	{
		AxisAligned.Polycube_deformation(inputmeshName, fn_in_normal_k, fn_in_patch_flag, fn_in_patch_k, iteration_polycube, step_polycube, normal_projection_edit);

	}
	else
	{
		AxisAligned.Polycube_deformation(inputmeshName, fn_in_normal_k, fn_in_patch_flag, fn_in_patch_k, iteration_polycube, step_polycube);

	}
	
	//for (int loopi = 0; loopi <cornerPoints.size(); loopi++)
	//{
	//	cout << cornerPoints[loopi] << "," ;
	//}
	
	AxisAligned.GetIntGeometry(fn_corner_out);





	

	if (fn_input_polycube_flag)
	{
		AxisAligned.ReadManualFile(fn_input_polycube_k, fn_corner_out);
		AxisAligned.Write(fn_out_vtk.c_str());
		AxisAligned.WriteKFileIndexPatch(fn_out_k.c_str());
	}
	else
	{
		AxisAligned.Write(fn_out_vtk.c_str());
		AxisAligned.WriteKFileIndexPatch(fn_out_k.c_str());
	}


#pragma region Generators for segmentation
	/*PolyCube GetGenerators;

	GetGenerators.CreateVertexID(inputmeshName);*/
#pragma endregion

	
	return 0;

}