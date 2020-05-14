#include "PolyCube.h"
#include "BasicDataStructure.h"
#include "StaticVars.h"
#include <sstream>
#include <iomanip>
#include <set> 


PolyCube::PolyCube(void)
{

}

PolyCube::PolyCube(string inputFileName)
{
}


PolyCube::~PolyCube(void)
{

}

bool PolyCube::ReadManualFile(string fn_input_polycube_k, string fn_corner_out)
{

	int i, j, tempElementSign;
	string oneLine;
	double x_value, y_value, z_value;
	
	fstream input(fn_input_polycube_k.c_str(), fstream::in);
	if (!input)
	{
		cout << "open input file error!!!" << endl;
		return false;
	}

	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);
	 cout << oneLine << endl;
	for (i = 0; i < elementNumber; i++)
	{
		input >> j >> tempElementSign >> j >> j >> j >> j >> j >> j >> j >> j;
		//cout << tempElementSign << endl;
		
	}
	getline(input, oneLine);
	//cout<< oneLine << endl;
	getline(input, oneLine);
	//cout << oneLine << endl;
	getchar();
	for (int i = 0; i < vertexNumber; i++)
	{
		input >> j >> x_value >> y_value >> z_value >> j >> j ;

		
		vertex[i][0] = x_value;
		vertex[i][1] = y_value;
		vertex[i][2] = z_value;
	/*	cout << fn_input_polycube_k << endl;
		cout << x_value << " " << y_value << " " << z_value << " " << endl;
		getchar();*/
	}
	input.close();

	for (int i = 0; i < vertexNumber; i++)
	{

		for (int j = 0; j < 3; j++)
		{
			//vertex[i][j] = Round_off(vertex[i][j], 3);
			vertex[i][j] = round(vertex[i][j] * 10) / 10.0; //todo 

		}

	}

	SetBoundaryVertexSign(0);
	for (int loopi_j = 0; loopi_j < ele_boundary.size(); loopi_j++)
	{
		//cout << ele_boundary[loopi_j][0] << "," << ele_boundary[loopi_j][1] << ",";
		vertexSign[ele_boundary[loopi_j][0]] = 1;
		vertexSign[ele_boundary[loopi_j][1]] = 1;
	}
	for (int loopi_j = 0; loopi_j < cornerPoints.size(); loopi_j++)
	{
		int index = cornerPoints[loopi_j];
		//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
		vertexSign[index] = 2;


	}
	SmoothEdge(10000);
	SmoothSurface(10000);

	ofstream outputCornerPoints;
	outputCornerPoints.open(fn_corner_out.c_str());

	for (int i = 0; i < cornerPoints.size(); i++)
	{
		int index = cornerPoints[i];

		//cout << index << ", " ;

		outputCornerPoints << index << " ";

		for (int j = 0; j < 3; j++)
		{
			outputCornerPoints << vertex[index][j] << " ";
		}

		outputCornerPoints << endl;
	}

	outputCornerPoints.close();


	return true;

}

bool PolyCube::GetIntGeometry(string fn_corner_out)
{
	for (int i = 0; i < vertexNumber ; i++)
	{

		for (int j = 0; j < 3; j++)
		{
			//vertex[i][j] = Round_off(vertex[i][j], 3);
			vertex[i][j] = round(vertex[i][j]*10)/10.0; //todo 

		}

	}

	SetBoundaryVertexSign(0);
	for (int loopi_j = 0; loopi_j < ele_boundary.size(); loopi_j++)
	{
		//cout << ele_boundary[loopi_j][0] << "," << ele_boundary[loopi_j][1] << "," ;
		vertexSign[ele_boundary[loopi_j][0]] = 1;
		vertexSign[ele_boundary[loopi_j][1]] = 1;
	}
	for (int loopi_j = 0; loopi_j < cornerPoints.size(); loopi_j++)
	{
		int index = cornerPoints[loopi_j];
		//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
		vertexSign[index] = 2;


	}
	SmoothEdge(10000);
	SmoothSurface(10000);

	ofstream outputCornerPoints;
	outputCornerPoints.open(fn_corner_out.c_str());

	for (int i = 0; i < cornerPoints.size(); i++)
	{
		int index = cornerPoints[i];

		//cout << index << ", " ;

		outputCornerPoints << index << " ";

		for (int j = 0; j < 3; j++)
		{
			outputCornerPoints << vertex[index][j] << " ";
		}

		outputCornerPoints << endl;
	}

	outputCornerPoints.close();


	return true;
}

double PolyCube::round(double N)
{
	return floor(N + 0.5);
}

double Round_off(double N, double n)
{
	int h;
	double l, a, b, c, d, e, i, j, m, f, g;
	b = N;
	c = floor(N);

	// Counting the no. of digits to the left of decimal point 
	// in the given no. 
	for (i = 0; b >= 1; ++i)
		b = b / 10;

	d = n - i;
	b = N;
	b = b * pow(10, d);
	e = b + 0.5;
	if ((float)e == (float)ceil(b)) {
		f = (ceil(b));
		h = f - 2;
		if (h % 2 != 0) {
			e = e - 1;
		}
	}
	j = floor(e);
	m = pow(10, d);
	j = j / m;
	//cout << "The number after rounding-off is " << j;
	return j;
}

bool PolyCube::Polycube_initial_simple(string inputFileName)
{

	inputName = inputFileName;

	string tempName;

	vector<array<double, 3>> pts;
	vector<array<int, 3>> cnct;

	if (INPUTFILE_RAW)
	{
		tempName = inputName + "_tri.raw";
		Read(tempName.c_str());

		for (int i = 0; i < vertexNumber; i++)
		{

			for (int j = 1; j < 3; j++)
			{

				vertex[i][j] = vertex[i][j];

			}

		}

		
		//cout << vertex[0][0] << endl;
		//cout << element[0][0] << endl;
	}
	else
	{
		//cout << "this step" << endl;
		ReadVtk_tri(inputName, pts, cnct);
		int nNode = pts.size();
		int nElement = cnct.size();

		 CreateNewMesh(TRIANGLE, nNode, nElement);

		for (int i = 0; i < nNode; i++)
		{

			for (int j = 0; j < 3; j++)
			{

				vertex[i][j] = pts[i][j];

			}
		
		}

		for (int i = 0; i < nElement; i++)
		{

			for (int j = 0; j < 3; j++)
			{

				element[i][j] = cnct[i][j];

			}

		}


	}
	
	//tempName = inputName + "_tri.vtk";
	//Write(tempName.c_str());


	//myio.GenerateDir(inputFileName);
	//ss.clear();
	//ss.str("");
	//ss << myio.mesh_write_dir << myio.mesh_name << "_whole_00000_tri.vtk";
	//Write(ss.str().c_str());
	///*tempName = inputName + "_0000_tri.vtk";
	//Write(tempName.str().c_str());*/

	//
	//ss.clear();
	//ss.str("");
	//ss << myio.mesh_write_dir << myio.mesh_name << "_deform_00000_tri.vtk";
	//Write(ss.str().c_str());
	////ReadK_Polycube(tempName.str().c_str()); 

	InitiateEdgeValence();
	InitiateElementValence();
	//SetBoundaryVertexSign(0);

	
	//cout << tempName << elementNumber << vertexNumber << endl;
	//getchar();

	GetRelativeCoordinate();
	//cout << vertexSign[8759]<<endl;
	return true;

}




bool PolyCube::Readtxt(vector<vector<int> >& Nt, string tempName)
{

	string oneLine;
	
	ifstream myFile(tempName);

	int tempElement;

	if (myFile.is_open())
	{

		while (getline(myFile, oneLine))
		{
			std::vector<int> Nt_temp;
			istringstream st(oneLine);
			while (st >> tempElement)
				Nt_temp.push_back(tempElement);

			Nt.push_back(Nt_temp);
		}

		myFile.close();

	}
	else
	{
		cout << "Unable to open file!";
	}

	return true;

}

bool PolyCube::SearchCornerandBoundary(vector<int>& cornerPoints)
{

	InitiateElementValence();

	int i, j;
	int index;
	int vertexCCW;

	//vector<int> cornerPoints; //It is defined as a variable of the class

	for (i = 0; i < vertexNumber; i++)
	{

		if (IsCornerPoint(i))
		{

			cornerPoints.push_back(i);

		}

	}


	//Eigen::VectorXi bnd;
	////Yuxuan_quad::boundary_loop(F, bnd); //longest
	//std::vector<std::vector<int> > bnd_all;
	//igl::boundary_loop(F, bnd_all); //all

	//ele << 1, 1;
	ele_boundary.clear();
	for (int k = 0; k < numberofpatch; k++)
	{

		vector<vector<int> > ele;
		for (i = 0; i < pcube_exyzelement.size(); i++)
		{
			if (pcube_exyzelement[i].indexPatch == k)
			{

				int foo[] = { element[i][0], element[i][1], element[i][2] };
				int n = sizeof(foo) / sizeof(foo[0]);

				sort(foo, foo + n);
				vector<int> vect1 = { foo[0], foo[1] };
				vector<int> vect2 = { foo[1], foo[2] };
				vector<int> vect3 = { foo[0], foo[2] };

				ele.push_back(vect1);
				ele.push_back(vect2);
				ele.push_back(vect3);
				//cout << ele << endl;
				//getchar();
			}


		}

		int freq = 1;

		for (i = 0; i < ele.size(); i++)
		{
			freq = 1;
			for (int j = 0; j < ele.size(); j++)
			{
				if (ele[i] == ele[j])
				{
					freq = freq + 1;
				}
				if (freq >= 3)
				{

					break;
				}
			}
			//cout << ele[i][0] << ele[i][1] << ele[i][2] << endl;
			if (freq == 2)
			{
				ele_boundary.push_back(ele[i]);
			}
		}
	}

	//cout << ele_boundary.size() << endl;
	sort(ele_boundary.begin(), ele_boundary.end());
	ele_boundary.erase(unique(ele_boundary.begin(), ele_boundary.end()), ele_boundary.end());
	//cout << ele_boundary.size() << endl;




	return true;

}


bool PolyCube::IsCornerPoint(int vertexID)
{

	//int i, j, index, indexCluster;
	//int kSize = 0;

	////bool counted = true;

	//vector<int> neiCluster;

	//for (i = 0; i < elementValenceNumber[vertexID]; i++)
	//{

	//	bool counted = true;

	//	index = elementValence[vertexID][i];

	//	indexCluster = pcube_exyzelement[index].indexCluster;

	//	for (j = 0; j < neiCluster.size(); j++)
	//	{

	//		if (neiCluster[j] == indexCluster)
	//		{

	//			counted = false;

	//		}

	//	}

	//	if (counted == true)
	//	{

	//		neiCluster.push_back(indexCluster);

	//		kSize++;

	//	}

	//}

	//if (kSize == 3)
	//{

	//	return true;

	//}
	//else
	//{

	//	return false;

	//}

	//For Bifurcation Cases, but may be also more general
	int i, j, index, indexPatch;
	int kSize = 0;

	vector<int> neiPatch;

	for (i = 0; i < elementValenceNumber[vertexID]; i++)
	{

		bool counted = true;

		index = elementValence[vertexID][i];

		indexPatch = pcube_exyzelement[index].indexPatch;

		for (j = 0; j < neiPatch.size(); j++)
		{

			if (neiPatch[j] == indexPatch)
			{

				counted = false;

			}

		}

		if (counted == true)
		{

			neiPatch.push_back(indexPatch);

			kSize++;

		}

	}

	if (kSize >= 3)
	{

		return true;

	}
	else
	{

		return false;

	}

}

bool PolyCube::Polycube_deformation(string inputFileName, string tempName_1, bool fn_in_patch_flag, string fn_in_patch_k, int loopn, double difference)
{
	string tempName;
	stringstream ss;
	int ifiledeform;

	int fileloop = 0;

	int labelMAX = -1;
	int label;
	int i, index;


	int nn = 3 * elementNumber;
	int ne = elementNumber;

	pcube_split.CreateNewMesh(pcube_split.TRIANGLE, nn, ne);

	pcube_exyzelement.resize(pcube_split.elementNumber);


	//if (READ_NORMAL==1)
	//{
	//	tempName_1 = inputName + "_normal_tri.k";
	ReadTargetNormal(tempName_1.c_str());
	//}
	//else
	//{
	//	tempName = inputName; //vtk
	//	ReadTargetNormal_VTK(tempName.c_str());
	//}

	//if (READ_PATCH ==1)
	//{
	//	tempName_2 = inputName + "_patch_tri.k";
	ReadPatch(tempName_1.c_str());   //give part
	for (i = 0; i < elementNumber; i++)
	{
		GetDirectNeighboringElementByElement(i);
	}


	RefinePatch(); //give real part
				   /*	string tempString = "test_indexPatch_write_mod.k";
				   WriteKFileIndexPatch(tempString.c_str());*/
				   //}
				   //else
				   //{
				   //	tempName = inputName;
				   //	ReadPatch_VTK(tempName.c_str());
				   //}
				   //
				   //getchar();
	SearchCornerandBoundary(cornerPoints);



	//for (i = 0; i < pcube_exyzelement.size(); i++)
	//{
	//	if (pcube_exyzelement[i].indexNormal == 4)
	//	{
	//		for (int j = 0; j < 3; j++)
	//		{

	//			int iVert = element[i][j];
	//			if (vertex[iVert][2] >= 0.984&&vertex[iVert][2] <= 1.10242)
	//			{
	//				vertex[iVert][2] = 1.10242;
	//			}
	//			if (vertex[iVert][2] >= 1.367&&vertex[iVert][2] <= 1.9653)
	//			{
	//				vertex[iVert][2] = 1.9653;
	//			}
	//		}
	//		
	//	}

	//	if (pcube_exyzelement[i].indexNormal == 5)
	//	{
	//		for (int j = 0; j < 3; j++)
	//		{

	//			int iVert = element[i][j];
	//			if (vertex[iVert][2] >= 0.867424&&vertex[iVert][2] <= 0.96)
	//			{
	//				vertex[iVert][2] = 0.867424;
	//			}
	//			if (vertex[iVert][2] >= 0.004576&&vertex[iVert][2] <= 0.605)
	//			{
	//				vertex[iVert][2] = 0.004576;
	//			}
	//		}

	//	}
	//	
	//}

	//cout << "test" << endl;
	//Matrix<int, Dynamic, 2> ele;
	vector<int> index_c;
	vector<vector<double> > distance_move;




	for (int k = 0; k < numberofpatch; k++)
	{

		vector<int>  point_temp;
		for (i = 0; i < pcube_exyzelement.size(); i++)
		{
			if (pcube_exyzelement[i].indexPatch == k)
			{

				//int foo[] = { element[i][0], element[i][1], element[i][2] };
				//int n = sizeof(foo) / sizeof(foo[0]);

				//sort(foo, foo + n);
				//vector<int> vect1 = { foo[0], foo[1] };
				//vector<int> vect2 = { foo[1], foo[2] };
				//vector<int> vect3 = { foo[0], foo[2] };

				point_temp.push_back(element[i][0]);
				point_temp.push_back(element[i][1]);
				point_temp.push_back(element[i][2]);
				//cout << ele << endl;
				//getchar();
			}


		}
		sort(point_temp.begin(), point_temp.end());
		point_temp.erase(unique(point_temp.begin(), point_temp.end()), point_temp.end());
		point_patch.push_back(point_temp);

	}


	//cout << ele.size()<<" "<< ele_boundary.size() << endl;
	//getchar();
	SetBoundaryVertexSign(0);
	for (int loopi_j = 0; loopi_j < ele_boundary.size(); loopi_j++)
	{
		//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
		vertexSign[ele_boundary[loopi_j][0]] = 1;
		vertexSign[ele_boundary[loopi_j][1]] = 1;


	}

	for (int loopi_j = 0; loopi_j < cornerPoints.size(); loopi_j++)
	{
		int index = cornerPoints[loopi_j];
		//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
		vertexSign[index] = 2;


	}


	//for (int loopi_j = 0; loopi_j < vertexNumber; loopi_j++)
	//{			
	//	//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
	//	if (vertexSign[loopi_j]==1)
	//	{
	//		cout << loopi_j << ", ";
	//	}
	//}
	//cout << "test" << endl;
	////cout << endl;
	////cout << endl;
	////getchar();
	//for (int loopi_j = 0; loopi_j < vertexNumber; loopi_j++)
	//{
	//	//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
	//	if (vertexSign[loopi_j] == 2)
	//	{
	//		cout << loopi_j << ", ";
	//	}
	//}




	//getchar();

	//SmoothSurface(2000);
	//cout << ele[0][0] << endl;



	/*tempName = "test_tri.vtk";
	Write(tempName.c_str());*/

	//for (i = 0; i < elementNumber; i++)
	//{
	//	label = pcube_exyzelement[i].indexPatch;

	//	if (label > labelMAX)
	//	{
	//		labelMAX = label;
	//	}
	//}

	//label = labelMAX + 1; // total number of labels
	//NUM_POLYCUBE_PATCH = label; // total number of separated polycube patches


	//cout << NUM_POLYCUBE_PATCH << endl;

	//SearchCornerandBoundary();

	//vector<int> Et;

	//Readtxt(Et);
	//
	//vector<int> Nt;

	//for (int i = 0; i < Et.size(); i++)
	//{
	//	Nt.push_back(element[Et[i]][0]);
	//	Nt.push_back(element[Et[i]][1]);
	//	Nt.push_back(element[Et[i]][2]);
	//}
	//set<int> s(Nt.begin(), Nt.end());
	//Nt.assign(s.begin(), s.end());

	for (int loopi = 0; loopi < loopn; loopi++)
	{

		//if (loopi<7500)
		//{
		//	for (int i = 0; i < Nt.size(); i++)
		//	{
		//		vertex[Nt[i]][0] = vertex[Nt[i]][0] + 5.0 / 3000.0;
		//	}
		//	
		//}
		//Smooth(1);
#pragma region exyzraw





		//if (other_normal_vector[0] == 7 && other_normal_vector[4] == 8)
		//{
		//	if (fileloop % 1890 == 0)
		//	{
		//		for (int elementpatch = 0; elementpatch < 8; i++)
		//		{
		//			SmoothbyPatch(20000, elementpatch);
		//		}

		//		//ss << myio.mesh_write_dir << myio.mesh_name << "_deform_" << setw(5) << setfill('0') << fileloop << "_tri.vtk";
		//		//Write(ss.str().c_str());
		//	}
		//}

		/*	if (loopi == 20000 - 1)
		{
		tempName =  "test_tri.vtk";
		Write(tempName.c_str());
		}*/
		//getchar();

		/*	for (int loopi_test = 0; loopi_test < distance_move.size(); loopi_test++)
		{
		cout << distance_move[loopi_test][0]*20000 << " " << distance_move[loopi_test][1] * 20000 << " " << distance_move[loopi_test][2] * 20000 << " " << endl;
		}
		getchar();*/






		//for (int i = 0; i < ne; i++)
		//{
		//	for (int j = 0; j < 3; j++)
		//	{
		//		pcube_split.element[i][j] = 3 * i + j;
		//		for (int k = 0; k < 3; k++)  //xyz coord
		//		{
		//			pcube_split.vertex[3 * i + j][k] = vertex[element[i][j]][k];
		//		}
		//	}
		//}


#pragma endregion 


		ComputeNormal(2);


#pragma region get_centerpoint

		


		//for (int loopi_normal_i = 0; loopi_normal_i < Normal_edit.size(); loopi_normal_i++)
		//{
		//	for (int loopi_normal_j = 0; loopi_normal_j < Normal_edit[loopi_normal_i].size(); loopi_normal_j++)
		//	{
		//		cout << Normal_edit[loopi_normal_i][loopi_normal_j]<<", ";
		//	}
		//	cout << endl;
		//}

		//

		//getchar();

		vector<vector<double> > center_point_patch;
		for (int k = 0; k < numberofpatch; k++)
		{
			bool normal_edit_flag = false;
			

			vector<double>  center_point_temp;
			center_point_temp.resize(3);
			center_point_temp[0] = 0.0; center_point_temp[1] = 0.0; center_point_temp[2] = 0.0;


			//centerpoint 
			/*for (int i = 0; i < point_patch[k].size(); ++i)
			{
			center_point_temp[0] += vertex[point_patch[k][i]][0];
			center_point_temp[1] += vertex[point_patch[k][i]][1];
			center_point_temp[2] += vertex[point_patch[k][i]][2];
			}

			center_point_temp[0] = center_point_temp[0] /(point_patch[k].size());
			center_point_temp[1] = center_point_temp[1] /(point_patch[k].size());
			center_point_temp[2] = center_point_temp[2] /(point_patch[k].size());
			*/

			switch (relation_patch_part[k][1])
			{
			case 0:
				center_point_temp[0] = -10000000.0;


				for (int i = 0; i < point_patch[k].size(); ++i)
				{
					if (center_point_temp[0] <= vertex[point_patch[k][i]][0])
					{
						center_point_temp[0] = vertex[point_patch[k][i]][0];
					}

					center_point_temp[1] += vertex[point_patch[k][i]][1];
					center_point_temp[2] += vertex[point_patch[k][i]][2];
				}

				if (normal_edit_flag)
				{
					center_point_temp[0] = 0.0;
					for (int i = 0; i < point_patch[k].size(); ++i)
					{
						center_point_temp[0] += vertex[point_patch[k][i]][0];
					}
					center_point_temp[0] = center_point_temp[0] / (point_patch[k].size());
				}


				center_point_temp[1] = center_point_temp[1] / (point_patch[k].size());
				center_point_temp[2] = center_point_temp[2] / (point_patch[k].size());
				break;
			case 1:
				center_point_temp[0] = 10000000.0;
				for (int i = 0; i < point_patch[k].size(); ++i)
				{
					if (center_point_temp[0] >= vertex[point_patch[k][i]][0])
					{
						center_point_temp[0] = vertex[point_patch[k][i]][0];
					}

					center_point_temp[1] += vertex[point_patch[k][i]][1];
					center_point_temp[2] += vertex[point_patch[k][i]][2];
				}

				if (normal_edit_flag)
				{
					center_point_temp[0] = 0.0;
					for (int i = 0; i < point_patch[k].size(); ++i)
					{
						center_point_temp[0] += vertex[point_patch[k][i]][0];
					}
					center_point_temp[0] = center_point_temp[0] / (point_patch[k].size());
				}

				center_point_temp[1] = center_point_temp[1] / (point_patch[k].size());
				center_point_temp[2] = center_point_temp[2] / (point_patch[k].size());
				break;
			case 2:
				center_point_temp[1] = -10000000.0;
				for (int i = 0; i < point_patch[k].size(); ++i)
				{
					if (center_point_temp[1] <= vertex[point_patch[k][i]][1])
					{
						center_point_temp[1] = vertex[point_patch[k][i]][1];
					}

					center_point_temp[0] += vertex[point_patch[k][i]][0];
					center_point_temp[2] += vertex[point_patch[k][i]][2];
				}

				if (normal_edit_flag)
				{
					center_point_temp[1] = 0.0;
					for (int i = 0; i < point_patch[k].size(); ++i)
					{
						center_point_temp[1] += vertex[point_patch[k][i]][1];
					}
					center_point_temp[1] = center_point_temp[1] / (point_patch[k].size());
				}

				center_point_temp[0] = center_point_temp[0] / (point_patch[k].size());
				center_point_temp[2] = center_point_temp[2] / (point_patch[k].size());
				break;
			case 3:
				center_point_temp[1] = 10000000.0;
				for (int i = 0; i < point_patch[k].size(); ++i)
				{
					if (center_point_temp[1] >= vertex[point_patch[k][i]][1])
					{
						center_point_temp[1] = vertex[point_patch[k][i]][1];
					}

					center_point_temp[0] += vertex[point_patch[k][i]][0];
					center_point_temp[2] += vertex[point_patch[k][i]][2];
				}

				if (normal_edit_flag)
				{
					center_point_temp[1] = 0.0;
					for (int i = 0; i < point_patch[k].size(); ++i)
					{
						center_point_temp[1] += vertex[point_patch[k][i]][1];
					}
					center_point_temp[1] = center_point_temp[1] / (point_patch[k].size());
				}

				center_point_temp[0] = center_point_temp[0] / (point_patch[k].size());
				center_point_temp[2] = center_point_temp[2] / (point_patch[k].size());
				break;
			case 4:
				center_point_temp[2] = -10000000.0;
				for (int i = 0; i < point_patch[k].size(); ++i)
				{
					if (center_point_temp[2] <= vertex[point_patch[k][i]][2])
					{
						center_point_temp[2] = vertex[point_patch[k][i]][2];
					}

					center_point_temp[0] += vertex[point_patch[k][i]][0];
					center_point_temp[1] += vertex[point_patch[k][i]][1];
				}

				if (normal_edit_flag)
				{
					center_point_temp[2] = 0.0;
					for (int i = 0; i < point_patch[k].size(); ++i)
					{
						center_point_temp[2] += vertex[point_patch[k][i]][2];
					}
					center_point_temp[2] = center_point_temp[2] / (point_patch[k].size());
				}

				center_point_temp[0] = center_point_temp[0] / (point_patch[k].size());
				center_point_temp[1] = center_point_temp[1] / (point_patch[k].size());
				break;
			case 5:
				center_point_temp[2] = 10000000.0;
				for (int i = 0; i < point_patch[k].size(); ++i)
				{
					if (center_point_temp[2] >= vertex[point_patch[k][i]][2])
					{
						center_point_temp[2] = vertex[point_patch[k][i]][2];
					}

					center_point_temp[0] += vertex[point_patch[k][i]][0];
					center_point_temp[1] += vertex[point_patch[k][i]][1];
				}

				if (normal_edit_flag)
				{
					center_point_temp[2] = 0.0;
					for (int i = 0; i < point_patch[k].size(); ++i)
					{
						center_point_temp[2] += vertex[point_patch[k][i]][2];
					}
					center_point_temp[2] = center_point_temp[2] / (point_patch[k].size());
				}

				center_point_temp[0] = center_point_temp[0] / (point_patch[k].size());
				center_point_temp[1] = center_point_temp[1] / (point_patch[k].size());

				break;
			default:
				for (int i = 0; i < point_patch[k].size(); ++i)
				{
					center_point_temp[0] += vertex[point_patch[k][i]][0];
					center_point_temp[1] += vertex[point_patch[k][i]][1];
					center_point_temp[2] += vertex[point_patch[k][i]][2];
				}
				center_point_temp[0] = center_point_temp[0] / (point_patch[k].size());
				center_point_temp[1] = center_point_temp[1] / (point_patch[k].size());
				center_point_temp[2] = center_point_temp[2] / (point_patch[k].size());

				break;
			}
			/*cout << relation_patch_part[k][1]<<": "<< center_point_temp[0] <<", " << center_point_temp[1] << ", " << center_point_temp[2] << endl;*/

			center_point_patch.push_back(center_point_temp);

		}

		//for (int i = 0; i < center_point_patch.size(); i++)
		//{
		//	
		//	cout<< center_point_patch[i][0]	<<", ";
		//	cout<< center_point_patch[i][1]	<<", ";
		//	cout<< center_point_patch[i][2]	<<endl;


		//}


		for (int k = 0; k < numberofpatch; k++)
		{
			Projection(center_point_patch, k, difference);

		}

		SmoothEdge(10);
		SmoothSurface(10);

#pragma endregion 



#pragma region Generateandimportnormal
		//read target normal

		pcube_split.ComputeNormal(2);

		double tempVec[3];

		for (int i = 0; i < pcube_split.elementNumber; i++)
		{


			for (int j = 0; j < 3; ++j)  //ien
			{
				tempVec[j] = pcube_split.elementNormal[i][j];


			}
			Normalize(tempVec);
			for (int j = 0; j < 3; j++)
			{

				pcube_exyzelement[i].normal[j] = tempVec[j];

			}


		}

		//cout << pcube_split.vertex[0][1] << endl;
#pragma endregion


		//#pragma region exyzraw_mv
		//		int tempien = 0;
		//
		//
		//		//visualize
		//
		//
		//		pcube_split_mv = pcube_split;
		//
		//		for (int i = 0; i < pcube_split.elementNumber; i++)
		//		{
		//			for (int j = 0; j < 3; j++)
		//			{
		//				tempien = pcube_split.element[i][j];
		//				for (int k = 0; k < 3; k++)  //xyz coord
		//				{
		//					pcube_split_mv.vertex[tempien][k] += 0.05*pcube_exyzelement[i].normal[k];
		//				}
		//
		//			}
		//		}
		//
		//
		//		fileloop++;
		//		//if (fileloop < 5)
		//		//{
		//		//	ss.clear();
		//		//	ss.str("");
		//		//	ss << myio.mesh_write_dir << myio.mesh_name << "_whole_" << setw(5) << setfill('0') << fileloop<< "_tri.vtk";
		//		//	
		//		//	pcube_split_mv.Write(ss.str().c_str());
		//		//}
		//		//
		//
		//		
		//#pragma endregion
		//
		//		//for (int i = 0; i < other_normal_vector.size(); i++)
		//		//{
		//		//	cout << other_normal_vector[i] << " ";
		//		//}
		//
		//#pragma region rotate
		//
		//
		//
		//		//center
		//
		//		int index;
		//		for (int i = 0; i < pcube_split.elementNumber; i++)
		//		{
		//			double centeOne[3] = { 0,0,0 };
		//
		//			for (int j = 0; j < 3; j++)
		//			{
		//
		//				index = pcube_split.element[i][j];
		//
		//				centeOne[0] += 1.0 / 3.0 * pcube_split.vertex[index][0];
		//				centeOne[1] += 1.0 / 3.0 * pcube_split.vertex[index][1];
		//				centeOne[2] += 1.0 / 3.0 * pcube_split.vertex[index][2];
		//
		//			}
		//
		//
		//			for (int j = 0; j < 3; j++)
		//			{
		//
		//				pcube_exyzelement[i].center[j] = centeOne[j];
		//
		//			}
		//		}
		//
		//		Vector3d Norm0, Norm1;
		//		Matrix3d Rt;
		//		Vector3d V0, V1;
		//
		//
		//		for (int i = 0; i < pcube_split.elementNumber; i++)
		//		{
		//			Norm0 << pcube_exyzelement[i].normal[0], pcube_exyzelement[i].normal[1], pcube_exyzelement[i].normal[2];
		//			Norm1 << pcube_exyzelement[i].normalNew[0], pcube_exyzelement[i].normalNew[1], pcube_exyzelement[i].normalNew[2];
		//			Rt = rotation_matrix(Norm0, Norm1);
		//			/*cout << Norm0 << endl;
		//			cout << Norm1 << endl;
		//			cout << Rt << endl;*/
		//			for (int j = 0; j < 3; j++)
		//			{
		//				index = pcube_split.element[i][j];
		//				V0 << pcube_split.vertex[index][0] - pcube_exyzelement[i].center[0],
		//					pcube_split.vertex[index][1] - pcube_exyzelement[i].center[1],
		//					pcube_split.vertex[index][2] - pcube_exyzelement[i].center[2];
		//				//cout << V0 << endl;
		//
		//				V1 = Rt*V0;
		//				//cout << V1 << endl;
		//				pcube_split.vertex[index][0] = V1[0] + pcube_exyzelement[i].center[0];
		//				pcube_split.vertex[index][1] = V1[1] + pcube_exyzelement[i].center[1];
		//				pcube_split.vertex[index][2] = V1[2] + pcube_exyzelement[i].center[2];
		//
		//			}
		//		}


		fileloop++;
		//if (fileloop < 5)
		//{
		//	ss.clear();
		//	ss.str("");
		//	ss << myio.mesh_write_dir << myio.mesh_name << "_whole_" <<setw(5) << setfill('0') << fileloop<< "_tri.vtk";
		//	
		//	pcube_split.Write(ss.str().c_str());
		//}




		//tempName = inputName + "_pcube_rotate_tri.vtk";
		//pcube_split.Write(tempName.c_str());
		//tempName = inputName + "_pcube_rotate_tri.raw";
		//pcube_split.Write(tempName.c_str());

		//Matrix3d rotation_matrix(Vector3d vector0, Vector3d vector1)

#pragma endregion

		//#pragma region Smooth
		//
		//		int i, j, k, nVert, iCount, ii;
		//		double center[3];
		//
		//		//find related points
		//		//for (ii = 0; ii<iLoop; ++ii)
		//		for (int i = 0; i < vertexNumber; ++i)
		//		{
		//			center[0] = 0.0; center[1] = 0.0; center[2] = 0.0;
		//			iCount = 0;
		//			for (j = 0; j < elementValenceNumber[i]; ++j)
		//			{
		//
		//				index = elementValence[i][j];
		//				for (k = 0; k < 3; k++)
		//				{
		//
		//					if (element[index][k] == i)
		//					{
		//
		//						break;
		//
		//					}
		//					//iCount = -1; //if no point on the element which makes iCount+1=0 not >0
		//				}
		//
		//				center[0] += pcube_split.vertex[3 * index + k][0];
		//				center[1] += pcube_split.vertex[3 * index + k][1];
		//				center[2] += pcube_split.vertex[3 * index + k][2];
		//				++iCount;
		//			}
		//
		//			if (iCount > 0)
		//			{
		//				vertex[i][0] += (center[0] / iCount - vertex[i][0])*difference;
		//				vertex[i][1] += (center[1] / iCount - vertex[i][1])*difference;
		//				vertex[i][2] += (center[2] / iCount - vertex[i][2])*difference;
		//			}
		//
		//			else
		//			{
		//				cout << "error" << endl;
		//			}
		//		}
		//
		//		//cout << corner_data[6][0] << "test" << endl;
		//		//cout << prin_normal[4][2];
		//		SetBoundaryVertexSign(0);
		//		//for (int loopi_test = 0; loopi_test < index_c.size(); loopi_test++)
		//		//{
		//		//	//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
		//		//	vertexSign[index_c[loopi_test]] = 1;
		//		//	
		//		//}
		//		for (int loopi_j = 0; loopi_j < ele_boundary.size(); loopi_j++)
		//		{
		//			//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
		//			vertexSign[ele_boundary[loopi_j][0]] = 1;
		//			vertexSign[ele_boundary[loopi_j][1]] = 1;
		//		}
		//		//if (other_normal_vector[0] == 7 && other_normal_vector[4] == 8)
		//		//{
		//		//	if (fileloop % 1890 == 0)
		//		//	{
		//		//		for (int elementpatch = 0; elementpatch < 8; i++)
		//		//		{
		//		//			SmoothbyPatch(20000, elementpatch);
		//		//		}
		//
		//		//		//ss << myio.mesh_write_dir << myio.mesh_name << "_deform_" << setw(5) << setfill('0') << fileloop << "_tri.vtk";
		//		//		//Write(ss.str().c_str());
		//		//	}
		//		//}
		//#pragma endregion		
		fileloop++;
		ss.clear();
		ss.str("");
		if (fileloop % 1000 == 0)
		{

			AddProgress((double)loopi / (double)loopn);

			//ss << myio.mesh_write_dir << myio.mesh_name << "_deform_" << setw(5) << setfill('0') << fileloop << "_tri.vtk";
			//Write(ss.str().c_str());
		}

		//cout << ifiledeform;

		//cout << tempName.c_str() << endl;


		//if (fileloop < 5)e
		//{
		//	ss.clear();
		//	ss.str("");
		//	ss << myio.mesh_write_dir << myio.mesh_name << "_whole_" << setw(5) << setfill('0') << fileloop << "_tri.vtk";

		//	Write(ss.str().c_str());
		//}



	}

	//string tempName;
	//if (other_normal_vector[0] == 7 && other_normal_vector[4] == 8)
	//{
	//	SetBoundaryVertexSign(0);
	//	for (int loopi_test = 0; loopi_test < index_c.size(); loopi_test++)
	//	{
	//		//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
	//		vertexSign[index_c[loopi_test]] = 1;

	//	}

	//	for (int i = 0; i < vertexNumber; ++i)
	//	{
	//		if (vertexSign[i] > 0)
	//			cout << i << ",";
	//	}
	//	getchar();
	//	for (int elementpatch = 0; elementpatch < 8; elementpatch++)
	//	{
	//		SmoothbyPatch(20000, elementpatch);
	//	}

	//SetBoundaryVertexSign(0);
	//for (int loopi_j = 0; loopi_j < ele_boundary.size(); loopi_j++)
	//{
	//	//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
	//	vertexSign[ele_boundary[loopi_j][0]] = 1;
	//	vertexSign[ele_boundary[loopi_j][1]] = 1;
	//}
	//SmoothSurface(20);

	//}
	if (fn_in_patch_flag)
	{
		ReadPatch(fn_in_patch_k.c_str());   //give part
		for (i = 0; i < elementNumber; i++)
		{
			GetDirectNeighboringElementByElement(i);
		}


		RefinePatch(); //give real part
					   /*	string tempString = "test_indexPatch_write_mod.k";
					   WriteKFileIndexPatch(tempString.c_str());*/
					   //}
					   //else
					   //{
					   //	tempName = inputName;
					   //	ReadPatch_VTK(tempName.c_str());
					   //}
					   //
					   //getchar();
		SearchCornerandBoundary(cornerPoints);
	}


	return true;

}


bool PolyCube::Polycube_deformation(string inputFileName,  string tempName_1, bool fn_in_patch_flag, string fn_in_patch_k, int loopn, double difference, string normal_projection_edit)
{
	string tempName;
	stringstream ss;
	int ifiledeform;
	
	int fileloop = 0;

	int labelMAX = -1;
	int label;
	int i, index;


	int nn = 3 * elementNumber;
	int ne = elementNumber;

	pcube_split.CreateNewMesh(pcube_split.TRIANGLE, nn, ne);

	pcube_exyzelement.resize(pcube_split.elementNumber);


	//if (READ_NORMAL==1)
	//{
	//	tempName_1 = inputName + "_normal_tri.k";
		ReadTargetNormal(tempName_1.c_str());
	//}
	//else
	//{
	//	tempName = inputName; //vtk
	//	ReadTargetNormal_VTK(tempName.c_str());
	//}

	//if (READ_PATCH ==1)
	//{
	//	tempName_2 = inputName + "_patch_tri.k";
		ReadPatch(tempName_1.c_str());   //give part
		for (i = 0; i < elementNumber; i++)
		{
			GetDirectNeighboringElementByElement(i);
		}
	

		RefinePatch(); //give real part
	/*	string tempString = "test_indexPatch_write_mod.k";
		WriteKFileIndexPatch(tempString.c_str());*/
	//}
	//else
	//{
	//	tempName = inputName;
	//	ReadPatch_VTK(tempName.c_str());
	//}
	//
		//getchar();
		SearchCornerandBoundary(cornerPoints);



		//for (i = 0; i < pcube_exyzelement.size(); i++)
		//{
		//	if (pcube_exyzelement[i].indexNormal == 4)
		//	{
		//		for (int j = 0; j < 3; j++)
		//		{

		//			int iVert = element[i][j];
		//			if (vertex[iVert][2] >= 0.984&&vertex[iVert][2] <= 1.10242)
		//			{
		//				vertex[iVert][2] = 1.10242;
		//			}
		//			if (vertex[iVert][2] >= 1.367&&vertex[iVert][2] <= 1.9653)
		//			{
		//				vertex[iVert][2] = 1.9653;
		//			}
		//		}
		//		
		//	}

		//	if (pcube_exyzelement[i].indexNormal == 5)
		//	{
		//		for (int j = 0; j < 3; j++)
		//		{

		//			int iVert = element[i][j];
		//			if (vertex[iVert][2] >= 0.867424&&vertex[iVert][2] <= 0.96)
		//			{
		//				vertex[iVert][2] = 0.867424;
		//			}
		//			if (vertex[iVert][2] >= 0.004576&&vertex[iVert][2] <= 0.605)
		//			{
		//				vertex[iVert][2] = 0.004576;
		//			}
		//		}

		//	}
		//	
		//}

		//cout << "test" << endl;
		//Matrix<int, Dynamic, 2> ele;
		vector<int> index_c; 
		vector<vector<double> > distance_move;


		

		for (int k = 0; k < numberofpatch; k++)
		{

			vector<int>  point_temp;
			for (i = 0; i < pcube_exyzelement.size(); i++)
			{
				if (pcube_exyzelement[i].indexPatch == k)
				{

					//int foo[] = { element[i][0], element[i][1], element[i][2] };
					//int n = sizeof(foo) / sizeof(foo[0]);

					//sort(foo, foo + n);
					//vector<int> vect1 = { foo[0], foo[1] };
					//vector<int> vect2 = { foo[1], foo[2] };
					//vector<int> vect3 = { foo[0], foo[2] };

					point_temp.push_back(element[i][0]);
					point_temp.push_back(element[i][1]);
					point_temp.push_back(element[i][2]);
					//cout << ele << endl;
					//getchar();
				}


			}
			sort(point_temp.begin(), point_temp.end());
			point_temp.erase(unique(point_temp.begin(), point_temp.end()), point_temp.end());
			point_patch.push_back(point_temp);

		}


		//cout << ele.size()<<" "<< ele_boundary.size() << endl;
		//getchar();
		SetBoundaryVertexSign(0);
		for (int loopi_j = 0; loopi_j < ele_boundary.size(); loopi_j++)
		{
			//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
			vertexSign[ele_boundary[loopi_j][0]] = 1;
			vertexSign[ele_boundary[loopi_j][1]] = 1;


		}
		
		for (int loopi_j = 0; loopi_j < cornerPoints.size(); loopi_j++)
		{
			int index = cornerPoints[loopi_j];
			//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
			vertexSign[index] = 2;


		}


		//for (int loopi_j = 0; loopi_j < vertexNumber; loopi_j++)
		//{			
		//	//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
		//	if (vertexSign[loopi_j]==1)
		//	{
		//		cout << loopi_j << ", ";
		//	}
		//}
		//cout << "test" << endl;
		////cout << endl;
		////cout << endl;
		////getchar();
		//for (int loopi_j = 0; loopi_j < vertexNumber; loopi_j++)
		//{
		//	//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
		//	if (vertexSign[loopi_j] == 2)
		//	{
		//		cout << loopi_j << ", ";
		//	}
		//}




		//getchar();

		//SmoothSurface(2000);
		//cout << ele[0][0] << endl;



		/*tempName = "test_tri.vtk";
		Write(tempName.c_str());*/

	//for (i = 0; i < elementNumber; i++)
	//{
	//	label = pcube_exyzelement[i].indexPatch;

	//	if (label > labelMAX)
	//	{
	//		labelMAX = label;
	//	}
	//}

	//label = labelMAX + 1; // total number of labels
	//NUM_POLYCUBE_PATCH = label; // total number of separated polycube patches


	//cout << NUM_POLYCUBE_PATCH << endl;

	//SearchCornerandBoundary();

	//vector<int> Et;

	//Readtxt(Et);
	//
	//vector<int> Nt;

	//for (int i = 0; i < Et.size(); i++)
	//{
	//	Nt.push_back(element[Et[i]][0]);
	//	Nt.push_back(element[Et[i]][1]);
	//	Nt.push_back(element[Et[i]][2]);
	//}
	//set<int> s(Nt.begin(), Nt.end());
	//Nt.assign(s.begin(), s.end());

	for (int loopi = 0; loopi < loopn; loopi++)
	{

		//if (loopi<7500)
		//{
		//	for (int i = 0; i < Nt.size(); i++)
		//	{
		//		vertex[Nt[i]][0] = vertex[Nt[i]][0] + 5.0 / 3000.0;
		//	}
		//	
		//}
		//Smooth(1);
#pragma region exyzraw

	


		
		//if (other_normal_vector[0] == 7 && other_normal_vector[4] == 8)
		//{
		//	if (fileloop % 1890 == 0)
		//	{
		//		for (int elementpatch = 0; elementpatch < 8; i++)
		//		{
		//			SmoothbyPatch(20000, elementpatch);
		//		}

		//		//ss << myio.mesh_write_dir << myio.mesh_name << "_deform_" << setw(5) << setfill('0') << fileloop << "_tri.vtk";
		//		//Write(ss.str().c_str());
		//	}
		//}

	/*	if (loopi == 20000 - 1)
		{
			tempName =  "test_tri.vtk";
			Write(tempName.c_str());
		}*/
		//getchar();

	/*	for (int loopi_test = 0; loopi_test < distance_move.size(); loopi_test++)
		{
			cout << distance_move[loopi_test][0]*20000 << " " << distance_move[loopi_test][1] * 20000 << " " << distance_move[loopi_test][2] * 20000 << " " << endl;
		}
		getchar();*/






		//for (int i = 0; i < ne; i++)
		//{
		//	for (int j = 0; j < 3; j++)
		//	{
		//		pcube_split.element[i][j] = 3 * i + j;
		//		for (int k = 0; k < 3; k++)  //xyz coord
		//		{
		//			pcube_split.vertex[3 * i + j][k] = vertex[element[i][j]][k];
		//		}
		//	}
		//}
		
		
#pragma endregion 


		ComputeNormal(2);


#pragma region get_centerpoint

		string normal_projection_edit_file=normal_projection_edit+".txt";
		vector<vector<int> > Normal_edit;
		Readtxt(Normal_edit, normal_projection_edit_file);


		//for (int loopi_normal_i = 0; loopi_normal_i < Normal_edit.size(); loopi_normal_i++)
		//{
		//	for (int loopi_normal_j = 0; loopi_normal_j < Normal_edit[loopi_normal_i].size(); loopi_normal_j++)
		//	{
		//		cout << Normal_edit[loopi_normal_i][loopi_normal_j]<<", ";
		//	}
		//	cout << endl;
		//}

		//

		//getchar();

		vector<vector<double> > center_point_patch;
		for (int k = 0; k < numberofpatch; k++)
		{
			bool normal_edit_flag = false;
			for (int loopi_normal_i = 0; loopi_normal_i < Normal_edit.size(); loopi_normal_i++)
			{
				if (k== Normal_edit[loopi_normal_i][0])
				{
					if (1==Normal_edit[loopi_normal_i][1])
					{
						normal_edit_flag = true;
					}
				}
			}


			vector<double>  center_point_temp;
			center_point_temp.resize(3);
			center_point_temp[0] = 0.0; center_point_temp[1] = 0.0; center_point_temp[2] = 0.0;


			//centerpoint 
			/*for (int i = 0; i < point_patch[k].size(); ++i)
			{
				center_point_temp[0] += vertex[point_patch[k][i]][0];
				center_point_temp[1] += vertex[point_patch[k][i]][1];
				center_point_temp[2] += vertex[point_patch[k][i]][2];
			}

			center_point_temp[0] = center_point_temp[0] /(point_patch[k].size());
			center_point_temp[1] = center_point_temp[1] /(point_patch[k].size());
			center_point_temp[2] = center_point_temp[2] /(point_patch[k].size());
			*/

			switch (relation_patch_part[k][1])
			{
			case 0:
				center_point_temp[0] = -10000000.0; 
				

				for (int i = 0; i < point_patch[k].size(); ++i)
				{
					if (center_point_temp[0] <= vertex[point_patch[k][i]][0])
					{
						center_point_temp[0] = vertex[point_patch[k][i]][0];
					}
					
					center_point_temp[1] += vertex[point_patch[k][i]][1];
					center_point_temp[2] += vertex[point_patch[k][i]][2];
				}

				if (normal_edit_flag)
				{
					center_point_temp[0] = 0.0;
					for (int i = 0; i < point_patch[k].size(); ++i)
					{
						center_point_temp[0] += vertex[point_patch[k][i]][0];
					}
					center_point_temp[0] = center_point_temp[0] / (point_patch[k].size());
				}


				center_point_temp[1] = center_point_temp[1] / (point_patch[k].size());
				center_point_temp[2] = center_point_temp[2] / (point_patch[k].size());
				break;
			case 1:
				center_point_temp[0] = 10000000.0;
				for (int i = 0; i < point_patch[k].size(); ++i)
				{
					if (center_point_temp[0] >= vertex[point_patch[k][i]][0])
					{
						center_point_temp[0] = vertex[point_patch[k][i]][0];
					}

					center_point_temp[1] += vertex[point_patch[k][i]][1];
					center_point_temp[2] += vertex[point_patch[k][i]][2];
				}

				if (normal_edit_flag)
				{
					center_point_temp[0] = 0.0;
					for (int i = 0; i < point_patch[k].size(); ++i)
					{
						center_point_temp[0] += vertex[point_patch[k][i]][0];
					}
					center_point_temp[0] = center_point_temp[0] / (point_patch[k].size());
				}

				center_point_temp[1] = center_point_temp[1] / (point_patch[k].size());
				center_point_temp[2] = center_point_temp[2] / (point_patch[k].size());
				break;
			case 2:
				center_point_temp[1] = -10000000.0;
				for (int i = 0; i < point_patch[k].size(); ++i)
				{
					if (center_point_temp[1] <= vertex[point_patch[k][i]][1])
					{
						center_point_temp[1] = vertex[point_patch[k][i]][1];
					}

					center_point_temp[0] += vertex[point_patch[k][i]][0];
					center_point_temp[2] += vertex[point_patch[k][i]][2];
				}

				if (normal_edit_flag)
				{
					center_point_temp[1] = 0.0;
					for (int i = 0; i < point_patch[k].size(); ++i)
					{
						center_point_temp[1] += vertex[point_patch[k][i]][1];
					}
					center_point_temp[1] = center_point_temp[1] / (point_patch[k].size());
				}

				center_point_temp[0] = center_point_temp[0] / (point_patch[k].size());
				center_point_temp[2] = center_point_temp[2] / (point_patch[k].size());
				break;
			case 3:
				center_point_temp[1] = 10000000.0;
				for (int i = 0; i < point_patch[k].size(); ++i)
				{
					if (center_point_temp[1] >= vertex[point_patch[k][i]][1])
					{
						center_point_temp[1] = vertex[point_patch[k][i]][1];
					}

					center_point_temp[0] += vertex[point_patch[k][i]][0];
					center_point_temp[2] += vertex[point_patch[k][i]][2];
				}

				if (normal_edit_flag)
				{
					center_point_temp[1] = 0.0;
					for (int i = 0; i < point_patch[k].size(); ++i)
					{
						center_point_temp[1] += vertex[point_patch[k][i]][1];
					}
					center_point_temp[1] = center_point_temp[1] / (point_patch[k].size());
				}

				center_point_temp[0] = center_point_temp[0] / (point_patch[k].size());
				center_point_temp[2] = center_point_temp[2] / (point_patch[k].size());
				break;
			case 4:
				center_point_temp[2] = -10000000.0;
				for (int i = 0; i < point_patch[k].size(); ++i)
				{
					if (center_point_temp[2] <= vertex[point_patch[k][i]][2])
					{
						center_point_temp[2] = vertex[point_patch[k][i]][2];
					}

					center_point_temp[0] += vertex[point_patch[k][i]][0];
					center_point_temp[1] += vertex[point_patch[k][i]][1];
				}

				if (normal_edit_flag)
				{
					center_point_temp[2] = 0.0;
					for (int i = 0; i < point_patch[k].size(); ++i)
					{
						center_point_temp[2] += vertex[point_patch[k][i]][2];
					}
					center_point_temp[2] = center_point_temp[2] / (point_patch[k].size());
				}

				center_point_temp[0] = center_point_temp[0] / (point_patch[k].size());
				center_point_temp[1] = center_point_temp[1] / (point_patch[k].size());
				break;
			case 5:
				center_point_temp[2] = 10000000.0;
				for (int i = 0; i < point_patch[k].size(); ++i)
				{
					if (center_point_temp[2] >= vertex[point_patch[k][i]][2])
					{
						center_point_temp[2] = vertex[point_patch[k][i]][2];
					}

					center_point_temp[0] += vertex[point_patch[k][i]][0];
					center_point_temp[1] += vertex[point_patch[k][i]][1];
				}

				if (normal_edit_flag)
				{
					center_point_temp[2] = 0.0;
					for (int i = 0; i < point_patch[k].size(); ++i)
					{
						center_point_temp[2] += vertex[point_patch[k][i]][2];
					}
					center_point_temp[2] = center_point_temp[2] / (point_patch[k].size());
				}

				center_point_temp[0] = center_point_temp[0] / (point_patch[k].size());
				center_point_temp[1] = center_point_temp[1] / (point_patch[k].size());

				break;
			default:
				for (int i = 0; i < point_patch[k].size(); ++i)
				{
					center_point_temp[0] += vertex[point_patch[k][i]][0];
					center_point_temp[1] += vertex[point_patch[k][i]][1];
					center_point_temp[2] += vertex[point_patch[k][i]][2];
				}
				center_point_temp[0] = center_point_temp[0] / (point_patch[k].size());
				center_point_temp[1] = center_point_temp[1] / (point_patch[k].size());
				center_point_temp[2] = center_point_temp[2] / (point_patch[k].size());

				break;
			}
			/*cout << relation_patch_part[k][1]<<": "<< center_point_temp[0] <<", " << center_point_temp[1] << ", " << center_point_temp[2] << endl;*/

			center_point_patch.push_back(center_point_temp);

		}

		//for (int i = 0; i < center_point_patch.size(); i++)
		//{
		//	
		//	cout<< center_point_patch[i][0]	<<", ";
		//	cout<< center_point_patch[i][1]	<<", ";
		//	cout<< center_point_patch[i][2]	<<endl;


		//}


		for (int k = 0; k < numberofpatch; k++)
		{
			Projection(center_point_patch, k, difference);

		}

		SmoothEdge(10);
		SmoothSurface(10);

#pragma endregion 



#pragma region Generateandimportnormal
		//read target normal

		pcube_split.ComputeNormal(2);

		double tempVec[3];

		for (int i = 0; i < pcube_split.elementNumber; i++)
		{


			for (int j = 0; j < 3; ++j)  //ien
			{
				tempVec[j] = pcube_split.elementNormal[i][j];


			}
			Normalize(tempVec);
			for (int j = 0; j < 3; j++)
			{

				pcube_exyzelement[i].normal[j] = tempVec[j];

			}


		}

		//cout << pcube_split.vertex[0][1] << endl;
#pragma endregion


//#pragma region exyzraw_mv
//		int tempien = 0;
//
//
//		//visualize
//
//
//		pcube_split_mv = pcube_split;
//
//		for (int i = 0; i < pcube_split.elementNumber; i++)
//		{
//			for (int j = 0; j < 3; j++)
//			{
//				tempien = pcube_split.element[i][j];
//				for (int k = 0; k < 3; k++)  //xyz coord
//				{
//					pcube_split_mv.vertex[tempien][k] += 0.05*pcube_exyzelement[i].normal[k];
//				}
//
//			}
//		}
//
//
//		fileloop++;
//		//if (fileloop < 5)
//		//{
//		//	ss.clear();
//		//	ss.str("");
//		//	ss << myio.mesh_write_dir << myio.mesh_name << "_whole_" << setw(5) << setfill('0') << fileloop<< "_tri.vtk";
//		//	
//		//	pcube_split_mv.Write(ss.str().c_str());
//		//}
//		//
//
//		
//#pragma endregion
//
//		//for (int i = 0; i < other_normal_vector.size(); i++)
//		//{
//		//	cout << other_normal_vector[i] << " ";
//		//}
//
//#pragma region rotate
//
//
//
//		//center
//
//		int index;
//		for (int i = 0; i < pcube_split.elementNumber; i++)
//		{
//			double centeOne[3] = { 0,0,0 };
//
//			for (int j = 0; j < 3; j++)
//			{
//
//				index = pcube_split.element[i][j];
//
//				centeOne[0] += 1.0 / 3.0 * pcube_split.vertex[index][0];
//				centeOne[1] += 1.0 / 3.0 * pcube_split.vertex[index][1];
//				centeOne[2] += 1.0 / 3.0 * pcube_split.vertex[index][2];
//
//			}
//
//
//			for (int j = 0; j < 3; j++)
//			{
//
//				pcube_exyzelement[i].center[j] = centeOne[j];
//
//			}
//		}
//
//		Vector3d Norm0, Norm1;
//		Matrix3d Rt;
//		Vector3d V0, V1;
//
//
//		for (int i = 0; i < pcube_split.elementNumber; i++)
//		{
//			Norm0 << pcube_exyzelement[i].normal[0], pcube_exyzelement[i].normal[1], pcube_exyzelement[i].normal[2];
//			Norm1 << pcube_exyzelement[i].normalNew[0], pcube_exyzelement[i].normalNew[1], pcube_exyzelement[i].normalNew[2];
//			Rt = rotation_matrix(Norm0, Norm1);
//			/*cout << Norm0 << endl;
//			cout << Norm1 << endl;
//			cout << Rt << endl;*/
//			for (int j = 0; j < 3; j++)
//			{
//				index = pcube_split.element[i][j];
//				V0 << pcube_split.vertex[index][0] - pcube_exyzelement[i].center[0],
//					pcube_split.vertex[index][1] - pcube_exyzelement[i].center[1],
//					pcube_split.vertex[index][2] - pcube_exyzelement[i].center[2];
//				//cout << V0 << endl;
//
//				V1 = Rt*V0;
//				//cout << V1 << endl;
//				pcube_split.vertex[index][0] = V1[0] + pcube_exyzelement[i].center[0];
//				pcube_split.vertex[index][1] = V1[1] + pcube_exyzelement[i].center[1];
//				pcube_split.vertex[index][2] = V1[2] + pcube_exyzelement[i].center[2];
//
//			}
//		}


		fileloop++;
		//if (fileloop < 5)
		//{
		//	ss.clear();
		//	ss.str("");
		//	ss << myio.mesh_write_dir << myio.mesh_name << "_whole_" <<setw(5) << setfill('0') << fileloop<< "_tri.vtk";
		//	
		//	pcube_split.Write(ss.str().c_str());
		//}
	
		


		//tempName = inputName + "_pcube_rotate_tri.vtk";
		//pcube_split.Write(tempName.c_str());
		//tempName = inputName + "_pcube_rotate_tri.raw";
		//pcube_split.Write(tempName.c_str());

		//Matrix3d rotation_matrix(Vector3d vector0, Vector3d vector1)

#pragma endregion

//#pragma region Smooth
//
//		int i, j, k, nVert, iCount, ii;
//		double center[3];
//
//		//find related points
//		//for (ii = 0; ii<iLoop; ++ii)
//		for (int i = 0; i < vertexNumber; ++i)
//		{
//			center[0] = 0.0; center[1] = 0.0; center[2] = 0.0;
//			iCount = 0;
//			for (j = 0; j < elementValenceNumber[i]; ++j)
//			{
//
//				index = elementValence[i][j];
//				for (k = 0; k < 3; k++)
//				{
//
//					if (element[index][k] == i)
//					{
//
//						break;
//
//					}
//					//iCount = -1; //if no point on the element which makes iCount+1=0 not >0
//				}
//
//				center[0] += pcube_split.vertex[3 * index + k][0];
//				center[1] += pcube_split.vertex[3 * index + k][1];
//				center[2] += pcube_split.vertex[3 * index + k][2];
//				++iCount;
//			}
//
//			if (iCount > 0)
//			{
//				vertex[i][0] += (center[0] / iCount - vertex[i][0])*difference;
//				vertex[i][1] += (center[1] / iCount - vertex[i][1])*difference;
//				vertex[i][2] += (center[2] / iCount - vertex[i][2])*difference;
//			}
//
//			else
//			{
//				cout << "error" << endl;
//			}
//		}
//
//		//cout << corner_data[6][0] << "test" << endl;
//		//cout << prin_normal[4][2];
//		SetBoundaryVertexSign(0);
//		//for (int loopi_test = 0; loopi_test < index_c.size(); loopi_test++)
//		//{
//		//	//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
//		//	vertexSign[index_c[loopi_test]] = 1;
//		//	
//		//}
//		for (int loopi_j = 0; loopi_j < ele_boundary.size(); loopi_j++)
//		{
//			//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
//			vertexSign[ele_boundary[loopi_j][0]] = 1;
//			vertexSign[ele_boundary[loopi_j][1]] = 1;
//		}
//		//if (other_normal_vector[0] == 7 && other_normal_vector[4] == 8)
//		//{
//		//	if (fileloop % 1890 == 0)
//		//	{
//		//		for (int elementpatch = 0; elementpatch < 8; i++)
//		//		{
//		//			SmoothbyPatch(20000, elementpatch);
//		//		}
//
//		//		//ss << myio.mesh_write_dir << myio.mesh_name << "_deform_" << setw(5) << setfill('0') << fileloop << "_tri.vtk";
//		//		//Write(ss.str().c_str());
//		//	}
//		//}
//#pragma endregion		
		fileloop++;
		ss.clear();
		ss.str("");
		if (fileloop % 1000 == 0)
		{
			
			AddProgress((double)loopi / (double)loopn);

			//ss << myio.mesh_write_dir << myio.mesh_name << "_deform_" << setw(5) << setfill('0') << fileloop << "_tri.vtk";
			//Write(ss.str().c_str());
		}
		
		//cout << ifiledeform;

		//cout << tempName.c_str() << endl;
		
	
		//if (fileloop < 5)e
		//{
		//	ss.clear();
		//	ss.str("");
		//	ss << myio.mesh_write_dir << myio.mesh_name << "_whole_" << setw(5) << setfill('0') << fileloop << "_tri.vtk";

		//	Write(ss.str().c_str());
		//}
		


	}

	//string tempName;
	//if (other_normal_vector[0] == 7 && other_normal_vector[4] == 8)
	//{
	//	SetBoundaryVertexSign(0);
	//	for (int loopi_test = 0; loopi_test < index_c.size(); loopi_test++)
	//	{
	//		//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
	//		vertexSign[index_c[loopi_test]] = 1;

	//	}

	//	for (int i = 0; i < vertexNumber; ++i)
	//	{
	//		if (vertexSign[i] > 0)
	//			cout << i << ",";
	//	}
	//	getchar();
	//	for (int elementpatch = 0; elementpatch < 8; elementpatch++)
	//	{
	//		SmoothbyPatch(20000, elementpatch);
	//	}

	//SetBoundaryVertexSign(0);
	//for (int loopi_j = 0; loopi_j < ele_boundary.size(); loopi_j++)
	//{
	//	//cout << ele_boundary[loopi][0] << "," << ele_boundary[loopi][1] << "," ;
	//	vertexSign[ele_boundary[loopi_j][0]] = 1;
	//	vertexSign[ele_boundary[loopi_j][1]] = 1;
	//}
	//SmoothSurface(20);

	//}
	if (fn_in_patch_flag)
	{
		ReadPatch(fn_in_patch_k.c_str());   //give part
		for (i = 0; i < elementNumber; i++)
		{
			GetDirectNeighboringElementByElement(i);
		}


		RefinePatch(); //give real part
					   /*	string tempString = "test_indexPatch_write_mod.k";
					   WriteKFileIndexPatch(tempString.c_str());*/
					   //}
					   //else
					   //{
					   //	tempName = inputName;
					   //	ReadPatch_VTK(tempName.c_str());
					   //}
					   //
					   //getchar();
		SearchCornerandBoundary(cornerPoints);
	}
	

	return true;

}


bool PolyCube::ReadTargetNormal(const char *inputName)
{
	int i, j, tempElementSign;
	string oneLine;

	string inputFileName = inputName;
	fstream input(inputFileName, fstream::in);
	if (!input)
	{
		cout << "open input file error!!!" << endl;
		return false;
	}

	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);

	for (i = 0; i < pcube_exyzelement.size(); i++)
	{
		input >> j >> tempElementSign >> j >> j >> j >> j >> j >> j >> j >> j;
		pcube_exyzelement[i].indexNormal = tempElementSign - 1;
	}

	input.close();

	return true;

}

bool PolyCube::ReadPatch(const char *inputName)
{
	int i, j, tempElementSign;
	string oneLine;

	string inputFileName = inputName;
	fstream input(inputFileName, fstream::in);
	if (!input)
	{
		cout << "open input file error!!!" << endl;
		return false;
	}

	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);
	getline(input, oneLine);

	for (i = 0; i < elementNumber; i++)
	{
		input >> j >> tempElementSign >> j >> j >> j >> j >> j >> j >> j >> j;
		pcube_exyzelement[i].indexPart = tempElementSign - 1;
	}

	input.close();

	return true;

}

bool PolyCube::RefinePatch(void)
{

		int i, j, k, index;
		int label = 0;

		vector<vector<int> > tempCluster(elementNumber);

		for (i = 0; i < elementNumber; i++)
		{

			tempCluster[i].push_back(-1);
			tempCluster[i].push_back(pcube_exyzelement[i].indexPart);

		}

		for (i = 0; i < elementNumber; i++)
		{

			if (0 > tempCluster[i][0])
			{

				tempCluster[i][0] = label;

				int count = 1;

				vector<int> tempIndex;

				tempIndex.push_back(i);
				for (int c = 0; c < count; c++)
				{

					for (int n = 0; n < pcube_exyzelement[tempIndex[c]].numDirectNei; n++)
					{

						index = pcube_exyzelement[tempIndex[c]].directNei[n];

						if (0 > tempCluster[index][0] && pcube_exyzelement[i].indexPart == pcube_exyzelement[index].indexPart
							)//&& pcube_exyzelement[i].indexPart == pcube_exyzelement[index].indexPart) //For bifurcation
						{

							tempCluster[index][0] = label;

							tempIndex.push_back(index);

							count++;

						}

					}

				}

				label++;

			}

		}

		numberofpatch = label; // total number of separated polycube patches
		//cout << "test"<<label << endl;

		for (i = 0; i < elementNumber; i++)
		{

			index = tempCluster[i][0];
			//polycubePatch[index].element.push_back(i);
			//polycubePatch[index].numElements++;	
			pcube_exyzelement[i].indexPatch = index;


		}

		sort(tempCluster.begin(), tempCluster.end());
		tempCluster.erase(unique(tempCluster.begin(), tempCluster.end()), tempCluster.end());

		relation_patch_part = tempCluster;




		//for (i = 0; i < NUM_POLYCUBE_PATCH; i++)
		//{

		//	index = polycubePatch[i].element[0];
		//	polycubePatch[i].indexCluster = tempCluster[index][1];

		//}

		//OutputPatchesVTKBifurcation("test.vtk");

		//if (READ_WRITE_K_INDEXPATCH == 0)
		//{
		//	string tempString;
		//	tempString = inputName + "_indexPatch_write.k";
		//	WriteKFileIndexPatch(tempString.c_str());
		//}

	

	return true;
}


bool PolyCube::WriteKFileIndexPatch(const char *outputName)
{


	int i, j;

	string outputFileName = outputName;

	fstream output(outputFileName, fstream::out);
	if (!output)
	{
		cout << "open input file error!!!" << endl;
		return false;
	}

	output << "$# LS-DYNA Keyword file created by LS-PrePost 4.2 (Beta)" << "\n";
	output << "$# Created on Jun-10-2015" << "\n";

	output << "*KEYWORD" << "\n";
	output << "*ELEMENT_SHELL" << "\n";
	for (i = 0; i < elementNumber; i++)
	{
		output << i + 1 << "," << pcube_exyzelement[i].indexPatch + 1 << "," << element[i][0] + 1 << "," << element[i][1] + 1 << "," << element[i][2] + 1 << "," << element[i][2] + 1 << "\n";
	}

	output << "*NODE" << "\n";
	for (i = 0; i < vertexNumber; i++)
	{
		output << i + 1 << "," << vertex[i][0] << "," << vertex[i][1] << "," << vertex[i][2] << "\n";
	}

	output << "*END" << "\n";

	output.close();
	return true;

}


bool PolyCube::GetDirectNeighboringElementByElement(int elementID)
{

	int i, iVert, j, jElem, k, kSize;

	vector<int> tempElements;

	for (i = 0; i < 3; i++)
	{

		iVert = element[elementID][i];

		for (j = 0; j < elementValenceNumber[iVert]; j++)
		{

			jElem = elementValence[iVert][j];

			kSize = tempElements.size();

			for (k = 0; k < kSize; k++)
			{

				if (jElem == elementID || jElem == tempElements[k])
				{

					break;

				}

			}
			if (k == kSize && jElem != elementID) //Very important when kSize = 0
			{

				tempElements.push_back(jElem);

			}

		}

	}

	//Find three direct neighboring elements
	kSize = 0;

	for (i = 0; i < tempElements.size(); i++)
	{

		for (j = 0; j < 3; j++)
		{

			iVert = element[tempElements[i]][j];

			for (k = 0; k < 3; k++)
			{

				if (iVert == element[elementID][k])
				{
					kSize++;
				}

			}

		}

		if (kSize == 2)
		{

			pcube_exyzelement[elementID].directNei.push_back(tempElements[i]);

		}

		kSize = 0;

	}

	pcube_exyzelement[elementID].numDirectNei = pcube_exyzelement[elementID].directNei.size();

	return true;

}



bool PolyCube::ReadTargetNormal_VTK(const char *inputName)
{
	string inputFileName = inputName;
	string fname(inputFileName + ".vtk"), stmp;
	int npts, neles, itmp,ntype,ndata;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		
		for (int i = 0; i<npts; i++)
		{
			getline(fin, stmp);
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		for (int i = 0; i<neles; i++)
		{
			getline(fin, stmp);
		}
		getline(fin, stmp);
		fin >> stmp >> ntype;
		for (int i = 0; i<ntype; i++)
		{
			getline(fin, stmp);
		}
		getline(fin, stmp);
		fin >> stmp >> ndata;
		getline(fin, stmp);
		getline(fin, stmp);
		getline(fin, stmp);
		for (int i = 0; i<ntype; i++)
		{
			fin >> pcube_exyzelement[i].indexNormal;
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}

	

	return true;

}

bool PolyCube::ReadPatch_VTK(const char *inputName)
{
	string inputFileName = inputName;
	string fname(inputFileName + ".vtk"), stmp;
	int npts, neles, itmp, ntype, ndata,data_temp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;

		for (int i = 0; i<npts; i++)
		{
			getline(fin, stmp);
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		for (int i = 0; i<neles; i++)
		{
			getline(fin, stmp);
		}
		getline(fin, stmp);
		fin >> stmp >> ntype;
		for (int i = 0; i<ntype; i++)
		{
			getline(fin, stmp);
		}
		getline(fin, stmp);
		fin >> stmp >> ndata;
		getline(fin, stmp);
		getline(fin, stmp);
		getline(fin, stmp);
		for (int i = 0; i<ntype; i++)
		{
			
			fin >> pcube_exyzelement[i].indexPart;
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}



	return true;

}

Matrix3d PolyCube::rotation_matrix(Vector3d vector0, Vector3d vector1)
{
	return Matrix3d(Quaterniond::FromTwoVectors(vector0, vector1));
}


bool PolyCube::GetRelativeCoordinate(void)
{
	for (int i = 0; i < elementNumber; ++i) {

		for (int j = 0; j < 3; ++j) {


			int n2 = element[i][j];
			int n0 = element[i][(j + 1) % 3];
			int n1 = element[i][(j + 2) % 3];

			Vector3d v0 ;
			Vector3d v1 ;
			Vector3d v2 ;

			v0 << vertex[n0][0], vertex[n0][1], vertex[n0][2];
			v1 << vertex[n1][0],vertex[n1][1],vertex[n1][2];
			v2 << vertex[n2][0],vertex[n2][1],vertex[n2][2];
		
			Vector3d v01(v1 - v0);
			Vector3d v01N(v01);
			
			v01N.normalize();
			



			/*Vector3d v01Rot90(v01(1), -v01(0));
			Vector3d v01Rot90N(v01Rot90);
			v01Rot90N.normalize();
*/

			Vector3d v02(v2 - v0);

			double fX = v02.dot(v01N);

			Vector3d v01Rot90(v02 - fX*v01N);

			Vector3d v01Rot90N(v01Rot90);

			v01Rot90N.normalize();

			double fY = v02.dot(v01Rot90N);


			Vector3d v2test(v0 + fX * v01 + fY * v01Rot90);
			Vector3d fLength(v2test - v2);
			//cout << "compare"<<fLength.norm() << endl;
			////t.vTriCoords[j] = Vector2d(fX, fY);
	
		}
	}
	return true;
}


bool PolyCube::CreateVertexID(string inputFileName)
{
	Polycube_initial_simple(inputFileName);
	double tempMin = 10000000.0f;
	double tempMax = -10000000.0f;
	vector <int> tempvertexid;
	int tempaxis, tempcount1,tempcount2;
	for (int j = 0; j < 3; j++)
	{
		tempMin = 10000000.0f;
		tempMax = -10000000.0f;
		tempaxis = j;
		for (int i = 0; i < vertexNumber; ++i)
		{		
			if (vertex[i][tempaxis] < tempMin)
			{
				tempMin = vertex[i][tempaxis];
				tempcount1 = i;
			}

			if (vertex[i][tempaxis] > tempMax)
			{
				tempMax = vertex[i][tempaxis];
				tempcount2 = i;
			}
		}

		tempvertexid.push_back(tempcount1);
		tempvertexid.push_back(tempcount2);
	}

	ofstream outputVertexID;
	outputVertexID.open(inputName + "_Output_MaxMinVertexID.txt");

	for (int i = 0; i < tempvertexid.size(); i++)
	{
		outputVertexID << tempvertexid[i]<< ",";
	}
	outputVertexID << endl;
	outputVertexID.close();
	return true;
}


bool AddProgress(double progress)
{
	int barWidth = 70;

	std::cout << "[";
	int pos = barWidth * progress;
	for (int i = 0; i < barWidth; ++i) {
		if (i < pos) std::cout << "=";
		else if (i == pos) std::cout << ">";
		else std::cout << " ";
	}
	std::cout << "] " << int(progress * 100.0) << " %\r";
	std::cout.flush();

	//	progress += 0.16; // for demonstration only

	std::cout << std::endl;
	return true;
}

void PolyCube::Projection(const vector<vector<double>>& center_point_patch, int patch_number,double step)
{
	//for (int i = 0; i < pcube_split.elementNumber; i++)
	//{
		double normVec2[3] = { 0,0,0 };

		switch (relation_patch_part[patch_number][1])
		{
		case 0:
			normVec2[0] = 1;
			break;
		case 1:
			normVec2[0] = -1;
			break;
		case 2:
			normVec2[1] = 1;
			break;
		case 3:
			normVec2[1] = -1;
			break;
		case 4:
			normVec2[2] = 1;
			break;
		case 5:
			normVec2[2] = -1;
			break;
		case 6:
			normVec2[0] = -1.0 / sqrt(2);
			normVec2[1] =0.0 / sqrt(2);
			normVec2[2] = 1.0 / sqrt(2);

			break;
		case 7:
			normVec2[0] = -1.0 / sqrt(2);
			normVec2[1] = 0 / sqrt(2);
			normVec2[2] = -1.0 / sqrt(2);
			//cout << i << normVec2[0] << " "<< normVec2[1] << " "<<normVec2[2] << endl;
			//getchar();
			break;
		default:
			break;
		}
	/*	for (int j = 0; j < 3; j++)
		{

			pcube_exyzelement[i].normalNew[j] = tempVec2[j];

		}*/
	//}

	//for (int i = 0; i < relation_patch_part.size(); i++)
	//{
	//	cout << relation_patch_part[i][0] << ", ";
	//	cout << relation_patch_part[i][1] << endl;
	//}

	//getchar();

		//cout << step << endl;

	for (int i = 0; i < point_patch[patch_number].size(); ++i)
	{
		double v2p[3] = { 0,0,0 };
		double dist;
		int index = point_patch[patch_number][i];
		center_point_patch;
		v2p[0] = vertex[index][0] - center_point_patch[patch_number][0];
		v2p[1] = vertex[index][1] - center_point_patch[patch_number][1];
		v2p[2] = vertex[index][2] - center_point_patch[patch_number][2];
		dist = v2p[0] * normVec2[0] + v2p[1] * normVec2[1] +v2p[2] * normVec2[2];
		vertex[index][0] = vertex[index][0]	- dist * normVec2[0] * step;
		vertex[index][1] = vertex[index][1] - dist * normVec2[1]  *step;
		vertex[index][2] = vertex[index][2] - dist * normVec2[2]  *step;
	}
}

void PolyCube::SmoothbyPatch(int iLoop, int elementpatch)
{

	int i, j, k, nVert, iCount, ii;
	double center[3];

	PushVertexSign();
	InitiateEdgeValence();
	int index;
	//find related points
	for (ii = 0; ii<iLoop; ++ii)
	{
		for (int i = 0; i < vertexNumber; ++i)
		{
			if (vertexSign[i] > 0)
				continue;
			center[0] = 0.0; center[1] = 0.0; center[2] = 0.0;
			iCount = 0;
			vector <int> neipoint;
			for (j = 0; j < elementValenceNumber[i]; ++j)
			{

				index = elementValence[i][j];

				//cout << index << endl;
				if (pcube_exyzelement[index].indexNormal != elementpatch)
				{

					continue;

				}
				neipoint.push_back(element[index][0]);
				neipoint.push_back(element[index][1]);
				neipoint.push_back(element[index][2]);
			}
				//iCount = -1; //if no point on the element which makes iCount+1=0 not >0
			sort(neipoint.begin(), neipoint.end());
			neipoint.erase(unique(neipoint.begin(), neipoint.end()), neipoint.end());
			//cout << "test:" << i << ":   ";
			for (int loop_sm = 0; loop_sm < neipoint.size(); loop_sm++)
			{

				if (neipoint[loop_sm]==i)
				{
					continue;
				}
				//cout << neipoint[loop_sm] << ",";
				center[0] += vertex[neipoint[loop_sm]][0];
				center[1] += vertex[neipoint[loop_sm]][1];
				center[2] += vertex[neipoint[loop_sm]][2];
				++iCount;
			}
			//cout << endl;

			//getchar();

			if (iCount > 0)
			{
				vertex[i][0] += (center[0] / iCount - vertex[i][0])*0.1;
				vertex[i][1] += (center[1] / iCount - vertex[i][1])*0.1;
				vertex[i][2] += (center[2] / iCount - vertex[i][2])*0.1;
			}

			//else
			//{
			//	cout << "Nothing error2" << endl;
			//}
		}
	}
	
}


void PolyCube::SmoothEdge(int iLoop)
{
	int i, j, nVert, iCount, ii;
	double center[3];

	PushVertexSign();
	InitiateEdgeValence();

	for (ii = 0; ii<iLoop; ++ii)
		for (i = 0; i<vertexNumber; ++i)
		{
			if (vertexSign[i] != 1)
				continue;

			center[0] = 0.0; center[1] = 0.0; center[2] = 0.0;
			iCount = 0;
			for (j = 0; j<ele_boundary.size(); ++j)
			{
				if (ele_boundary[j][0]==i)
				{
					center[0] += vertex[ele_boundary[j][1]][0];
					center[1] += vertex[ele_boundary[j][1]][1];
					center[2] += vertex[ele_boundary[j][1]][2];
					++iCount;
				}
				if (ele_boundary[j][1] == i)
				{
					center[0] += vertex[ele_boundary[j][0]][0];
					center[1] += vertex[ele_boundary[j][0]][1];
					center[2] += vertex[ele_boundary[j][0]][2];
					++iCount;
				}
				
			}

			if (iCount ==2)
			{
				vertex[i][0] += (center[0] / iCount - vertex[i][0])*0.02;
				vertex[i][1] += (center[1] / iCount - vertex[i][1])*0.02;
				vertex[i][2] += (center[2] / iCount - vertex[i][2])*0.02;
			}
			else
			{
				/*cout << "error edge:" << i << endl;*/
			}
		}

	PopVertexSign();
	return;
}

