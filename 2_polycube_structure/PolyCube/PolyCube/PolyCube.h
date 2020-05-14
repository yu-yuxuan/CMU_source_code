#pragma once

#include <vector>
#include <fstream>
#include "rawmesh.h"
#include "myio.h"
#include "StaticVars.h"
#include "Eigen/Dense"
#include "Eigen/LU"
#include "Eigen/SparseCore"
#include "Eigen/SparseCholesky"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/Geometry"


using namespace Eigen;
using namespace std;

const double EPSILON = 1e-6; // 1e-8 is too small sometimes


class PolyCube :
	public RawMesh
{
public:


	Iodir myio;
	PolyCube(void);
	PolyCube(string inputFileName);

	~PolyCube(void);

	bool ReadManualFile(string fn_input_polycube_k, string fn_corner_out);

	static const int NUM_CLUSTER = 6;
	static const int NUM_NEI_CLUSTER = 12;
	static const int READ_PATCH = 1;
	static const int READ_NORMAL = 1;
	static const int INPUTFILE_RAW = 1;

	int NUM_POLYCUBE_PATCH;

	string inputName;


	struct Element_Stru
	{
		int indexNormal;
		int index;
		//int ne;
		int indexPatch;
		int indexPart;
		double normal[3];
		double normalNew[3];
		double center[3];
		vector<int> directNei;
		int numDirectNei;
		
	};
	int numberofpatch;
	vector<vector<int> > ele_boundary;
	vector<vector<int> > relation_patch_part;
	vector<int> cornerPoints;
	//vector<Element_Stru> elementArray;
	vector<vector<int> > point_patch;
	vector<Element_Stru> pcube_exyzelement; //same with pcube_split_exyzelement
	Matrix3d rotation_matrix(Vector3d vector0, Vector3d vector1);

	RawMesh pcube_split;
	RawMesh pcube_split_mv;
	bool GetIntGeometry(string fn_corner_out);
	double round(double N);
	// for visualization
	bool Polycube_initial_simple(string inputFileName);
	void SmoothEdge(int iLoop);
	bool Readtxt(vector<vector<int> >& Nt, string tempName);
	bool SearchCornerandBoundary(vector<int>& cornerPoints);
	bool IsCornerPoint(int vertexID);
	bool Polycube_deformation(string inputFileName, string tempName_1, bool fn_in_patch_flag, string fn_in_patch_k, int loopn, double difference);
	//bool Polycube_deformation(string inputFileName, string tempName_1,  int loopn, double difference, string normal_projection_edit);
	bool Polycube_deformation(string inputFileName, string tempName_1, bool fn_in_patch_flag, string fn_in_patch_k, int loopn, double difference, string normal_projection_edit);
	//bool Polycube_deformation(string inputFileName,int loopn, double difference);
	bool ReadTargetNormal(const char *inputName);
	bool ReadPatch(const char *inputName);
	bool RefinePatch(void);
	bool WriteKFileIndexPatch(const char * outputName);
	bool GetDirectNeighboringElementByElement(int elementID);
	bool ReadTargetNormal_VTK(const char *inputName);
	bool ReadPatch_VTK(const char *inputName);


	bool CreateVertexID(string inputFileName);
	void Projection(const vector<vector<double>>& center, int patch_number, double step);
	void SmoothbyPatch(int iLoop, int elementpatch);
	//bool EnforceLabelConnectivity();
	//bool CheckLabelConnectivity();
	bool GetRelativeCoordinate(void);
	
protected:

	stringstream ss;
	fstream fs;


private:
	/*bool SearchCornerandBoundary(void);
	bool IsCornerPoint(int vertexID);
	bool IsBoundaryPoint(int vertexID);
	bool IsBoundaryEdge(int vertexIDone, int vertexIDtwo);*/

private:


};

bool AddProgress(double progress);