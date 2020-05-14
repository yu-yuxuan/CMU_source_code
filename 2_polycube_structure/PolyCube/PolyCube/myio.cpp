#include "myio.h"

Iodir::Iodir()
{

}

Iodir::~Iodir()
{

}

bool	Iodir::GenerateDir(string filename)
{
	if (filename.length() < 5)
	{
		cerr<<"Invalid mesh file name!!!" <<endl;
		return false;
	}
	mesh_name = ExtractMeshname(filename);
	mkdir();
	
	return true;
}


string	Iodir::ExtractMeshname(string filename)
{
	string s = filename;
	mesh_read_dir = DIRECTORY;
	s.erase(s.size()-4, s.size());
	//if (s.find("_")) s.erase(s.size()-4, s.size());

	return s;
}
bool	Iodir::mkdir()
{
	mesh_write_dir = DIRECTORY;
	mesh_write_dir.append(mesh_name).append("\\");
	_mkdir(mesh_write_dir.c_str());

	return true;

}
 