#ifndef _MYIO_H_
#define _MYIO_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <direct.h> 

#define DIRECTORY	".\\"

using namespace std;

class Iodir
{
public:
	
	string mesh_name;
	string mesh_read_dir;
	string mesh_write_dir;

	Iodir();
	~Iodir();

	bool	GenerateDir(string filename);

protected:
	
	
private:

	string	ExtractMeshname(string filename);
	bool	mkdir();
};
#endif