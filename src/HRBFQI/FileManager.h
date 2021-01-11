#ifndef FILE_MANAGER_H
#define FILE_MANAGER_H


#include "../DataStructure/PolygonalMesh.h"
#include "../DataStructure/PointSet.h"
#include "stdio.h"

#include <string>

class FileManager  
{
protected:
	void saveWaveFront(PolygonalMesh *mesh, float oriCt[3], float scale = 1.0f);

	FILE* file;
	std::string file_name;
	std::string file_ext;

public:
	void save(PolygonalMesh *mesh, float oriCt[3], float scale = 1.0f);
	void setFile(FILE* file, std::string file_name, std::string file_ext);
	void open(PointSet* ps);
	FileManager();
	virtual ~FileManager();

};


#endif