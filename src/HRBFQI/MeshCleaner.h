#ifndef MESH_CLEANER_H
#define MESH_CLEANER_H


#include "../DataStructure/PolygonalMesh.h"
#include "HRBF.h"

class CMeshCleaner
{
public:
	CMeshCleaner(void);
	virtual ~CMeshCleaner(void);

	PolygonalMesh * removeLowConfidenceGeometry(PolygonalMesh *oriMesh, HRBF *func, float confThreshold);
	PolygonalMesh * removeSmallIsolatedComponent(PolygonalMesh *oriMesh, int componentSize);
};


#endif // !MESH_CLEANER_H