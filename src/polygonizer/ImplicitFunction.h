#ifndef IMPLICIT_FUNCTION_H
#define IMPLICIT_FUNCTION_H


class ImplicitFunction  
{
public:
	virtual void asignValueToVoxels(float ***&v, bool ***&isIn, float o[], float space[], int dim[]);
	virtual float value(float x, float y, float z, bool &isValid);
	virtual void gradient(float g[3], float x, float y, float z);
	virtual float value(float x, float y, float z);
	ImplicitFunction();
	virtual ~ImplicitFunction();

};

#endif
