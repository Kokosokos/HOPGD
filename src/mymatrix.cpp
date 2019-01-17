#include "mymatrix.h"

SimpleMatrix2D::SimpleMatrix2D()
{
	dim=2;
}
//---------------------------------------------------------------------------------------------------
SimpleMatrix2D::SimpleMatrix2D(int in_size1, int in_size2)
{
	setSize(in_size1, in_size2);
}
//---------------------------------------------------------------------------------------------------
void SimpleMatrix2D::setSize(int in_size1, int in_size2)
{
	dim=2;
	size = new int[dim];
	size[0]=in_size1;
	size[1]=in_size2;
	value=new double[size[0]*size[1]];
//#include <iostream>
//	if (value == nullptr)
//	{
//	    cout << "Error: memory could not be allocated";
//	}

}
//---------------------------------------------------------------------------------------------------
int SimpleMatrix2D::columns()
{
	return size[0];
}
//---------------------------------------------------------------------------------------------------
int SimpleMatrix2D::rows()
{
	return size[1];
}
//---------------------------------------------------------------------------------------------------
double SimpleMatrix2D::getElement(int i, int j)
{
	return value[i*size[0]+j];
}
//---------------------------------------------------------------------------------------------------
void SimpleMatrix2D::setElement(int i, int j, double val)
{
	value[i*size[0]+j]=val;
}
//---------------------------------------------------------------------------------------------------
SimpleMatrix2D::~SimpleMatrix2D()
{
	delete[] size;
	delete[] value;
}
//---------------------------------------------------------------------------------------------------
