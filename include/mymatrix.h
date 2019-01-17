/*
 * matrix.h
 *
 *  Created on: 6 Nov 2018
 *      Author: ikriuchevs
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <blaze/Blaze.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/Submatrix.h>
#include <blaze/math/Row.h>
#include <blaze/math/Rows.h>
#include <blaze/math/Column.h>
#include <blaze/math/DynamicVector.h>
//typedef SimpleMatrix2D Matrix;
typedef blaze::DynamicMatrix<double,blaze::rowMajor> Matrix;
typedef blaze::DynamicMatrix<double,blaze::columnMajor> cMatrix;
typedef blaze::DynamicVector<double,blaze::rowMajor> Vector;

#include <opencv2/opencv.hpp>

//#define opencvMatrixType CV_64FC1
#define opencvMatrixType CV_32FC1
class SimpleMatrix2D
{
public:
	SimpleMatrix2D();
	SimpleMatrix2D(int in_size1, int in_size2);
	void setSize(int in_size1, int in_size2);
	int columns();
	int rows();
	double getElement(int i, int j);
	void setElement(int i, int j, double val);
	~SimpleMatrix2D();
private:
	int dim;
	int* size;
	double* value;
};

//class BlazeMatrix2D{}

#endif /* MATRIX_H_ */
