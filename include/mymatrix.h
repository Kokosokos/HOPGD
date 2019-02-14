/*
 * matrixVector.h
 *
 *  Created on: 6 Nov 2018
 *      Author: ikriuchevs
 */

#ifndef MATRIXVECTOR_H_
#define MATRIXVECTOR_H_

#include <blaze/Blaze.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/Submatrix.h>
#include <blaze/math/Row.h>
#include <blaze/math/Rows.h>
#include <blaze/math/Column.h>
#include <blaze/math/DynamicVector.h>
#include <blaze/math/Subvector.h>

typedef blaze::DynamicMatrix<double,blaze::rowMajor> Matrix;
typedef blaze::DynamicMatrix<double,blaze::columnMajor> cMatrix;
typedef blaze::DynamicVector<double,blaze::rowMajor> Vector;
typedef blaze::DynamicVector<int,blaze::rowMajor> intVector;

#include <opencv2/opencv.hpp>
#include "opencv2/core/cuda.hpp"
#include "opencv2/cudaarithm.hpp"
#include "opencv2/cudaimgproc.hpp"

#include <vector>
//#define opencvMatrixType CV_64FC1
#define opencvMatrixType CV_32FC1
#define NOindexArray
#define NOCHECKS
void blaze2opencv(const Matrix & source, cv::Mat &dest);
void opencv2blaze(const cv::Mat source, Matrix &dest);


//Add copy operator?!
struct InputData
{
	//number of result files:
	int		timeDegreOfFreedom;
	//number of rows in each result file
	int		spaceDegreOfFreedom;
	//number of different Temperature(or another param) points = #of folders
	int		param1DegreOfFreedom;
	float	error;
	Vector	param1;

	Matrix* A;
	InputData(){};
	void init(){A=new Matrix[param1DegreOfFreedom];};
	~InputData(){ delete[] A;};
};



struct Model
{
	int nbrModes;
	Matrix F1;
	Matrix F2;
	Matrix F3;

	Vector param1;
};

struct cvModel
{
	int nbrModes;
	cv::Mat F1;
	cv::Mat F2;
	cv::Mat F3;

	Vector param1;
};

struct NModel
{
	int nbrModes;
	int dim;
//	Matrix* F; //size dim
//	NModel(int d){dim=d; F=new Matrix[dim];nbrModes=0;};
//	~NModel(){delete[] F;};

	std::vector<Matrix> F;
	Vector param1;
//	Vector* params; //must be Vector *params
	std::vector<Vector> params;
	NModel(){};
	NModel(int d){dim=d; F=std::vector<Matrix>(dim);nbrModes=0;};
//	NModel operator = (NModel nm){ NModel temp(dim); temp.nbrModes=nm.nbrModes;for (int i=0;i<dim;++i)temp.F[i]=nm.F[i];return temp;};

};
struct cvNModel
{
	int nbrModes;
	int dim;
	std::vector<cv::Mat> F;
	std::vector<Vector> params;

	cvNModel(){};
	cvNModel(int d){dim=d; F=std::vector<cv::Mat>(dim);nbrModes=0;};
};



#endif /* MATRIXVECTOR_H_ */
