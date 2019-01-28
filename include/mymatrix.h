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
//typedef myMatrix Matrix;
typedef blaze::DynamicMatrix<double,blaze::rowMajor> Matrix;
typedef blaze::DynamicMatrix<double,blaze::columnMajor> cMatrix;
typedef blaze::DynamicVector<double,blaze::rowMajor> Vector;
typedef blaze::DynamicVector<int,blaze::rowMajor> intVector;

#include <opencv2/opencv.hpp>
#include "opencv2/core/cuda.hpp"
#include "opencv2/cudaarithm.hpp"
#include "opencv2/cudaimgproc.hpp"

//#define opencvMatrixType CV_64FC1
#define opencvMatrixType CV_32FC1

void blaze2opencv(const Matrix & source, cv::Mat &dest);
void opencv2blaze(const cv::Mat source, Matrix &dest);


//Add copy operator!!
struct inputData
{
	//number of result files:
	int   timeDegreOfFreedom;
	//number of rows in each result file
	int  spaceDegreOfFreedom;
	//number of different Temperature(or another param) points = #of folders
	int param1DegreOfFreedom;
	float error;
	Vector param1;

	Matrix* A;
	inputData(){};
	void init(){A=new Matrix[param1DegreOfFreedom];};
	~inputData(){ delete[] A;};
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
	Matrix* F; //size dim

	Vector params; //must be Vector *params
	NModel(int d){dim=d; F=new Matrix[dim];nbrModes=0;};
	~NModel(){delete[] F;};
};
//===================================================================================================
//---------------------------------------------------------------------------------------------------

class myMatrix
{
public:
	myMatrix();
	myMatrix(int dim, intVector sizes);
	void setSize(int dim, intVector sizes);
	intVector sizes() const;

	int dimensionality() const;

	bool next();
	bool ifEndPosition();
	bool getNext(double& val);
	bool setNext(const double& val);

	bool goToIdx(intVector idx);
	bool goToPos(int pos);
	void resetCurrentPosition();
	double getCurrentVal();
	int getCurrentPosition();

	void getCurrentIdx(intVector& idx);


	double getElement(int pos);
	double getElement(const intVector idx);
	void setElement(int pos, double val);
	void setElement(intVector idx, double val);

	double dist( myMatrix m2);

	myMatrix operator - (const myMatrix& m2);

	myMatrix operator + (const myMatrix& m2);

	~myMatrix();
private:
	bool ifCompatible(const myMatrix& m2);
	bool getPosFromIdx(const intVector& idx, int& pos);
	bool getIdxFromPos(intVector& idx, int pos);
	bool checkPos(int pos);


	int m_dim;
	intVector m_sizes;
	int m_totSize; // m_sizes[0]*m_sizes[1]....
	Vector m_values;

	int m_current_position;
};

//class BlazeMatrix2D{}

#endif /* MATRIXVECTOR_H_ */
