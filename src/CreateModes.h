/*
 * CreateModes.h
 *
 *  Created on: 3 Dec 2018
 *      Author: ikriuchevs
 */

#ifndef CREATEMODES_H_
#define CREATEMODES_H_

#include "fileManager.h"
#include "mymatrix3.h"
//typedef blaze::DynamicVector<Matrix,blaze::rowMajor> MVector;


class CreateModes
{
public:
	CreateModes():c_nmax(500),m_ifCudaInit(false)
	{
	}
	//3D method
	CreateModes(const inputData& inData,int nmax=500);
	//ND method
	CreateModes(const NinputData& input);
	CreateModes(const NinputData3& input);

	//using triple nested loop
	void fitNew_Loops(const double& newParam1, Matrix& result) const;

	//using blaze matrix multiplication
	void fitNewND(const Vector& newParam1, Matrix& result) const;
	void fitNewND2(const Vector& newParam1, Matrix& result) const;

	void fitNew(const double& newParam1, Matrix& result) const;
	void fitNew(const double& newParam1, Matrix& result, int nModes) const;

	//using opencv::cuda matrix multiplication
	void cudaInit(int nModes=-1);
	void fitNewCuda(const double& newParam1, cv::Mat &result) const;



	Model model;
	cvModel cvmodel;
	NModel nmodel;

private:
	void findF(const int& dimId, myMatrix& R, Vector* F);
	void findApprox( myMatrix& R, Vector* F);

	void findApprox( myMatrix3& R, Vector* F);
	void findF(const int& dimId, myMatrix3& R, Vector* F);

	void HOPGD( const myMatrix& M, double ec);
	void HOPGD( const myMatrix3& M, double ec);
	const int c_nmax;
	bool m_ifCudaInit;
};

#endif /* CREATEMODES_H_ */
