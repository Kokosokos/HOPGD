/*
 * CreateModes.h
 *
 *  Created on: 3 Dec 2018
 *      Author: ikriuchevs
 */

#ifndef CREATEMODES_H_
#define CREATEMODES_H_

#include "fileManager.h"
//typedef blaze::DynamicVector<Matrix,blaze::rowMajor> MVector;


class CreateModes
{
public:
	CreateModes():c_nmax(500),m_ifCudaInit(false)
	{
	}
	CreateModes(const inputData& inData,int nmax=500);

	//using triple nested loop
	void fitNew_Loops(const double& newParam1, Matrix& result) const;

	//using blaze matrix multiplication
	void fitNew(const double& newParam1, Matrix& result) const;
	void fitNew(const double& newParam1, Matrix& result, int nModes) const;

	//using opencv::cuda matrix multiplication
	void cudaInit(int nModes=-1);
	void fitNewCuda(const double& newParam1, cv::Mat &result) const;



	Model model;
	cvModel cvmodel;

private:
	const int c_nmax;
	bool m_ifCudaInit;
};

#endif /* CREATEMODES_H_ */
