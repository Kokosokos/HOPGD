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
	CreateModes():c_nmax(500)
	{
	}
	CreateModes(const inputData& inData,int nmax=500);

	void fitNew(const double& newParam1, Matrix& result) const;
	void fitNew2(const double& newParam1, Matrix& result) const;
	void fitNew(const double& newParam1, Matrix& result, int nModes) const;
	void fitNew3(const double& newParam1, Matrix& result, int nModes) const;
	 //aproximation rank, # of Modes
	//MVector Model;

	Model model;

private:
	const int c_nmax;
};

#endif /* CREATEMODES_H_ */
