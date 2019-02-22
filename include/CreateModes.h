/*
 * CreateModes.h
 *
 *  Created on: 3 Dec 2018
 *      Author: ikriuchevs
 */

#ifndef CREATEMODES_H_
#define CREATEMODES_H_

#include "fileManager.h"
#include "NDMatrix.h"
//typedef blaze::DynamicVector<Matrix,blaze::rowMajor> MVector;

/**
 * @brief  A model creation class. Creates, stores and uses the model to fit new parameters.
 */
class CreateModes
{
public:
	CreateModes():m_nmax(500),m_ifCudaInit(false)
	{
	}

	/**
	 * @brief N Dimensional HOPGD method. Creates model and saves it in CreateModes.nmodel.
	 * @param[in] input input data structure
	 * @see NinputData3
	 */
	CreateModes(const NinputData3& input);


	/**
	 * @brief N Dimensional HOPGD method. Recreates model and saves it in CreateModes.nmodel.
	 * @param[in] input input data structure
	 * @see NinputData3
	 */
	void create(const NinputData3& input);

	/**
	 * @brief N Dimensional method that uses the model to produce the approximated result for a new parameter.
	 * @param[in] newParam1 The vector of the parameters values.
	 * @param[out] result The resulting approximation.
	 */
	bool fitNewND(const Vector& newParam1, Matrix& result, int nModes) const;

	/**
	 * @brief N Dimensional cuda initialization routine.  Copies nmodel to cvnmodel.
	 * @param[in] nModes Specifies the maximum number of modes used in approximation. Default =-1: use the same number of modes as in model.
	 */
	void cudaInitND(int nModes=-1);

	/**
	 * @brief N Dimensional method that uses the model to produce the approximated result for a new parameter.
	 * Uses openCV::Cuda for matrix-matrix multiplication.
	 * @param[in] newParam1 The vector of the parameters values.
	 * @param[out] result The resulting approximation stored in opencv::Mat class.
	 */
	bool fitNewNDCuda(const Vector& newParam1, cv::Mat &result, int nModes) const;

	/**
	 * The ND model structure that uses opencv::Mat.
	 */
	NModel nmodel;
	/**
	 * The ND model structure that uses opencv::Mat.
	 */
	cvNModel cvnmodel;

private:

	void findApprox( NDMatrix& R, Vector* F);
	/**
	 * @brief Finds the basis function F[dimId] assuming that all the rest F[*] are known.  \f$ R= \prod_{d=1}^{N}  F_d \f$
	 * @param dimId The dimension index of the basis vector that is unknown.
	 * @param R Residuals matrix that is being decomposed.
	 * @param F the vector of basis vectors.
	 */
	void findF(const int& dimId, NDMatrix& R, Vector* F);

	/**
	 * @brief N Dimensional HOPGD routine.
	 * @param M N Dimensional to be decomposed.
	 * @param ec Precision.
	 */
	void HOPGD( const NDMatrix& M, double ec);

	/**
	 * The maximum number of modes used in model.
	 */
	int m_nmax;
	bool m_ifCudaInit;
};

#endif /* CREATEMODES_H_ */
