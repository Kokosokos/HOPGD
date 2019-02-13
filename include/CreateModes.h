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

/**
 * A model creation class. Creates, stores and uses the model to fit new parameters.
 */
class CreateModes
{
public:
	CreateModes():c_nmax(500),m_ifCudaInit(false)
	{
	}
	/**
	 * @brief 3 Dimensional HOPGD method. Creates model and saves it in CreateModes.model.
	 * @param[in] inData input data structure
	 * @param[in] nmax the maximum number of modes (500 by default)
	 * @see inputData
	 */
	CreateModes(const inputData& inData,int nmax=500);

	/**
	 * @brief N Dimensional HOPGD method. Creates model and saves it in CreateModes.nmodel.
	 * @param[in] input input data structure
	 * @see NinputData3
	 */
	CreateModes(const NinputData3& input);


	/**
	 * @brief 3 Dimensional method that uses the model to produce the approximated result for a new parameter. Uses triple nested loop, very slow.
	 * @param[in] newParam1 Parameters value.
	 * @param[out] result The resulting approximation.
	 */
	void fitNew_Loops(const double& newParam1, Matrix& result) const;
	/**
	 * @brief 3 Dimensional method that uses the model to produce the approximated result for a new parameter. Uses vectors/matrix operations instead of triple loop.
	 * @param[in] newParam1 Parameters value.
	 * @param[out] result The resulting approximation.
	 */
	void fitNew(const double& newParam1, Matrix& result) const;
	/**
	* @brief 3 Dimensional method that uses the model to produce the approximated result for a new parameter.
	* Uses vectors/matrix operations instead of triple loop. Uses limited number of modes.
	* @param[in] newParam1 Parameters value.
	* @param[out] result The resulting approximation.
	* @param[in] nModes Specifies the maximum number of modes used in approximation.
	*/
	void fitNew(const double& newParam1, Matrix& result, int nModes) const;

	/**
	 * @brief N Dimensional method that uses the model to produce the approximated result for a new parameter.
	 * @param[in] newParam1 The vector of the parameters values.
	 * @param[out] result The resulting approximation.
	 */
	void fitNewND(const Vector& newParam1, Matrix& result) const;



	/**
	 * @brief 3 Dimensional cuda initialization routine. Copies model to cvmodel.
	 * @param[in] nModes Specifies the maximum number of modes used in approximation. Default =-1: use the same number of modes as in model.
	 */
	void cudaInit(int nModes=-1);

	/**
	 * @brief N Dimensional cuda initialization routine.  Copies nmodel to cvnmodel.
	 * @param[in] nModes Specifies the maximum number of modes used in approximation. Default =-1: use the same number of modes as in model.
	 */
	void cudaInitND(int nModes=-1);

	/**
	 * @brief 3 Dimensional method that uses the model to produce the approximated result for a new parameter.
	 * Uses openCV::Cuda for matrix-matrix multiplication.
	 * @param[in] newParam1 Parameters value.
	 * @param[out] result The resulting approximation stored in opencv::Mat class.
	 */
	void fitNewCuda(const double& newParam1, cv::Mat &result) const;
	/**
	 * @brief N Dimensional method that uses the model to produce the approximated result for a new parameter.
	 * Uses openCV::Cuda for matrix-matrix multiplication.
	 * @param[in] newParam1 The vector of the parameters values.
	 * @param[out] result The resulting approximation stored in opencv::Mat class.
	 */
	void fitNewNDCuda(const Vector& newParam1, cv::Mat &result) const;




	/**
	 * The 3D model structure.
	 */
	Model model;
	/**
	 * The 3D model structure that uses opencv::Mat.
	 */
	cvModel cvmodel;
	/**
	 * The ND model structure that uses opencv::Mat.
	 */
	NModel nmodel;
	/**
	 * The ND model structure that uses opencv::Mat.
	 */
	cvNModel cvnmodel;

private:

	void findApprox( myMatrix3& R, Vector* F);
	/**
	 * @brief Finds the basis function F[dimId] assuming that all the rest F[*] are known.  \f$ R= \prod_{d=1}^{N}  F_d \f$
	 * @param dimId The dimension index of the basis vector that is unknown.
	 * @param R Residuals matrix that is being decomposed.
	 * @param F the vector of basis vectors.
	 */
	void findF(const int& dimId, myMatrix3& R, Vector* F);

	/**
	 * @brief N Dimensional HOPGD routine.
	 * @param M N Dimensional to be decomposed.
	 * @param ec Precision.
	 */
	void HOPGD( const myMatrix3& M, double ec);

	/**
	 * The maximum number of modes used in model.
	 */
	const int c_nmax;
	bool m_ifCudaInit;
};

#endif /* CREATEMODES_H_ */
