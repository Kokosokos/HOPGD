/*
 * matrixVector.h
 *
 *  Created on: 6 Nov 2018
 *      Author: ikriuchevs
 */

#ifndef MATRIXVECTOR3_H_
#define MATRIXVECTOR3_H_
#include "mymatrix.h"



//===================================================================================================
//---------------------------------------------------------------------------------------------------
/**
 * @brief N-dimensional matrix(tensor) container class.
 *
 * Considers first 2 dimensions as largest (space & time) and all the rest as additional parameters.
 * Stores N-dimensional matrix as a vector of 2D matrices, but can return an N-dimensional index.
 *
 */
class NDMatrix
{
public:
	NDMatrix();
	/**
	 * @brief Constructor compatible with "old" 3D InputData structure.
	 */
	NDMatrix(const InputData& data);
	~NDMatrix();
	/**
	 * @brief Default Constructor to be used.
	 * @param dim Dimensionality of the matrix.
	 * @param sizes Vector of sizes for each dimension.
	 */
	NDMatrix(int dim, intVector sizes);

	/**
	 * @brief This function sets the size of the matrix and all its elements to zero.
	 * Used to reinitialize already existing matrix.
	 * @param dim Dimensionality of the matrix.
	 * @param sizes Vector of sizes for each dimension.
	 */
	void		setSize(int dim, intVector sizes);
	/**
	 * @brief Returns the vector that contains the sizes of all dimensions.
	 * @return NDMatrix.m_sizes
	 */
	intVector	sizes() const;

	/**
	 * @brief Returns matrix dimensionality.
	 * @return NDMatrix.m_dim
	 */
	int		dimensionality() const;

	/**
	 * @brief Increments the NDMatrix.m_current_param_position.
	 */
	bool	next();

	void	incrementIndex(int d);
	bool	ifEndPosition();
	bool	ifNotEndPosition();
	bool	getNext(Matrix& val);
	bool	setNext(const Matrix& val);
	bool	append(const Vector& m, int paramId);



	void  resetCurrentPosition();
	double	getCurrentVal();
	int		getCurrentPosition();
	void	getCurrentIdx(intVector& idx);
	bool	goToIdx(intVector idx);
	bool	goToPos(int pos);


//	double	getElement(int pos);
//	double	getElement(const intVector idx);
//	void	setElement(int pos, double val);
//	void	setElement(intVector idx, double val);

	/**
	 * @brief Returns the "distance" between 2 matrices.
	 * @param m2 The matrix to which the "distance" is calculated.
	 * @return \f$ max\left(\frac{\left|\left| m_1(i,j,p_1,p_2,...) \right|\right|_{ij}}{\left|\left| m_2(i,j,p_1,p_2...) \right|\right|_{ij}}\right) \f$
	 */
	double	dist( NDMatrix m2);

	/**
		 * @brief ND Matrix subtraction operator.
		 */
	NDMatrix& operator - (const NDMatrix& m2);

	/**
		 * @brief ND Matrix summation operator.
		 */
	NDMatrix& operator + (const NDMatrix& m2);

//private:
	bool	ifCompatible(const NDMatrix& m2);
	bool	getPosFromIdx(const intVector& idx, int& pos);
	bool	getIdxFromPos(intVector& idx, int pos);
	bool	checkPos(int pos);



	/**
	 * @brief Dimensionality of the matrix.
	 */
	int			m_dim;
	/**
	 * @brief The vector containing the sizes of each dimension.
	 */
	intVector	m_sizes;
	/**
		 * @brief The vector containing the cumulative sizes of each dimension.
		 *
		 * m_cumul_sizes[0]=m_sizes[0],
		 * m_cumul_sizes[n]=m_cumul_sizes[n-1]*m_sizes[n].
		 * Used to speed-up the ND index calculation.
		 */
	intVector	m_cumul_sizes;

	/**
		 * @brief The total number of elements in the ND-matrix.
		 *
		 * m_totSize=m_sizes[0]*m_sizes[1]*...m_sizes[N-1]
		 */
	int			m_totSize;
	/**
		 * @brief The number of elements in the 2D-matrix.
		 *
		 * m_matrix_size=m_sizes[0]*m_sizes[1]
		 */
	int			m_matrix_size;
	/**
			 * @brief The total size of the parameters space.
			 *
			 * m_paramSize=m_sizes[2]*m_sizes[3]*...m_sizes[N-1]
			 */
	int			m_paramSize;

	/**
				 * @brief The actual values in the N-dimensional matrix (as a vector of 2D matrices).
				 */
	std::vector<Matrix> m_values;

	int			m_current_position;

	/**
		 * @brief ND Matrix current parameter position.
		 */
	int			m_current_param_position;
	intVector	m_current_index;
};

struct NinputData3
{

	float	error;
	int nModesMax;
	std::vector<Vector> params;

	NDMatrix A;
	NinputData3(){};
	NinputData3(InputData data):A(data),error(data.error){params.push_back(data.param1);};
};

#endif /* MATRIXVECTOR2_H_ */
