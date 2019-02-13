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

class NDMatrix
{
public:
	NDMatrix();
	NDMatrix(const inputData& data);
	~NDMatrix();
	NDMatrix(int dim, intVector sizes);

	void		setSize(int dim, intVector sizes);
	intVector	sizes() const;

	int		dimensionality() const;

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

	double	dist( NDMatrix m2);

	NDMatrix& operator - (const NDMatrix& m2);
	NDMatrix& operator + (const NDMatrix& m2);

//private:
	bool	ifCompatible(const NDMatrix& m2);
	bool	getPosFromIdx(const intVector& idx, int& pos);
	bool	getIdxFromPos(intVector& idx, int pos);
	bool	checkPos(int pos);



	int			m_dim;
	intVector	m_sizes;
	intVector	m_cumul_sizes;

	int			m_totSize; // m_sizes[0]*m_sizes[1]....
	int			m_matrix_size; // m_sizes[0]*m_sizes[1]....
	int			m_paramSize;

	std::vector<Matrix> m_values;

	int			m_current_position;
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
	NinputData3(inputData data):A(data),error(data.error){params.push_back(data.param1);};
};

#endif /* MATRIXVECTOR2_H_ */
