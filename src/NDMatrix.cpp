#include "NDMatrix.h"


//===================================================================================================

//---------------------------------------------------------------------------------------------------

NDMatrix::NDMatrix()
{
}
//---------------------------------------------------------------------------------------------------
NDMatrix::NDMatrix(int dim, intVector sizes)
{
	m_current_position=0;
	setSize(dim, sizes);

}
NDMatrix::NDMatrix(const InputData& data)
{

	intVector sizes(3);

	sizes[0]=data.spaceDegreOfFreedom;
	sizes[1]=data.timeDegreOfFreedom;
	sizes[2]=data.param1DegreOfFreedom;
	setSize(3, sizes);
	for(int p=0;p<m_paramSize;++p)
	{
		m_values[p]=data.A[p];
	}
}
//---------------------------------------------------------------------------------------------------
void NDMatrix::setSize(int dim, intVector sizes)
{
	m_current_position=0;
	m_current_param_position=0;
	m_dim=dim;
	m_sizes=sizes;
	m_totSize=1;
	m_cumul_sizes=intVector(dim,1);
	m_current_index=intVector(dim,0);



	for(int i=0; i< m_dim;++i)
	{
		m_totSize*=m_sizes[i];
		m_cumul_sizes[i]=m_totSize;
	}
	m_matrix_size=(m_sizes[0]*m_sizes[1]);
	m_paramSize=m_totSize/m_matrix_size;
	for(int i=0;i<m_paramSize;++i)
	{
			m_values.push_back(Matrix(m_sizes[0],m_sizes[1],0));
	}

}
//---------------------------------------------------------------------------------------------------
intVector NDMatrix::sizes() const
{
	return m_sizes;
}
//---------------------------------------------------------------------------------------------------
int NDMatrix::dimensionality() const
{
	return m_dim;
}
//---------------------------------------------------------------------------------------------------
bool  NDMatrix::next() //Next parameter
{
	if(m_current_param_position < m_paramSize)
	{
		++m_current_param_position;
		m_current_position+=m_matrix_size;
		return true;
	}

	return false;
}
//---------------------------------------------------------------------------------------------------

void	NDMatrix::incrementIndex(int d)
{

	m_current_index[d]++;
	if(m_current_index[d]==m_sizes[d])
	{
		m_current_index[d]=0;
		incrementIndex(d+1);
	}

}
//---------------------------------------------------------------------------------------------------
bool  NDMatrix::ifEndPosition()
{
	return 	m_current_position >=m_totSize;
}
//---------------------------------------------------------------------------------------------------
bool  NDMatrix::ifNotEndPosition()
{
	return 	m_current_position < m_totSize;
}
//---------------------------------------------------------------------------------------------------

bool NDMatrix::getNext(Matrix& val)
{
	if(m_current_param_position < m_paramSize)
	{
		val=m_values[m_current_param_position];

		return next();
	}
	return false;
}
//---------------------------------------------------------------------------------------------------

bool NDMatrix::setNext(const Matrix& val)
{
	if(m_current_param_position < m_paramSize)
	{
		m_values[m_current_param_position]=val;
//		m_current_position++;
		next();
		return true;
	}
	return false;
}

//---------------------------------------------------------------------------------------------------
bool NDMatrix::goToIdx(intVector idx)
{
	return getPosFromIdx(idx, m_current_position);
}
//---------------------------------------------------------------------------------------------------

bool NDMatrix::goToPos(int pos)
{
#ifndef NOCHECKS
	if( !checkPos(pos) )
	{
		std::cerr<<"ERROR: NDMatrix::goToPos: pos= "<<pos<<" and maxPos="<<m_totSize<<std::endl;
		return false;
	}
#endif
	m_current_position=pos;
	m_current_param_position=pos;
	return true;
}

//---------------------------------------------------------------------------------------------------

void NDMatrix::resetCurrentPosition()
{
	m_current_index=intVector(m_dim,0);
	if(!goToPos(0))
		std::cerr<<"ERROR: NDMatrix::resetCurrentPosition:"<<std::endl;

}
//---------------------------------------------------------------------------------------------------
//double NDMatrix::getCurrentVal()
//{
//#ifndef NOCHECKS
//	if(m_current_position >= m_totSize)
//	{
//		std::cerr<<"Warning: NDMatrix::getCurrentVal: reading out of array bounds="<<std::endl;
//			return m_values[0];
//	}
//#endif
//
//	return m_values[m_current_position];
//}
//---------------------------------------------------------------------------------------------------

int NDMatrix::getCurrentPosition()
{
	return m_current_position;
}
//---------------------------------------------------------------------------------------------------

void NDMatrix::getCurrentIdx(intVector& idx)
{

#ifndef NOCHECKS
	if((m_current_position >= m_totSize))
	{
		std::cerr<<"ERROR: NDMatrix::getCurrentIdx:\n";
		idx=intVector(m_dim,0);
	}
#endif

	getIdxFromPos(idx, m_current_position);

}
//---------------------------------------------------------------------------------------------------
//double NDMatrix::getElement(const intVector idx)
//{
//	int pos=0;
//	getPosFromIdx(idx, pos);
//	return m_values[pos];
//}
////---------------------------------------------------------------------------------------------------
//void NDMatrix::setElement(intVector idx, double val)
//{
//	int pos=0;
//	getPosFromIdx(idx, pos);
//	m_values[pos]=val;
//}
//---------------------------------------------------------------------------------------------------
//double NDMatrix::getElement(int pos)
//{
//#ifndef NOCHECKS
//	if( !checkPos(pos) )
//	{
//		std::cerr<<"ERROR: NDMatrix::getElement: pos= "<<pos<<" and maxPos="<<m_totSize<<std::endl;
//	}
//#endif
//	return m_values[pos];
//}
//---------------------------------------------------------------------------------------------------
//void NDMatrix::setElement(int pos, double val)
//{
//#ifndef NOCHECKS
//	if( !checkPos(pos) )
//	{
//		std::cerr<<"ERROR: NDMatrix::setElement: pos= "<<pos<<" and maxPos="<<m_totSize<<std::endl;
//	}
//
//#endif
//	m_values[pos]=val;
//}
//---------------------------------------------------------------------------------------------------
bool NDMatrix::getPosFromIdx(const intVector& idx, int& pos)
{
#ifndef NOCHECKS
	if(idx[0] >= m_sizes[0])
	{
		std::cerr<<"ERROR: NDMatrix::getPosFromIdx: index is out of bounds!!!\n";
		return false;
	}
#endif

	pos=idx[0];
	int temp_size=1;
	for(int i=1; i< m_dim;++i)
	{
		temp_size*=m_sizes[i-1];
		if(idx[i] >= m_sizes[i])
		{
			std::cerr<<"ERROR: NDMatrix::getPosFromIdx: index is out of bounds!!!\n";
			return false;
		}
		pos+=idx[i]*temp_size;
	}
	return true;
}
//---------------------------------------------------------------------------------------------------
bool NDMatrix::getIdxFromPos(intVector& idx, int pos)
{
#ifndef NOCHECKS
	if( !checkPos(pos) )
	{
		std::cerr<<"ERROR: NDMatrix::getIdxFromPos: position is out of bounds. pos="<<pos<<std::endl;
		return false;
	}
#endif

	for(int i=m_dim-1; i>0;--i)
	{
		idx[i]=int(pos/m_cumul_sizes[i-1]);
		pos-=idx[i]*m_cumul_sizes[i-1];
	}
	idx[0]=pos;

	return true;


}
bool NDMatrix::checkPos(int pos)
{
	return pos<m_totSize;
}
//---------------------------------------------------------------------------------------------------
bool NDMatrix::ifCompatible(const NDMatrix& m2)
{
	if(m2.m_dim != m_dim)
	{
			std::cerr<<"ERROR: NDMatrix::ifCompatible: matrices are not compatible: different dimensions:"<<m_dim<<"!="<<m2.m_dim<<std::endl;
			return false;
	}

	if (sum(m_sizes-m2.m_sizes) !=0)
	{
			std::cerr<<"ERROR: NDMatrix::ifCompatible: matrices are not compatible: different sizes."<<std::endl;
			return false;
	}
	return true;
}
//---------------------------------------------------------------------------------------------------
//"Distance" between two matrices
double NDMatrix::dist( NDMatrix m2)
{
#ifndef NOCHECKS
	if(!ifCompatible(m2))
	{
		std::cerr<<"ERROR: NDMatrix:: cannot calculate distance between incompatible matrices."<<std::endl;
		return -1.0;
	}
#endif
	Vector norms=Vector(m_paramSize,0);
	for(int i=0; i<m_paramSize;++i )
	{
		norms[i]=blaze::norm(m_values[i])/blaze::norm(m2.m_values[i]);
	}
	return double(max(norms));
}
//---------------------------------------------------------------------------------------------------
NDMatrix& NDMatrix::operator - (const NDMatrix& m2)
{
#ifndef NOCHECKS
	if(!ifCompatible(m2))
	{
		std::cerr<<"ERROR: NDMatrix:: cannot subtract incompatible matrices."<<std::endl;
		return *this;
	}
#endif
	for(int i=0;i<m_paramSize;++i)
		m_values[i] = m_values[i] - m2.m_values[i];
	return *this;
}
//---------------------------------------------------------------------------------------------------

//NDMatrix NDMatrix::operator + (const NDMatrix& m2)
NDMatrix& NDMatrix::operator + (const NDMatrix& m2)
{
#ifndef NOCHECKS
	if(!ifCompatible(m2))
	{
		std::cerr<<"ERROR: NDMatrix:: cannot subtract incompatible matrices."<<std::endl;
		return *this;
	}
#endif

	for(int i=0;i<m_paramSize;++i)
		m_values[i] = m_values[i] + m2.m_values[i];
	return *this;
}
//---------------------------------------------------------------------------------------------------
NDMatrix::~NDMatrix()
{

}
//---------------------------------------------------------------------------------------------------
