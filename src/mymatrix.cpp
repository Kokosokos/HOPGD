#include "mymatrix.h"

void blaze2opencv(const Matrix & source, cv::Mat &dest)
{
	int xsize=source.columns();
	int ysize=source.rows();
	dest.create(ysize,xsize,opencvMatrixType);
	for (int i=0;i<ysize;++i)
		for (int j=0;j<xsize;++j)
		{
			if (opencvMatrixType==CV_32FC1)
				dest.at<float>(i,j,0)=float(source(i,j));
			else if (opencvMatrixType==CV_64FC1)
				dest.at<double>(i,j,0)=source(i,j);
		}

}

void opencv2blaze(const cv::Mat source, Matrix &dest)
{
	int xsize=source.cols;
	int ysize=source.rows;
	//std::cout<<"xsize="<<xsize<<" ysize="<<ysize<<"\n";
	dest=Matrix(ysize,xsize);
	for (int i=0;i<ysize;++i)
		for (int j=0;j<xsize;++j)
		{
			if (opencvMatrixType==CV_32FC1)
				dest(i,j)=double(source.at<float>(i,j,0));
			else if (opencvMatrixType==CV_64FC1)
				dest(i,j)=source.at<double>(i,j,0);
		}
}
//===================================================================================================
//---------------------------------------------------------------------------------------------------
myMatrix::myMatrix()
{
}
//---------------------------------------------------------------------------------------------------
myMatrix::myMatrix(int dim, intVector sizes)
{
	m_current_position=0;
	setSize(dim, sizes);
}
//---------------------------------------------------------------------------------------------------
void myMatrix::setSize(int dim, intVector sizes)
{
	m_dim=dim;
	m_sizes=sizes;
	m_totSize=1;
	for(int i=0; i< m_dim;++i)
		m_totSize*=m_sizes[i];
	m_values=Vector(m_totSize,0);

}
//---------------------------------------------------------------------------------------------------
intVector myMatrix::sizes() const
{
	return m_sizes;
}
//---------------------------------------------------------------------------------------------------
int myMatrix::dimensionality() const
{
	return m_dim;
}
//---------------------------------------------------------------------------------------------------
bool  myMatrix::next()
{
	if(m_current_position < m_totSize)
	{
		m_current_position++;
		return true;
	}
	return false;
}
//---------------------------------------------------------------------------------------------------
bool  myMatrix::ifEndPosition()
{
	return 	m_current_position >=m_totSize;
}
//---------------------------------------------------------------------------------------------------

bool myMatrix::getNext(double& val)
{
	if(m_current_position < m_totSize)
	{
		val=m_values[m_current_position];
		m_current_position++;
		return true;
	}
	return false;
}
//---------------------------------------------------------------------------------------------------

bool myMatrix::setNext(const double& val)
{
	if(m_current_position < m_totSize)
	{
		m_values[m_current_position]=val;
		m_current_position++;
		return true;
	}
	return false;
}
//---------------------------------------------------------------------------------------------------
bool myMatrix::goToIdx(intVector idx)
{
	return getPosFromIdx(idx, m_current_position);
}
//---------------------------------------------------------------------------------------------------

bool myMatrix::goToPos(int pos)
{
	if( !checkPos(pos) )
	{
		std::cerr<<"ERROR: myMatrix::goToPos: pos= "<<pos<<" and maxPos="<<m_totSize<<std::endl;
		return false;
	}
	m_current_position=pos;
	return true;
}

//---------------------------------------------------------------------------------------------------

void myMatrix::resetCurrentPosition()
{
	goToPos(0);
}
//---------------------------------------------------------------------------------------------------
double myMatrix::getCurrentVal()
{
	if(m_current_position < m_totSize)
	{
		return m_values[m_current_position];
	}
	std::cerr<<"Warning: myMatrix::getCurrentVal: reading out of array bounds="<<std::endl;
	return m_values[0];
}
//---------------------------------------------------------------------------------------------------

int myMatrix::getCurrentPosition()
{
	return m_current_position;
}
//---------------------------------------------------------------------------------------------------

void myMatrix::getCurrentIdx(intVector& idx)
{
	if (!getIdxFromPos(idx, m_current_position))
		std::cerr<<"ERROR: myMatrix::getCurrentIdx:\n";

}

//---------------------------------------------------------------------------------------------------
double myMatrix::getElement(const intVector idx)
{
	int pos=0;
	getPosFromIdx(idx, pos);
	return m_values[pos];
}
//---------------------------------------------------------------------------------------------------
void myMatrix::setElement(intVector idx, double val)
{
	int pos=0;
	getPosFromIdx(idx, pos);
	m_values[pos]=val;
}
//---------------------------------------------------------------------------------------------------
double myMatrix::getElement(int pos)
{

	if( !checkPos(pos) )
	{
		std::cerr<<"ERROR: myMatrix::getElement: pos= "<<pos<<" and maxPos="<<m_totSize<<std::endl;
	}
	return m_values[pos];
}
//---------------------------------------------------------------------------------------------------
void myMatrix::setElement(int pos, double val)
{
	if( !checkPos(pos) )
	{
		std::cerr<<"ERROR: myMatrix::setElement: pos= "<<pos<<" and maxPos="<<m_totSize<<std::endl;
	}
	m_values[pos]=val;
}
//---------------------------------------------------------------------------------------------------
bool myMatrix::getPosFromIdx(const intVector& idx, int& pos)
{
	if(idx[0] >= m_sizes[0])
	{
		std::cerr<<"ERROR: myMatrix::getPosFromIdx: index is out of bounds!!!\n";
		return false;
	}

	pos=idx[0];
	int temp_size=1;
	for(int i=1; i< m_dim;++i)
	{
		temp_size*=m_sizes[i-1];
		if(idx[i] >= m_sizes[i])
		{
			std::cerr<<"ERROR: myMatrix::getPosFromIdx: index is out of bounds!!!\n";
			return false;
		}
		pos+=idx[i]*temp_size;
	}
	return true;
}
//---------------------------------------------------------------------------------------------------
bool myMatrix::getIdxFromPos(intVector& idx, int pos)
{
	if( !checkPos(pos) )
	{
		std::cerr<<"ERROR: myMatrix::getIdxFromPos: position is out of bounds. pos="<<pos<<std::endl;
		return false;
	}
		int temp_size=m_totSize;
		for(int i=m_dim-1; i>0;--i)
		{
			temp_size/=m_sizes[i];
			idx[i]=int(pos/temp_size);
			pos-=idx[i]*temp_size;
		}
		idx[0]=pos;

	return true;

}
bool myMatrix::checkPos(int pos)
{
	return pos<m_totSize;
}
//---------------------------------------------------------------------------------------------------
bool myMatrix::ifCompatible(const myMatrix& m2)
{
	if(m2.m_dim != m_dim)
	{
			std::cerr<<"ERROR: myMatrix::ifCompatible: matrices are not compatible: different dimensions:"<<m_dim<<"!="<<m2.m_dim<<std::endl;
			return false;
	}

	if (sum(m_sizes-m2.m_sizes) !=0)
	{
			std::cerr<<"ERROR: myMatrix::ifCompatible: matrices are not compatible: different sizes."<<std::endl;
			return false;
	}
	return true;
}
//---------------------------------------------------------------------------------------------------
//"Distance" between two matrices
double myMatrix::dist( myMatrix m2)
{
	if(!ifCompatible(m2))
	{
		std::cerr<<"ERROR: myMatrix:: cannot calculate distance between incompatible matrices."<<std::endl;
		return -1.0;
	}

	Vector temp=m_values;
	temp = temp*temp;
	int size2d=(m_sizes[m_dim-1]*m_sizes[m_dim-2]);
	Vector norms=Vector(m_totSize/size2d,0);
	for(int i=0; i<norms.size()-1;++i )
	{
		auto sv1 = subvector(temp, i*size2d,size2d);
		auto sv2 = subvector(m2.m_values*m2.m_values, i*size2d,size2d);
		norms[i]=sqrt(blaze::sum(sv1))/sqrt(blaze::sum(sv2));
	}
	return double(max(norms));
}
//---------------------------------------------------------------------------------------------------

myMatrix myMatrix::operator - (const myMatrix& m2)
{
	if(!ifCompatible(m2))
	{
		std::cerr<<"ERROR: myMatrix:: cannot subtract incompatible matrices."<<std::endl;
		return *this;
	}

	myMatrix temp=m2;
	temp.m_values = m_values - m2.m_values;

	return temp;
}
//---------------------------------------------------------------------------------------------------

myMatrix myMatrix::operator + (const myMatrix& m2)
{
	if(!ifCompatible(m2))
	{
		std::cerr<<"ERROR: myMatrix:: cannot subtract incompatible matrices."<<std::endl;
		return *this;
	}

	myMatrix temp=m2;
	temp.m_values = m_values + m2.m_values;

	return temp;
}
//---------------------------------------------------------------------------------------------------
myMatrix::~myMatrix()
{

}
//---------------------------------------------------------------------------------------------------
