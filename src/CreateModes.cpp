/*
 * CreateModes.cpp
 *
 *  Created on: 3 Dec 2018
 *      Author: ikriuchevs
 */

#include "CreateModes.h"
#include <iostream>
#include <omp.h>
//#include "test_opencv.cpp"
const double	c_ebr  = 0.000001;
const int		c_vmax = 2400;
const double	c_eBc  = 0.000001;
//---------------------------------------------------------------------------------------------------
//vector p must be sorted
int searchUpperNearest(const Vector& p,double v)
{
	for(uint i =0;i<p.size();++i)
	{
		if(v<=p[i])
			return i;
	}
	return p.size();

}
//---------------------------------------------------------------------------------------------------

void CreateModes::HOPGD( const NDMatrix& M, double ec)
{
	intVector sizes = M.sizes();
	int dim = M.dimensionality();

	NModel model(dim);

	//~f^n - n-th mode of approximation. M = sum over n of f^n
	NDMatrix f(dim, sizes);

	//Residuals matrix
	NDMatrix R = M;

	//"Unknown functions of n-th mode". f^n = F[0] * F[1] * ....
	Vector *F;
	F = new Vector[dim];
	std::cout<<std::endl;
	//Loop over modes
	for (int n = 0; n < m_nmax; ++n)
	{
		// The product of all F norms except one
		// normBF = ||F[0]||*...||F[m-1]||*||F[m+1]||*...*||F[N]||
		double normBF = 1;
		// Squared  F Norm of the previous iteration: ||F[m-1]||
		double prevNorm = 1;

		// Convergence check. If the "distance" between R and M is small -> exit
		if(R.dist(M) < ec )
		{
			break;
		}

		// The norms of the F[m]
		// Current step norm
		Vector bt1  = Vector(dim,0);
		// Previous step norm
		Vector bta1 = Vector(dim,0);
		// 1st step norm
		Vector br1  = Vector(dim,1);

		//Set initial values
		for(int d = 0; d < dim; ++d)
		{
			F[d] = Vector(sizes[d], 1);
			normBF *= blaze::sqrNorm(F[d]);
		}
		//Converge n-th mode
		for (int v = 0; v < c_vmax; ++v)
		{
			std::cout<< "\rn: " << n << "/"<<m_nmax<< "; v:" << v << "/"<<c_vmax<<"                "<<std::flush;
			if(v == c_vmax-1)
			{
				std::cout<<"WARNING: inner loop v==vmax;"<<std::endl;
			}

			bool e1 = true;

			//Loop over dimensions
			for(int m = 0; m < dim; ++m)
			{
				// Squared F Norm of the current iteration: ||F_{m-1}||^2
				double currNorm = blaze::sqrNorm(F[m]);
				normBF /= currNorm;
				normBF *= prevNorm;

				//Finds F[m] assuming that all the rest F[k] (k!=m) are known.
				findF(m,R,F);
				F[m] = F[m] / normBF;

				// Recalculate the norm of updated F[m]
				bt1[m]   = blaze::norm(F[m]);
				prevNorm = bt1[m]*bt1[m];

				if(v==0)
				{
					br1[m] = bt1[m];
				}

				//Convergence criterion of n-th mode
				e1 = e1 && ((bt1[m] - bta1[m]) / br1[m] < c_eBc);
			}

			if(e1)
			{
				//Calculates f from all F[m]. f = F[0] * F[1] * ....
				findApprox(f, F);

				//Update residuals
				R = R - f;

				// Save all F[m] of the n-th mode in the model container
				for(int m = 0; m < dim; ++m)
				{
					model.F[m].resize(F[m].size(), model.nbrModes + 1);
					column(model.F[m], model.nbrModes) = F[m];
				}

				model.nbrModes++;


				break;
			}
			else
			{
				bta1 = bt1;
			}

		}


	}
	std::cout<<"\rDone                                       "<<std::endl;
	delete[] F;
	nmodel = model;

}

//---------------------------------------------------------------------------------------------------


void CreateModes::findF(const int& dimId, NDMatrix& R, Vector* F)
{
	R.resetCurrentPosition();
	int dim     = R.dimensionality();
	F[dimId]    = Vector(F[dimId].size(), 1);
	Vector newF = Vector(F[dimId].size(), 0);
	intVector idx(dim, 0);

	double fMult;

	if(dimId == 0)
	{
		for(int k = 0; k < R.m_paramSize; ++k)
		{
			R.getCurrentIdx(idx);

			fMult = 1;
			for(int d = 2; d < dim; ++d)
			{
				fMult *= F[d][idx[d]];
			}
			newF += R.m_values[k] * F[1] * fMult;
			R.next();
		}
	}
	else if(dimId == 1)
	{
		for(int k = 0; k < R.m_paramSize; ++k)
		{
			R.getCurrentIdx(idx);
			fMult = 1;
			for(int d = 2; d < dim; ++d)
			{
				fMult *= F[d][idx[d]];
			}
			newF += trans(R.m_values[k]) * F[0] * fMult;
			R.next();
		}

	}
	else
	{
		for(int k = 0;k < R.m_paramSize; ++k)
		{
			R.getCurrentIdx(idx);
			fMult = 1;
			for(int d = 2; d < dim; ++d)
				fMult *= F[d][idx[d]];
			newF[idx[dimId]] += (trans(F[0]) *R.m_values[k] *F[1]) *fMult;
			R.next();
		}
	}

	F[dimId]=newF;
}
//---------------------------------------------------------------------------------------------------
void CreateModes::findApprox( NDMatrix& R, Vector* F)
{

	int dim = R.dimensionality();
	intVector idx(dim, 0);
	R.resetCurrentPosition();
	Matrix M(R.sizes()[0], R.sizes()[1]);
	for(int k = 0; k < R.m_paramSize; ++k)
	{
		R.getCurrentIdx(idx);
		double fMult = 1;
		for(int d = 2; d < dim; d++)
			fMult *= F[d][idx[d]];
		M = F[0] * fMult * trans(F[1]);
		R.setNext(M);
	}

}
//---------------------------------------------------------------------------------------------------
CreateModes::CreateModes(const NinputData3& input)
{
	create(input);
}

void CreateModes::create(const NinputData3& input)
{
	m_nmax=input.nModesMax;
	HOPGD(input.A,input.error);
	nmodel.params = input.params;
}

//---------------------------------------------------------------------------------------------------
bool CreateModes::fitNewND(const Vector& newParam1, Matrix& result, int nModes) const
{

		Vector Vinter(nmodel.nbrModes,1);
		for(int d = 2; d < nmodel.dim; ++d)
		{
			int idx = searchUpperNearest(nmodel.params[d - 2], newParam1[d - 2]);
			if( (idx > 0 && idx < nmodel.params[d - 2].size()) || (idx ==0 && newParam1==nmodel.params[d-2][idx]))
			{
				int idxlow=idx==0?idx:idx-1;
				Vector low=column(trans(nmodel.F[d]),idxlow);
				Vector up=column(trans(nmodel.F[d]),idx);

				Vector Vinter2;
				if(idx ==0 && newParam1==nmodel.params[d-2][idx])
				{
					Vinter2=low;
				}
				else
				{


					Vinter2=low+(up-low)*(1.0/(nmodel.params[d-2][idx]-nmodel.params[d-2][idx-1])*(newParam1[d-2]-nmodel.params[d-2][idx-1]));
				}
				Vinter=Vinter*Vinter2;
			}
			else
			{
				std::cout<<"Enrichment is needed\n";
				return false;
			}
		}
		result.resize(nmodel.F[0].rows(),nmodel.F[1].rows(),false);
		result=0;
		Matrix result2=result;
		Matrix F1temp=submatrix(nmodel.F[1],0,0,nmodel.F[1].rows(),nModes);//nmodel.F[1];

		for(uint mode=0;mode<nModes;++mode)
		{
			double temp=Vinter[mode];
			column(F1temp,mode)=temp*column(F1temp,mode);
		}
		result=submatrix(nmodel.F[0],0,0,nmodel.F[0].rows(),nModes)*trans(F1temp);
		return true;
}

//---------------------------------------------------------------------------------------------------

void CreateModes::cudaInitND(int nModes)
{
	cvnmodel=cvNModel(nmodel.dim);
//	cvnmodel.dim=nmodel.dim;

	int rows=nmodel.F[0].rows();
	int cols=nmodel.F[1].rows();

	if (nModes<0 || nModes > nmodel.nbrModes)
	{
		nModes=nmodel.nbrModes;
		for(int d=0;d<cvnmodel.dim;++d)
			blaze2opencv(nmodel.F[d],cvnmodel.F[d]);

	}
	else
	{
		for(int d=0;d<cvnmodel.dim;++d)
			blaze2opencv(submatrix(nmodel.F[d],0,0,nmodel.F[d].rows(),nModes),cvnmodel.F[d]);
	}
	cvnmodel.nbrModes=nModes;
	int CudaDevice = cv::cuda::getDevice();
	cv::cuda::setDevice(CudaDevice);
	cv::Mat m = cv::Mat::ones(rows, cols, opencvMatrixType);
	cv::cuda::GpuMat d_src;
	cv::cuda::GpuMat d_src2;
	cv::cuda::GpuMat d_src3(cols,cols,opencvMatrixType);
	d_src.upload(m);
	d_src2.upload(m.t());


	cv::cuda::gemm(d_src2, d_src, 1.0, cv::cuda::GpuMat(), 0.0, d_src3);
	m_ifCudaInit=true;
	d_src.download(m);
}


//---------------------------------------------------------------------------------------------------

bool CreateModes::fitNewNDCuda(const Vector& newParam1, cv::Mat &result, int nModes) const
{

		Vector Vinter(nmodel.nbrModes,1);
		for(int d=2;d<nmodel.dim;++d)
		{
			int idx=searchUpperNearest(nmodel.params[d-2],newParam1[d-2]);
			if( (idx > 0 && idx < nmodel.params[d-2].size()) || (idx ==0 && newParam1==nmodel.params[d-2][idx]))
			{
				int idxlow=idx==0?idx:idx-1;
				Vector low=column(trans(nmodel.F[d]),idxlow);
				Vector up=column(trans(nmodel.F[d]),idx);

				Vector Vinter2;
				if(idx ==0 && newParam1==nmodel.params[d-2][idx])
				{
					Vinter2=low;
				}
				else
				{


					Vinter2=low+(up-low)*(1.0/(nmodel.params[d-2][idx]-nmodel.params[d-2][idx-1])*(newParam1[d-2]-nmodel.params[d-2][idx-1]));
				}
				Vinter=Vinter*Vinter2;
			}
			else
			{
				std::cout<<"Enrichment is needed\n";
				return false;
			}
		}

		cv::Mat F1temp=cvnmodel.F[1].clone();

		for(uint mode=0;mode<cvnmodel.nbrModes;++mode)
		{
			double temp=Vinter[mode];
			F1temp.col(mode)=temp*F1temp.col(mode);
		}

		cv::cuda::GpuMat dummy;
		cv::cuda::GpuMat gf1;
		gf1.upload(cvnmodel.F[0]);
		cv::cuda::GpuMat gf2;
		gf2.upload(F1temp);
		cv::cuda::transpose(gf2,gf2);
		cv::cuda::GpuMat gres;
		cv::cuda::gemm(gf1, gf2, 1.0, cv::cuda::GpuMat(), 0.0, gres);
		gres.download(result);
		return true;
}
