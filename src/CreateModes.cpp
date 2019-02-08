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
const double c_ebr=0.000001;

const int c_vmax=500;

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

void CreateModes::HOPGD( const myMatrix3& M, double ec)
{
	//	const int c_nmax=100;
	double start,fin;
	//	start = omp_get_wtime();
	const double eBc=0.000001;
	const int vmax=2400;
	intVector sizes=M.sizes();
	int dim=M.dimensionality();

	NModel model(dim);
	//	myMatrix3 fna(dim, sizes);
	myMatrix3 f(dim, sizes);
	myMatrix3 R=M;
	Vector *F;
	F=new Vector[dim];
	//	std::cout<< "Prep time: "<<omp_get_wtime()-start<<std::endl;
	double findFtime=0;
	double findFnorm=0;
	for (int n=0;n<c_nmax;++n)
	{
		double normBF=1;
		double prevNorm=1;

		double start,fin;
		std::cout<<"\n";
		std::cout<<n<<": v: ";
		std::cout.flush();
		//		start = omp_get_wtime();

		if(R.dist(M) < ec )
		{
			std::cout<<"!!!!CONGRATULATIONS!!!!!\n";
			break;
		}

		Vector bta1=Vector(dim,0);
		Vector br1=Vector(dim,1);
		Vector bt1=Vector(dim,0);

		for(int d =0;d<dim;d++)
		{
			F[d]=Vector(sizes[d],1); //Must come from M or model
			normBF*=blaze::sqrNorm(F[d]);

		}

		for (int v=0;v<vmax;++v)
		{
			std::cout<<v<<" ";
			std::cout.flush();
			bool e1=true;

			//			start = omp_get_wtime();
			for(int m=0;m<dim;++m)
			{

				double currNorm=blaze::sqrNorm(F[m]);
				normBF/=currNorm;
				normBF*=prevNorm;
				findF(m,R,F);
				F[m]=F[m]/normBF;

				bt1[m]=blaze::norm(F[m]);
				prevNorm=bt1[m]*bt1[m];
				//				std::cout<<(trans(F[m])*F[m])<<"  "<<bt1[m]<<std::endl;
				if(v==0)
				{
					br1[m]=bt1[m];
				}

				e1=e1&&((bt1[m]-bta1[m])/br1[m]<eBc);
			}
			//			fin = omp_get_wtime();
			//			findFtime=findFtime+fin-start;
			if(e1)
			{
				findApprox( f, F);

				R=R-f;
				//				fna=fna+f;


				//Store n-mode functions in NModel structure
				for(int m=0;m<dim;++m)
				{
					model.F[m].resize(F[m].size(), model.nbrModes+1);
					column(model.F[m],model.nbrModes)=F[m];
				}
				model.nbrModes++;


				break;
			}
			else
			{
				bta1=bt1;
			}

		}


	}
	delete[] F;
	//	std::cout<<"findFtime: "<<findFtime<<std::endl;
	//	std::cout<<"findFnorm: "<<findFnorm<<std::endl;
	nmodel=model;

}
//---------------------------------------------------------------------------------------------------
void CreateModes::HOPGD( const myMatrix& M, double ec)
{

	double start,fin;
	//	start = omp_get_wtime();
	const double eBc=0.000001;
	const int vmax=200;
	intVector sizes=M.sizes();
	int dim=M.dimensionality();
	std::cout<<dim<<"\n";
	NModel model(dim);
	myMatrix fna(dim, sizes);
	myMatrix f(dim, sizes);
	myMatrix R=M;
	Vector *F;
	F=new Vector[dim];
	//	std::cout<< "Prep time: "<<omp_get_wtime()-start<<std::endl;
	double findFtime=0;
	double findFnorm=0;

	for (int n=0;n<c_nmax;++n)
	{
		double normBF=1;
		double prevNorm=1;

		double start,fin;
		std::cout<<n<<"\n";
		//		start = omp_get_wtime();

		if(R.dist(M) < ec )
		{
			//			std::cout<<"!!!!CONGRATULATIONS!!!!!\n";
			break;
		}
		//		fin = omp_get_wtime();
		//		findFnorm=findFnorm+fin-start;
		//Initialial values
		Vector bta1=Vector(dim,0);
		Vector br1=Vector(dim,1);
		Vector bt1=Vector(dim,0);

		for(int d =0;d<dim;d++)
		{
			F[d]=Vector(sizes[d],1); //Must come from M or model
			normBF*=blaze::sqrNorm(F[d]);

		}

		for (int v=0;v<vmax;++v)
		{
			bool e1=true;

			//			start = omp_get_wtime();
			for(int m=0;m<dim;++m)
			{

				double currNorm=blaze::sqrNorm(F[m]);
				normBF/=currNorm;
				normBF*=prevNorm;
				findF(m,R,F);
				F[m]=F[m]/normBF;

				bt1[m]=blaze::norm(F[m]);
				prevNorm=bt1[m]*bt1[m];
				//				std::cout<<(trans(F[m])*F[m])<<"  "<<bt1[m]<<std::endl;
				if(v==0)
				{
					br1[m]=bt1[m];
				}

				e1=e1&&((bt1[m]-bta1[m])/br1[m]<eBc);
			}
			//			fin = omp_get_wtime();
			//			findFtime=findFtime+fin-start;
			if(e1)
			{
				findApprox( f, F);

				R=R-f;
				//				fna=fna+f;


				//Store n-mode functions in NModel structure
				for(int m=0;m<dim;++m)
				{
					model.F[m].resize(F[m].size(), model.nbrModes+1);
					column(model.F[m],model.nbrModes)=F[m];
				}
				model.nbrModes++;


				break;
			}
			else
			{
				bta1=bt1;
			}

		}


	}
	delete[] F;
	//	std::cout<<"findFtime: "<<findFtime<<std::endl;
	//	std::cout<<"findFnorm: "<<findFnorm<<std::endl;
	nmodel=model;

}

//---------------------------------------------------------------------------------------------------
#ifndef NOindexArray
void CreateModes::findF(const int& dimId, myMatrix& R, Vector* F)
{
	//	R.resetCurrentPosition();
	F[dimId]=Vector(F[dimId].size(),1);
	Vector newF=Vector(F[dimId].size(),0);
	int dim=R.dimensionality();
	intVector idx(dim,0);

	int imax=R.m_totSize;

	for(int i=0;i <imax;++i)
	{
		double fMult=1;

		double r=R.getElement(i);
		idx=R.m_index[i];
		for(int d=0;d<dim;d++)
			fMult*=F[d][idx[d]];
		newF[idx[dimId]]+=r*fMult;
	}
	F[dimId]=newF;



}
#else
void CreateModes::findF(const int& dimId, myMatrix& R, Vector* F)
{
	R.resetCurrentPosition();
	F[dimId]=Vector(F[dimId].size(),1);
	Vector newF=Vector(F[dimId].size(),0);
	int dim=R.dimensionality();
	intVector idx(dim,0);
	double fMult;
	while(!R.ifEndPosition())
		//	while(R.ifNotEndPosition())
	{
		fMult=1;

		R.getCurrentIdx(idx);

		//		for(int d=0;d<dim;d++)
		//			fMult*=F[d][idx[d]];

		for(int d=0;d<dimId;++d)
			fMult*=F[d][idx[d]];
		for(int d=dimId+1;d<dim;++d)
			fMult*=F[d][idx[d]];
		newF[idx[dimId]]+=R.getCurrentVal()*fMult;
		R.next();
	}
	F[dimId]=newF;
}


#endif

//---------------------------------------------------------------------------------------------------

void CreateModes::findApprox( myMatrix& R, Vector* F)
{

	int dim=R.dimensionality();
	intVector idx(dim,0);
	R.resetCurrentPosition();
	while(!R.ifEndPosition())
	{
		R.getCurrentIdx(idx);
		double fMult=1;
		for(int d=0;d<dim;d++)
			fMult*=F[d][idx[d]];
		R.setNext(fMult);
	}

}
//---------------------------------------------------------------------------------------------------


void CreateModes::findF(const int& dimId, myMatrix3& R, Vector* F)
{
	R.resetCurrentPosition();
	F[dimId]=Vector(F[dimId].size(),1);
	Vector newF=Vector(F[dimId].size(),0);
	int dim=R.dimensionality();
	intVector idx(dim,0);
	double fMult;

	if(dimId ==0) //0 or 1 (space/time, or 2 max size dimensions)
	{
		for(int k=0;k<R.m_paramSize;++k)
		{
			R.getCurrentIdx(idx);
			fMult=1;
			for(int d=2;d<dim;++d)
			{
				fMult*=F[d][idx[d]];
			}
			newF+=R.m_values[k]*F[1]*fMult;
			R.next();
		}
		//		newF[idx[dimId]]+=(trans(F[0])*R.m_values[k]*F[1])*fMult;
	}
	else if(dimId==1)
	{
		for(int k=0;k<R.m_paramSize;++k)
		{
			R.getCurrentIdx(idx);
			fMult=1;
			for(int d=2;d<dim;++d)
			{
				fMult*=F[d][idx[d]];
			}
			newF+=trans(R.m_values[k])*F[0]*fMult;
			R.next();
		}

	}
	else
	{
		for(int k=0;k<R.m_paramSize;++k)
		{
			R.getCurrentIdx(idx);
			fMult=1;
			for(int d=2;d<dim;++d)
				fMult*=F[d][idx[d]];
			newF[idx[dimId]]+=(trans(F[0])*R.m_values[k]*F[1])*fMult;
			R.next();
		}
	}

	double normBF=1;
	//	for(int d =0;d<dim;d++)
	//	{
	//		if(d!=dimId)
	//		normBF*=blaze::sqrNorm(F[d]);
	//
	//	}
	F[dimId]=newF/normBF;
}
//---------------------------------------------------------------------------------------------------
void CreateModes::findApprox( myMatrix3& R, Vector* F)
{

	int dim=R.dimensionality();
	intVector idx(dim,0);
	R.resetCurrentPosition();
	Matrix M(R.sizes()[0],R.sizes()[1]);
	for(int k=0;k<R.m_paramSize;++k)
	{

		R.getCurrentIdx(idx);
		double fMult=1;
		for(int d=2;d<dim;d++)
			fMult*=F[d][idx[d]];
		M=F[0]*fMult*trans(F[1]);
		R.setNext(M);
	}

}
////---------------------------------------------------------------------------------------------------

//void CreateModes::findF(int dimId, myMatrix& R, Vector* F)
//{
//	F[dimId]=Vector(F[dimId].size(),1);
//	Vector temp=R.m_values;
//	Vector newF=Vector(F[dimId].size(),0);
//	int dim=R.dimensionality();
//
//	intVector cum_sizes(dim+1,0);
//	cum_sizes[0]=1;
//	for(int d=0;d<dim;++d)
//	{
//		cum_sizes[d+1]=cum_sizes[d]*R.m_sizes[d];
//
//	}
//	for(int d=0;d<dim;d++)
//	{
//		for(int j=0; j<R.m_totSize/cum_sizes[d+1];++j)
//		{
//			for (int f =0;f<F[d].size();++f)
//			{
//				auto sv2 = subvector( temp, j*cum_sizes[d+1], cum_sizes[d] );
//				sv2=sv2*F[d][f];
//			}
//		}
//
//	}
//	for(int j=0; j<R.m_totSize/cum_sizes[dimId+1];++j)
//	{
//		for (int f =0;f<F[dimId].size();++f)
//		{
//			auto sv2 = subvector( temp, j*cum_sizes[dimId+1], cum_sizes[dimId] );
//			newF[f]+=blaze::sum(sv2);
//		}
//	}
//	double ssss=blaze::sum(temp);
//
//	double normBF=1;
//	for(int d=0;d<dim;d++)
//	{
//		if (d!=dimId)
//			normBF*=(trans(F[d])*F[d]);
//	}
//	F[dimId]=newF/normBF;
//
//}
//---------------------------------------------------------------------------------------------------
CreateModes::CreateModes(const NinputData& input):c_nmax(input.nModesMax)
{
	HOPGD( input.A,input.error);
	nmodel.params=input.params;
	if(nmodel.dim==3)
	{
		model.nbrModes=nmodel.nbrModes;
		model.F1.resize(nmodel.F[0].rows(),nmodel.F[0].columns());
		model.F1=nmodel.F[0];
		model.F2=nmodel.F[1];
		model.F3=nmodel.F[2];
		model.param1=input.params[0];
	}
}
//---------------------------------------------------------------------------------------------------
CreateModes::CreateModes(const NinputData3& input):c_nmax(input.nModesMax)
{
	HOPGD( input.A,input.error);
	nmodel.params=input.params;
	if(nmodel.dim==3)
	{
		model.nbrModes=nmodel.nbrModes;
		model.F1.resize(nmodel.F[0].rows(),nmodel.F[0].columns());
		model.F1=nmodel.F[0];
		model.F2=nmodel.F[1];
		model.F3=nmodel.F[2];
		model.param1=input.params[0];
	}
}
//---------------------------------------------------------------------------------------------------
CreateModes::CreateModes(const inputData& inData,int nmax):c_nmax(nmax)
{
	model.nbrModes=0;
	int K=inData.param1DegreOfFreedom;
	model.param1=inData.param1;

	Vector F1(inData.spaceDegreOfFreedom,  1.0f);
	Vector F2(inData.timeDegreOfFreedom,   1.0f);
	Vector F3(inData.param1DegreOfFreedom, 1.0f);


	double bta1=0;
	double bta2=0;
	double bta3=0;
	Matrix* R;
	R=new Matrix[inData.param1DegreOfFreedom];
	//~f^m - mth mode of approximation
	Matrix* f;
	f=new Matrix[inData.param1DegreOfFreedom];
	//~f^n - nth rank approximation of X = sum over m of f^m
	Matrix* fna;
	fna=new Matrix[inData.param1DegreOfFreedom];
	Vector norm_X(K); //norm of inData.A

	//std::cout<<"Norms: ";
	for (int k =0;k < K;++k)
	{
		R[k]=Matrix(inData.spaceDegreOfFreedom,inData.timeDegreOfFreedom);
		f[k]=Matrix(inData.spaceDegreOfFreedom,inData.timeDegreOfFreedom);
		fna[k]=Matrix(inData.spaceDegreOfFreedom,inData.timeDegreOfFreedom);
		fna[k]=0;// by default?
		R[k]=inData.A[k];
		f[k]=0;
		norm_X[k]=blaze::norm(inData.A[k]);
		//std::cout<<inData.param1[k]<<" "<<norm_X[k]<<"\n";
	}

	//-----------------------------------------------------------------------
	//for n=1:500
	for(int n =0;n<c_nmax;++n)
	{
		std::cout<<"mode: "<<n<<std::endl;
		Vector norm_R(K); //norm of residuals

		for(int k=0; k < K; ++k)
		{
			norm_R[k]=blaze::norm(R[k]);
		}
		double e=blaze::norm(norm_R/norm_X); //e(n)
		if (e < inData.error)
			break;

		Vector M=Vector(inData.spaceDegreOfFreedom,0);
		for(int k=0; k < K; ++k)
		{
			double D=F3[k];
			M=M+R[k]*F2*D;
		}
		F1=M/((trans(F2)*F2)*(trans(F3)*F3)) ;

		M=Vector(inData.timeDegreOfFreedom,0);
		for(int k=0; k < K; ++k)
		{
			double D=F3[k];
			M=M+trans(R[k])*F1*D;
		}
		F2=M/((trans(F1)*F1)*(trans(F3)*F3)) ;

		for(int k=0; k < K; ++k)
		{
			F3[k]=(trans(F1)*R[k]*F2)/((trans(F2)*F2)*(trans(F1)*F1));
		}

		double br1=blaze::norm(F1);
		double br2=blaze::norm(F2);
		double br3=blaze::norm(F3);

		for (int v=0; v<c_vmax; ++v)
		{
			Vector M=Vector(inData.spaceDegreOfFreedom,0);

			for(int k=0; k < K; ++k)
			{
				double D=F3[k];
				M=M+R[k]*F2*D;
			}
			F1=M/((trans(F2)*F2)*(trans(F3)*F3)) ;

			M=Vector(inData.timeDegreOfFreedom,0);
			for(int k=0; k < K; ++k)
			{
				double D=F3[k];
				M=M+trans(R[k])*F1*D;
			}
			F2=M/((trans(F1)*F1)*(trans(F3)*F3)) ;

			for(int k=0; k < K; ++k)
			{
				F3[k]=(trans(F1)*R[k]*F2)/((trans(F2)*F2)*(trans(F1)*F1));
			}
			//std::cout<<std::endl;
			//same
			double bt1=blaze::norm(F1);
			double bt2=blaze::norm(F2);
			double bt3=blaze::norm(F3);

			//		double e1=blaze::norm(bt1-bta1)/br1;
			//		double e2=blaze::norm(bt2-bta2)/br2;
			//		double e3=blaze::norm(bt3-bta3)/br3;
			double e1=fabs(bt1-bta1)/br1;
			double e2=fabs(bt2-bta2)/br2;
			double e3=fabs(bt3-bta3)/br3;


			//std::cout<<"v: "<<v<<"; err: "<<e1<<" "<<e2<<" "<<e3<<std::endl;
			if(e1 < c_ebr && e2 < c_ebr && e3 < c_ebr)
			{
				for(int k=0; k < K; ++k)
				{
					double D=F3[k];
					f[k]=F1*D*trans(F2);
					fna[k]=fna[k]+f[k];
					R[k]=inData.A[k]-fna[k]; // same R[k]-=f[k]
				}

				//std::cout<<"nbrModes "<<nbrModes<<"\n";
				model.F1.resize(F1.size(), model.nbrModes+1);
				model.F2.resize(F2.size(), model.nbrModes+1);
				model.F3.resize(F3.size(), model.nbrModes+1);
				//std::cout<<"model size: "<<model.F1.rows()<<" "<<model.F1.columns()<<std::endl;
				column(model.F1,model.nbrModes)=F1;
				column(model.F2,model.nbrModes)=F2;
				column(model.F3,model.nbrModes)=F3;
				model.nbrModes++;

				bta1=0;
				bta2=0;
				bta3=0;
				F1=Vector(inData.spaceDegreOfFreedom,  1.0f);
				F2=Vector(inData.timeDegreOfFreedom,   1.0f);
				F3=Vector(inData.param1DegreOfFreedom, 1.0f);

				break;
			}
			else
			{
				bta1=bt1;
				bta2=bt2;
				bta3=bt3;
			}
			if(v==c_vmax-1)
			{
				std::cout<<"WARNING: inner loop v==vmax; e1,e2,e3:"<<e1<<" "<<e2<<" "<<e3<<std::endl;
			}

		}

	}
	//-----------------------------------------------------------------------


	delete[] f;
	delete[] fna;
}

//---------------------------------------------------------------------------------------------------

void CreateModes::fitNew_Loops(const double& newParam1, Matrix& result) const
{
	Vector Vinter(model.nbrModes,0);
	//param1.find(newParam1);

	int idx=searchUpperNearest(model.param1,newParam1);
	//	std::cout<<idx<<" "<<newParam1<<"< "<<param1[idx]<<"\n";
	if( idx > 0 && idx < model.param1.size())
	{
		//		std::cout<<"model.F3 size"<<model.F3.rows()<<std::endl;
		Vector low=column(trans(model.F3),idx-1);
		Vector up=column(trans(model.F3),idx);
		Vector Vinter=low+(up-low)*(1.0/(model.param1[idx]-model.param1[idx-1])*(newParam1-model.param1[idx-1]));

		result.resize(model.F1.rows(),model.F2.rows(),false);
		result=0;

		for( uint i=0; i<model.F1.rows();++i)
			for (uint j=0; j<model.F2.rows();++j)
				for(uint mode=0;mode<model.nbrModes;++mode)
					result(i,j) = result(i,j)+  model.F1(i,mode)*model.F2(j,mode)*Vinter[mode];

	}
	else
	{
		std::cout<<"Enrichment is needed\n";
	}


}

//---------------------------------------------------------------------------------------------------
void CreateModes::fitNewND(const Vector& newParam1, Matrix& result) const
{
	if(nmodel.dim==3)
	{
		fitNew(newParam1[0], result);
	}
	else
	{
		Vector Vinter(nmodel.nbrModes,1);
		for(int d=2;d<nmodel.dim;++d)
		{
			int idx=searchUpperNearest(nmodel.params[d-2],newParam1[d-2]);
			if( idx >= 0 && idx < nmodel.params[d-2].size())
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
			}
		}
		result.resize(nmodel.F[0].rows(),nmodel.F[1].rows(),false);
		result=0;
		Matrix result2=result;
		Matrix F1temp=nmodel.F[1];

		for(uint mode=0;mode<nmodel.nbrModes;++mode)
		{
			double temp=Vinter[mode];
			column(F1temp,mode)=temp*column(F1temp,mode);
		}
		result=nmodel.F[0]*trans(F1temp);

	}
}

//---------------------------------------------------------------------------------------------------

void CreateModes::fitNew(const double& newParam1, Matrix& result) const
{
	Vector Vinter(model.nbrModes,0);
	//param1.find(newParam1);

	int idx=searchUpperNearest(model.param1,newParam1);
	//	std::cout<<idx<<" "<<newParam1<<"< "<<param1[idx]<<"\n";
	if( idx > 0 && idx < model.param1.size())
	{
		Vector low=column(trans(model.F3),idx-1);
		Vector up=column(trans(model.F3),idx);
		Vector Vinter=low+(up-low)*(1.0/(model.param1[idx]-model.param1[idx-1])*(newParam1-model.param1[idx-1]));



		result.resize(model.F1.rows(),model.F2.rows(),false);
		result=0;
		Matrix result2=result;
		Matrix F2temp=model.F2;

		for(uint mode=0;mode<model.nbrModes;++mode)
		{
			double temp=Vinter[mode];
			column(F2temp,mode)=temp*column(F2temp,mode);
		}
		std::cout<<"\n";
		result=model.F1*trans(F2temp);

	}
	else
	{
		std::cout<<"Enrichment is needed\n";
	}
}

//---------------------------------------------------------------------------------------------------

void CreateModes::fitNew(const double& newParam1, Matrix& result, int nModes) const
{
	double start, fin;
	start = omp_get_wtime();

	Vector Vinter(model.nbrModes,0);
	//param1.find(newParam1);

	int idx=searchUpperNearest(model.param1,newParam1);
	std::cout<< omp_get_wtime() - start<<" ";//1
	//	std::cout<<idx<<" "<<newParam1<<"< "<<param1[idx]<<"\n";
	if( idx > 0 && idx < model.param1.size())
	{
		start = omp_get_wtime();
		Vector low=column(trans(model.F3),idx-1);
		Vector up=column(trans(model.F3),idx);
		Vector Vinter=low+(up-low)*(1.0/(model.param1[idx]-model.param1[idx-1])*(newParam1-model.param1[idx-1]));
		std::cout<< omp_get_wtime() - start<<" ";//2
		start = omp_get_wtime();


		result.resize(model.F1.rows(),model.F2.rows(),false);
		result=0;

		Matrix result2=result;
		Matrix F2temp=model.F2;
		F2temp.resize(model.F2.rows(),nModes);
		std::cout<< omp_get_wtime() - start<<" ";//3
		start = omp_get_wtime();
		for(uint mode=0;mode<nModes;++mode)
		{
			double temp=Vinter[mode];
			column(F2temp,mode)=temp*column(F2temp,mode);
		}
		std::cout<< omp_get_wtime() - start<<" ";//4
		start = omp_get_wtime();
		result= submatrix(model.F1,0,0,model.F1.rows(),nModes)*trans(F2temp);
		std::cout<< omp_get_wtime() - start<<" ";//5
	}
	else
	{
		std::cout<<"Enrichment is needed\n";
	}


}

//Just no cout
//void CreateModes::fitNew(const double& newParam1, Matrix& result, int nModes) const
//{
//	Vector Vinter(model.nbrModes,0);
//	//param1.find(newParam1);
//
//	int idx=searchUpperNearest(model.param1,newParam1);
////	std::cout<<idx<<" "<<newParam1<<"< "<<param1[idx]<<"\n";
//	if( idx > 0 && idx < model.param1.size())
//	{
//
//		Vector low=column(trans(model.F3),idx-1);
//		Vector up=column(trans(model.F3),idx);
//		Vector Vinter=low+(up-low)*(1.0/(model.param1[idx]-model.param1[idx-1])*(newParam1-model.param1[idx-1]));
//
//
//
//		result.resize(model.F1.rows(),model.F2.rows(),false);
//		result=0;
//		Matrix result2=result;
//		Matrix F2temp=model.F2;
//		F2temp.resize(model.F2.rows(),nModes);
//		for(uint mode=0;mode<nModes;++mode)
//		{
//			double temp=Vinter[mode];
//			column(F2temp,mode)=temp*column(F2temp,mode);
//		}
//
//		result= submatrix(model.F1,0,0,model.F1.rows(),nModes)*trans(F2temp);
//	}
//	else
//	{
//		std::cout<<"Enrichment is needed\n";
//	}
//
//
//}
//

//---------------------------------------------------------------------------------------------------

void CreateModes::cudaInit(int nModes)
{
	int rows=model.F1.rows();
	int cols=model.F2.rows();

	if (nModes<0 || nModes > model.nbrModes)
	{
		nModes=model.nbrModes;
		blaze2opencv(model.F1,cvmodel.F1);
		blaze2opencv(model.F2,cvmodel.F2);
		blaze2opencv(model.F3,cvmodel.F3);
	}
	else
	{
		blaze2opencv(submatrix(model.F1,0,0,model.F1.rows(),nModes),cvmodel.F1);
		blaze2opencv(submatrix(model.F2,0,0,model.F2.rows(),nModes),cvmodel.F2);
		blaze2opencv(submatrix(model.F3,0,0,model.F3.rows(),nModes),cvmodel.F3);
	}

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

void CreateModes::fitNewCuda(const double& newParam1,cv::Mat &result) const
{

	double start, fin;
	start = omp_get_wtime();
	int nbrModes=cvmodel.F1.cols;
	Vector Vinter(nbrModes,0);
	//
	int idx=searchUpperNearest(model.param1,newParam1);
	std::cout<< omp_get_wtime() - start<<" ";//1

	//		std::cout<<idx<<" "<<newParam1<<"< "<<param1[idx]<<"\n";
	if(( idx >= 0 && idx < model.param1.size()) )
	{
		start = omp_get_wtime();
		int idxlow=idx==0?idx:idx-1;
		Vector low = column(trans(model.F3), idxlow);
		Vector up  = column(trans(model.F3), idx);

		Vector Vinter;
		if(idx ==0 && newParam1==model.param1[idx])
		{
			Vinter=low;
		}
		else
		{


			Vinter=low+(up-low)*(1.0/(model.param1[idx]-model.param1[idx-1])*(newParam1-model.param1[idx-1]));
		}

		std::cout<< omp_get_wtime() - start<<" ";//2
		start = omp_get_wtime();

		//result.create(F1.rows,F2.rows,opencvMatrixType);
		cv::Mat F2temp=cvmodel.F2.clone();
		std::cout<< omp_get_wtime() - start<<" ";//3
		start = omp_get_wtime();
		for(uint mode=0;mode<nbrModes;++mode)
		{
			double temp=Vinter[mode];
			F2temp.col(mode)=temp*F2temp.col(mode);
		}
		std::cout<< omp_get_wtime() - start<<" ";//4

		start = omp_get_wtime();
		cv::cuda::GpuMat dummy;
		cv::cuda::GpuMat gf1;
		gf1.upload(cvmodel.F1);
		cv::cuda::GpuMat gf2;
		gf2.upload(F2temp);
		cv::cuda::transpose(gf2,gf2);
		cv::cuda::GpuMat gres;
		//			cv::gemm(F1, F2temp.t(), 1.0, cv::Mat(), 0.0, result);
		std::cout<< omp_get_wtime() - start<<" ";//5

		start = omp_get_wtime();
		cv::cuda::gemm(gf1, gf2, 1.0, cv::cuda::GpuMat(), 0.0, gres);
		std::cout<< omp_get_wtime() - start<<" ";//6

		start = omp_get_wtime();
		gres.download(result);
		std::cout<< omp_get_wtime() - start<<" ";//7

	}
	else
	{
		std::cout<<"Enrichment is needed: "<<idx<<std::endl;
	}

}
//---------------------------------------------------------------------------------------------------

void CreateModes::fitNewNDCuda(const Vector& newParam1, cv::Mat &result) const
{
	if(nmodel.dim==3)
	{
		fitNewCuda(newParam1[0], result);
	}
	else
	{
		Vector Vinter(nmodel.nbrModes,1);
		for(int d=2;d<nmodel.dim;++d)
		{
			int idx=searchUpperNearest(nmodel.params[d-2],newParam1[d-2]);
			if( idx >= 0 && idx < nmodel.params[d-2].size())
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
			}
		}

		cv::Mat F1temp=cvnmodel.F[1].clone();

		for(uint mode=0;mode<nmodel.nbrModes;++mode)
		{
			double temp=Vinter[mode];
			F1temp.col(mode)=temp*F1temp.col(mode);
		}

//		double start = omp_get_wtime();
		cv::cuda::GpuMat dummy;
		cv::cuda::GpuMat gf1;
		gf1.upload(cvnmodel.F[0]);
		cv::cuda::GpuMat gf2;
		gf2.upload(F1temp);
		cv::cuda::transpose(gf2,gf2);
		cv::cuda::GpuMat gres;

//		std::cout<< omp_get_wtime() - start<<" ";//5

//		start = omp_get_wtime();
		cv::cuda::gemm(gf1, gf2, 1.0, cv::cuda::GpuMat(), 0.0, gres);
//		std::cout<< omp_get_wtime() - start<<" ";//6

//		start = omp_get_wtime();
		gres.download(result);
//		std::cout<< omp_get_wtime() - start<<" ";//7

	}
}
