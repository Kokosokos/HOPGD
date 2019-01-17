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


//vector p must be sorted
int searchUpperNearest(const Vector& p,double v)
{
	for(uint i =0;i<p.size();++i)
	{
		if(v<p[i])
			return i;
	}
	return p.size();

}
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
					R[k]=inData.A[k]-fna[k];
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
		model.cF2=trans(model.F2);
		model.cF3=trans(model.F3);
	}
	//-----------------------------------------------------------------------


	delete[] f;
	delete[] fna;
}

void CreateModes::fitNew2(const double& newParam1, Matrix& result) const
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
		#pragma omp parallel for
		        for (uint j=0; j<model.F2.rows();++j)
		            for(uint mode=0;mode<model.nbrModes;++mode)
		            	result(i,j) = result(i,j)+  model.F1(i,mode)*model.F2(j,mode)*Vinter[mode];

	}
	else
	{
		std::cout<<"Enrichment is needed\n";
	}


}

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

void CreateModes::fitNew3(const double& newParam1, Matrix& result, int nModes) const
{
	Vector Vinter(model.nbrModes,0);
	//param1.find(newParam1);

	int idx=searchUpperNearest(model.param1,newParam1);
//	std::cout<<idx<<" "<<newParam1<<"< "<<param1[idx]<<"\n";
	if( idx > 0 && idx < model.param1.size())
	{

		Vector low=column(model.cF3,idx-1);
		Vector up=column(model.cF3,idx);
		Vector Vinter=low+(up-low)*(1.0/(model.param1[idx]-model.param1[idx-1])*(newParam1-model.param1[idx-1]));



		result.resize(model.F1.rows(),model.F2.rows(),false);
		result=0;
		Matrix result2=result;
		cMatrix F2temp=model.cF2;
		F2temp.resize(nModes,model.cF2.columns());
		for(uint mode=0;mode<nModes;++mode)
		{
			double temp=Vinter[mode];
			row(F2temp,mode)=temp*row(F2temp,mode);
		}

		result= submatrix(model.F1,0,0,model.F1.rows(),nModes)*F2temp;
	}
	else
	{
		std::cout<<"Enrichment is needed\n";
	}


}
