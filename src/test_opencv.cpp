/*
 * test_opencv.cpp
 *
 *  Created on: 19 Dec 2018
 *      Author: ikriuchevs
 */


#include <opencv2/opencv.hpp>
#include "opencv2/core/cuda.hpp"
#include "opencv2/cudaarithm.hpp"
#include "opencv2/cudaimgproc.hpp"
#include "mymatrix.h"
#include "CreateModes.h"

#include <omp.h>


//Dummy function
//Trying to exclude first (very slow) cuda malloc
void cudaInit(int rows,int cols)
{
	int CudaDevice = cv::cuda::getDevice();
	cv::cuda::setDevice(CudaDevice);
	cv::Mat m = cv::Mat::ones(rows, cols, opencvMatrixType);
	cv::cuda::GpuMat d_src;
	cv::cuda::GpuMat d_src2;
	cv::cuda::GpuMat d_src3(cols,cols,opencvMatrixType);
	d_src.upload(m);
	d_src2.upload(m.t());


	cv::cuda::gemm(d_src2, d_src, 1.0, cv::cuda::GpuMat(), 0.0, d_src3);

//	d_src.download(m);
}
void blaze2opencv(const Matrix & source, cv::Mat &dest)
{
	int xsize=source.columns();
	int ysize=source.rows();
	//std::cout<<"xsize="<<xsize<<" ysize="<<ysize<<"\n";
	dest.create(ysize,xsize,opencvMatrixType);
	for (int i=0;i<ysize;++i)
		for (int j=0;j<xsize;++j)
		{
//			dest.at<double>(i,j,0)=source(i,j);
			if (opencvMatrixType==CV_32FC1)
				dest.at<float>(i,j,0)=float(source(i,j));
			else if (opencvMatrixType==CV_64FC1)
				dest.at<double>(i,j,0)=source(i,j);
		}

	//	namedWindow( "Display window", WINDOW_AUTOSIZE );// Create a window for display.
	//	imshow( "Display window", dest );                   // Show our image inside it.
	//
	//	waitKey(0);
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

	//	namedWindow( "Display window", WINDOW_AUTOSIZE );// Create a window for display.
	//	imshow( "Display window", dest );                   // Show our image inside it.
	//
	//	waitKey(0);
}

int testMult()
{
	cv::Mat sourceMat = cv::Mat::ones(1080, 1920, CV_32FC3);

	//This is the color Matrix
	float matrix[3][3] = {
			{ 1.057311, -0.204043, 0.055648 }
			, { 0.041556, 1.875992, -0.969256 }
			, { -0.498535, -1.537150, 3.240479 }
	};

	cv::Mat colorMatrixMat = cv::Mat(3, 3, CV_32FC1, matrix).t();

	cv::Mat linearSourceMat = sourceMat.reshape(1, 1080 * 1920);
	cv::Mat multipliedMatrix = linearSourceMat * colorMatrixMat;

	try {
		cv::Mat dummy, gpuMultipliedMatrix;

		// Regular gemm
		cv::gemm(linearSourceMat, colorMatrixMat, 1.0, dummy, 0.0, gpuMultipliedMatrix);
		// CUDA gemm
		// cv::cuda::gemm(linearSourceMat, colorMatrixMat, 1.0, dummy, 0.0, gpuMultipliedMatrix);

		std::cout << (cv::countNonZero(multipliedMatrix != gpuMultipliedMatrix) == 0);
	} catch (cv::Exception& e) {
		std::cerr << e.what();
		return -1;
	}
}
int searchUpperNearest(const Vector& p,double v);
void fitNew(const double& newParam1,cv::Mat &result, const cv::Mat &F1, const cv::Mat &F2, const Matrix& F3, const Vector& param1)
{

	double start, fin;
	start = omp_get_wtime();
	int nbrModes=F1.cols;
	Vector Vinter(nbrModes,0);
	//
	int idx=searchUpperNearest(param1,newParam1);
	std::cout<< omp_get_wtime() - start<<" ";//1

	//		std::cout<<idx<<" "<<newParam1<<"< "<<param1[idx]<<"\n";
	if( idx > 0 && idx < param1.size())
	{
		start = omp_get_wtime();
		Vector low=column(trans(F3),idx-1);
		Vector up=column(trans(F3),idx);
		Vector Vinter=low+(up-low)*(1.0/(param1[idx]-param1[idx-1])*(newParam1-param1[idx-1]));
		std::cout<< omp_get_wtime() - start<<" ";//2
		start = omp_get_wtime();

		//result.create(F1.rows,F2.rows,opencvMatrixType);
		cv::Mat F2temp=F2.clone();
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
		gf1.upload(F1);
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

		//
	}
	else
	{
		std::cout<<"Enrichment is needed\n";
	}

}
