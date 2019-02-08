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

