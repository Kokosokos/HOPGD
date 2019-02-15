//============================================================================
// Name        : melted.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include "mymatrix.h"
#include "fileManager.h"
#include <boost/filesystem.hpp>
#include "CreateModes.h"
#include <omp.h>



using namespace boost::filesystem;

using namespace std;


int main(int argc, char * argv[])
{
	int nthreads = 8;				//number of threads for openmp blaze-lib parallelization
	FileManager file_manager;		//input files reader class
	NinputData3 nData;				//input data container

	string snapshotsFolder = 		/*snapshots folder*/
			"/home/ikriuchevs/workspace/melted/tests/Ndim//Jobs_bash_test/";

	string modelFolder =			/*folder to store the model files*/
			"/home/ikriuchevs/workspace/melted/tests/Ndim//Jobs_bash_test/181model/";


	std::cout<<"Reading Snapshots: ";
	double start = omp_get_wtime();
	file_manager.readFolder(snapshotsFolder, nData );
	std::cout<< omp_get_wtime() - start <<" sec"<<"\n";

	nData.error 	= 0.001;		//Model precision
	nData.nModesMax = 20;			//Maximum number of modes

//	std::cout<<"Creating modes: ";
//	start = omp_get_wtime();
//	CreateModes modesCreator3(nData);
//	std::cout<< omp_get_wtime() - start <<" sec"<<"\n";
//
//	std::cout<<"Saving the model: ";
//	start = omp_get_wtime();
//	file_manager.saveModel(modelFolder, modesCreator3.nmodel);
//	std::cout<< omp_get_wtime() - start <<" sec"<<"\n";
//


//	CreateModes modesCreator3;
//	file_manager.loadModel(modelFolder, modesCreator3.nmodel);
//	cout<<modesCreator3.nmodel.F[0].rows()<<" "<<modesCreator3.nmodel.F[1].rows()<<endl;

//
//	string folder="/home/ikriuchevs/workspace/melted/tests/Ndim/Validation/";
//	Matrix M(modesCreator3.nmodel.F[0].rows(),modesCreator3.nmodel.F[1].rows());
//	Matrix Mfit(modesCreator3.nmodel.F[0].rows(),modesCreator3.nmodel.F[1].rows());
//
//	double start = omp_get_wtime();
//	double fin;
//
//	FILE * outFile;
//	string outFileName=folder+"/nmax"+std::to_string( modesCreator3.nmodel.nbrModes)+".omp"+std::to_string(nthreads)+".dat";
//	outFile = fopen(outFileName.c_str(), "w");
//	fprintf (outFile, "#nbr of modes = %d\n",modesCreator3.nmodel.nbrModes);
//	fprintf (outFile, "#Folder Error Time CudaError CudaTime\n");
//	blaze::setNumThreads( nthreads );
//	//CUDA
//	modesCreator3.cudaInitND(20);
//	cv::Mat res;
//	res.create(modesCreator3.nmodel.F[0].rows(),modesCreator3.nmodel.F[1].rows(),opencvMatrixType);
//
//	for(int i=1;i<=8;++i)
//	{
//		string filename=folder+std::to_string(i)+"/"+c_odb_filename;//<--odb_filename defined in
//		cout<<filename<<endl;
//		file_manager.readODB2(filename, M);
//
//		string paramname=folder+std::to_string(i)+"/params.dat";
////		cout<<paramname<<endl;
//		std::ifstream inFile;
//		inFile.open(paramname, std::ios_base::in);
//		double temp=0;
//		Vector param(dim-2,0);
//		for(int d=0;d<dim-2;++d)
//		{
//			inFile>>temp;
//			if(d>=0)
//			{
//				param[dim-1-2-(d)]=temp;
//			}
//		}
//
//		inFile.close();
//		start = omp_get_wtime();
//		modesCreator3.fitNewND(param,Mfit);
//		fin=omp_get_wtime()-start;
//
//		double err=fabs(100.0*(1.0-blaze::norm(M)/blaze::norm(Mfit)));
//
//		cout<<blaze::trans(param)<<" : OMP: error: "<<err<<"%; time = "<<fin<<" sec;";
//		fprintf (outFile, "%d %f %lf ", i, err, fin );
//
//		start = omp_get_wtime();
//		modesCreator3.fitNewNDCuda(param,res);
//		fin = omp_get_wtime()-start;
//
//		opencv2blaze(res,Mfit);
//
//		err=fabs(100.0*(1.0-blaze::norm(M)/blaze::norm(Mfit)));
//		std::cout<<" CUDA: error="<<err<<"%; time = "<<fin  << " sec"<<std::endl;
//		fprintf (outFile, "%f %lf\n", err, fin );
//
//	}
//	fclose (outFile);

	return(0);
}




