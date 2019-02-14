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


#include <odb_API.h>

int main(int argc, char * argv[])
{
	odb_initializeAPI();

	int dim=6;
	int nthreads=8;
	FileManager file_manager;
	NinputData3 ndata;
	intVector sizes(dim-2,0);

	//lARGE test
//	std::cout<<"Small TEST 1-81: ";
//	sizes[0]=3; //Temperature
//	sizes[1]=3; //Specific heat
//	sizes[2]=3; //Density
//	sizes[3]=3; //Density
//
//	string folder="/home/ikriuchevs/workspace/melted/tests/Ndim/Jobs/";
//	string outfolder="/home/ikriuchevs/workspace/melted/tests/Ndim/Jobs/181model/";
//
//	double start = omp_get_wtime();
//	std::cout<<"Reading Files: ";
//
//	file_manager.readFolder(folder,ndata,dim,sizes);
//	std::cout<< omp_get_wtime() - start <<" sec"<<"\n";
//	ndata.error=0.001;
//	ndata.nModesMax=20;
//
//	std::cout<<"Creating modes: ";
//	CreateModes modesCreator3(ndata);
//	std::cout<< omp_get_wtime() - start <<" sec"<<"\n";
//
//	std::cout<<"Saving the model: \n";
//	file_manager.saveModel(outfolder, modesCreator3.nmodel);

	string outfolder="/home/ikriuchevs/workspace/melted/tests/Ndim/Jobs/181model/";
	CreateModes modesCreator3;
	file_manager.loadModel(outfolder, modesCreator3.nmodel);
	cout<<modesCreator3.nmodel.F[0].rows()<<" "<<modesCreator3.nmodel.F[1].rows()<<endl;
	string folder="/home/ikriuchevs/workspace/melted/tests/Ndim/Validation/";
	Matrix M(modesCreator3.nmodel.F[0].rows(),modesCreator3.nmodel.F[1].rows());
	Matrix Mfit(modesCreator3.nmodel.F[0].rows(),modesCreator3.nmodel.F[1].rows());

//	cout<<modesCreator3.nmodel.params[0]<<endl;
//	cout<<modesCreator3.nmodel.F[1]<<endl;
//	cout<<modesCreator3.nmodel.F[1]<<endl;
	double start = omp_get_wtime();
	double fin;

	FILE * outFile;
	string outFileName=folder+"/nmax"+std::to_string( modesCreator3.nmodel.nbrModes)+".omp"+std::to_string(nthreads)+".dat";
	outFile = fopen(outFileName.c_str(), "w");
	fprintf (outFile, "#nbr of modes = %d\n",modesCreator3.nmodel.nbrModes);
	fprintf (outFile, "#Folder Error Time CudaError CudaTime\n");
	blaze::setNumThreads( nthreads );
	//CUDA
	modesCreator3.cudaInitND(20);
	cv::Mat res;
	res.create(modesCreator3.nmodel.F[0].rows(),modesCreator3.nmodel.F[1].rows(),opencvMatrixType);

	for(int i=1;i<=8;++i)
	{
		string filename=folder+std::to_string(i)+c_odb_filename;//<--odb_filename defined in
//		cout<<filename<<endl;
		file_manager.readODB2(filename, M);

		string paramname=folder+std::to_string(i)+"/params.dat";
//		cout<<paramname<<endl;
		std::ifstream inFile;
		inFile.open(paramname, std::ios_base::in);
		double temp=0;
		Vector param(dim-2,0);
		for(int d=0;d<dim-2;++d)
		{
			inFile>>temp;
			if(d>=0)
			{
				param[dim-1-2-(d)]=temp;
			}
		}

		inFile.close();
		start = omp_get_wtime();
		modesCreator3.fitNewND(param,Mfit);
		fin=omp_get_wtime()-start;

		double err=fabs(100.0*(1.0-blaze::norm(M)/blaze::norm(Mfit)));

		cout<<blaze::trans(param)<<" : OMP: error: "<<err<<"%; time = "<<fin<<" sec;";
		fprintf (outFile, "%d %f %lf ", i, err, fin );

		start = omp_get_wtime();
		modesCreator3.fitNewNDCuda(param,res);
		fin = omp_get_wtime()-start;

		opencv2blaze(res,Mfit);

		err=fabs(100.0*(1.0-blaze::norm(M)/blaze::norm(Mfit)));
		std::cout<<" CUDA: error="<<err<<"%; time = "<<fin  << " sec"<<std::endl;
		fprintf (outFile, "%f %lf\n", err, fin );

	}
	fclose (outFile);



	odb_finalizeAPI();

	return(0);
}




