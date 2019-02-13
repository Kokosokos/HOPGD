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
#include <ctime>
#include <omp.h>



using namespace boost::filesystem;

using namespace std;


void validation(string valJobsFolder,string outFolderName,const CreateModes& modesCreator, int nModes)
{
	if (nModes > modesCreator.model.nbrModes || nModes < 1)
	{
		//	cout<<"Warning: nModes > total # of modes in the model.\n Forced to be the same.\n";
		nModes=modesCreator.model.nbrModes;
	}
	//	std::clock_t start;
	//	std::clock_t fin;
	double start;
	double fin;
	start = omp_get_wtime();

	//Validation
	//-----------------------------------------------------------------------
	Matrix M;
	string outFileName= outFolderName+"nmax."+std::to_string(nModes)+".omp8.dat";
	std::cout<<outFolderName<<std::endl;
	FILE * outFile;

	std::cout<<"Reading validation files: ";
	outFile = fopen(outFileName.c_str(), "w");

	//start = std::clock();
	start = omp_get_wtime();
	FileManager validation_files;

	validation_files.readFolder(valJobsFolder);
	std::cout<< omp_get_wtime() - start <<" sec"<<"\n";

	start = omp_get_wtime();
	fprintf (outFile, "#nbr of modes = %d\n", nModes);
	fprintf (outFile, "#value error time\n");

	for (int val_idx =0;val_idx<validation_files.inData.param1DegreOfFreedom;++val_idx)
	{
		//		start = std::clock();
		start = omp_get_wtime();
		modesCreator.fitNew(validation_files.inData.param1[val_idx],M, nModes);
		//		fin = std::clock();
		fin = omp_get_wtime();
		double err=fabs(100.0*(1.0-blaze::norm(M)/blaze::norm(validation_files.inData.A[val_idx])));

		std::cout<<"Validation: C="<<validation_files.inData.param1[val_idx]<<"; err="<<err<<"%; "<<( fin - start ) << "sec"<<std::endl;
		//		fprintf (outFile, "%.1f %f %f\n", validation_files.inData.param1[val_idx], err,( fin - start ) / (double) CLOCKS_PER_SEC);
		fprintf (outFile, "%.1f %f %lf\n", validation_files.inData.param1[val_idx], err,( fin - start ) );
	}

	fclose (outFile);
}
void validation(string valJobsFolder,string outFolderName,const CreateModes& modesCreator)
{
	validation(valJobsFolder, outFolderName, modesCreator, modesCreator.model.nbrModes);
}
void TESTconstNumberOfmodes(int nmax,string snapshots_dir,string valJobsFolder,string outFolderName, bool useValidation=true)
{
	double start;
	//Reading snapshots
	//-----------------------------------------------------------------------
	FileManager file_manager;
	std::cout<<"Reading files: ";
	start = omp_get_wtime();
	file_manager.readFolder(snapshots_dir);
	inputData &inData=file_manager.inData;
	inData.error=0.00001;
	std::cout<< omp_get_wtime() - start <<" sec"<<"\n";


	//Create model
	//-----------------------------------------------------------------------
	std::cout<<"Creating modes: ";
	start = omp_get_wtime();
	CreateModes modesCreator(inData,nmax);
	std::cout<< omp_get_wtime() - start <<" sec"<<"\n";
	file_manager.saveModel(outFolderName,modesCreator.model);

	std::cout<<"#modes = "<<modesCreator.model.nbrModes<<std::endl;
	//Validation
	//----------------------------------------------------------------------
	if(useValidation)
	{
		validation(valJobsFolder,outFolderName, modesCreator);
	}
}

void TESTonlinePart(string valJobsFolder, string outFolderName, int nModes=0)
{
	CreateModes modesCreator;
	FileManager file_manager;
	file_manager.loadModel(outFolderName, modesCreator.model);
	validation(valJobsFolder,outFolderName, modesCreator,nModes);
}

int mainCuda(int argc, char * argv[])
{
	int nthreadsInit=8;
	int nthreadsMax=8;
	//	int nnModes=7;
	//	int nModesArray[]={1,10,20,30,50,100,150};
	int nnModes=6;
	int nModesArray[]={1,10,20,30,50,90};
	//	int nnModes=1;
	//	int nModesArray[]={20};
	int nModes=50;

	cv::Mat F1cv;
	cv::Mat F2cv,F3cv,res;
	double start;
	double fin;


	string snapshots_dir ="/home/ikriuchevs/workspace/melted/tests/369/Snapshots/TxtOutPutFiles/";
	string valJobsFolder="/home/ikriuchevs/workspace/melted/tests/369/Snapshots/ValidationJobs/";
	string outFolder="/home/ikriuchevs/workspace/melted/tests/369/testsCpp/";

	snapshots_dir ="/home/ikriuchevs/workspace/melted/tests/11256/Snapshots/TxtOutPutFiles/";
	valJobsFolder="/home/ikriuchevs/workspace/melted/tests/11256/Snapshots/ValidationJobs/";
	outFolder="/home/ikriuchevs/workspace/melted/tests/11256/testsCpp/";
	//
	snapshots_dir ="/home/ikriuchevs/workspace/melted/tests/128320/Snapshots/TxtOutPutFiles/";
	valJobsFolder="/home/ikriuchevs/workspace/melted/tests/128320/Snapshots/ValidationJobs/";
	outFolder="/home/ikriuchevs/workspace/melted/tests/128320/testsCpp/";

	snapshots_dir ="/home/ikriuchevs/workspace/melted/tests/1069188/Snapshots/JobsGarbage/";
	valJobsFolder="/home/ikriuchevs/workspace/melted/tests/1069188/Snapshots/JobsGarbage/";
	outFolder="/home/ikriuchevs/workspace/melted/tests/1069188/testsCpp/";

	CreateModes modesCreator;
	FileManager file_manager;
	file_manager.loadModel(outFolder, modesCreator.model);


	FileManager validation_files;

	validation_files.readFolder(valJobsFolder);
	for (int nn=0;nn<nnModes;++nn)
	{
		nModes=nModesArray[nn];
		for (int nthreads=nthreadsInit;nthreads<=nthreadsMax;nthreads++)
		{

			modesCreator.cudaInit(nModes);

			blaze::setNumThreads( nthreads );

			string outFileName= outFolder+"nmax."+std::to_string(nModes)+".omp"+std::to_string(nthreads)+".cuda.dat";
			FILE * outFile;
			outFile = fopen(outFileName.c_str(), "w");
			fprintf (outFile, "#nbr of modes = %d\n", nModes);
			fprintf (outFile, "#value error time cudaError cudaTime\n");


			std::cout<<"#modes: "<<nModes<<std::endl;
			Matrix M,ress;
			cv::Mat Mcv;
			res.create(F1cv.rows,F2cv.rows,opencvMatrixType);

			for (int val_idx =0;val_idx<validation_files.inData.param1DegreOfFreedom;++val_idx)
			{
				std::cout<<"Validation: C="<<validation_files.inData.param1[val_idx]<<std::endl;

				start = omp_get_wtime();
				modesCreator.fitNew(validation_files.inData.param1[val_idx],M, nModes);
				fin = omp_get_wtime();

				double err=fabs(100.0*(1.0-blaze::norm(M)/blaze::norm(validation_files.inData.A[val_idx])));
				std::cout<<"err="<<err<<"%; "<<( fin - start ) << "sec"<<std::endl;
				fprintf (outFile, "%.1f %f %lf ", validation_files.inData.param1[val_idx], err,( fin - start ) );

				start = omp_get_wtime();
				modesCreator.fitNewCuda(validation_files.inData.param1[val_idx],res);
				fin = omp_get_wtime();

				opencv2blaze(res,M);

				err=fabs(100.0*(1.0-blaze::norm(M)/blaze::norm(validation_files.inData.A[val_idx])));
				std::cout<<"err="<<err<<"%; "<<( fin - start ) << "sec"<<std::endl;
				fprintf (outFile, "%f %lf\n", err,( fin - start ) );

			}
			fclose (outFile);
		}
	}
	std::cout<<"END!!!"<<std::endl;
	return 1;
}

int main3(int argc, char * argv[])
{

	blaze::setNumThreads( 8 );
	//Tests
	//---------------------------------------------------------------------------------------------------

	//#1: different number of modes
	string snapshots_dir ="/home/ikriuchevs/workspace/melted/tests/369/Snapshots/TxtOutPutFiles/";
	string valJobsForlder="/home/ikriuchevs/workspace/melted/tests/369/Snapshots/ValidationJobs/";
	string outFolder="/home/ikriuchevs/workspace/melted/tests/369/testsCpp/";
	TESTconstNumberOfmodes(100,snapshots_dir,valJobsForlder,outFolder);
	//	TESTonlinePart(valJobsForlder,outFolder,1);
	//	TESTonlinePart(valJobsForlder,outFolder,5);
	//	TESTonlinePart(valJobsForlder,outFolder,10);
	//	TESTonlinePart(valJobsForlder,outFolder,40);
	//	TESTonlinePart(valJobsForlder,outFolder,20);
	//	TESTonlinePart(valJobsForlder,outFolder,50);
	//	TESTonlinePart(valJobsForlder,outFolder,100);
	//---------------------------------------------------------------------------------------------------

	//#2: different number of reference snapshots


	//---------------------------------------------------------------------------------------------------
	//#3: Larger mesh
	snapshots_dir ="/home/ikriuchevs/workspace/melted/tests/11256/Snapshots/TxtOutPutFiles/";
	valJobsForlder="/home/ikriuchevs/workspace/melted/tests/11256/Snapshots/ValidationJobs/";
	outFolder="/home/ikriuchevs/workspace/melted/tests/11256/testsCpp/";

	//		TESTconstNumberOfmodes(100,snapshots_dir,valJobsForlder,outFolder);
	//	TESTonlinePart(valJobsForlder,outFolder,1);
	//	TESTonlinePart(valJobsForlder,outFolder,5);
	//	TESTonlinePart(valJobsForlder,outFolder,10);
	//		TESTonlinePart(valJobsForlder,outFolder,40);
	//	TESTonlinePart(valJobsForlder,outFolder,20);
	//	TESTonlinePart(valJobsForlder,outFolder,50);
	//	TESTonlinePart(valJobsForlder,outFolder,100);
	snapshots_dir ="/home/ikriuchevs/workspace/melted/tests/128320/Snapshots/TxtOutPutFiles/";
	valJobsForlder="/home/ikriuchevs/workspace/melted/tests/128320/Snapshots/ValidationJobs/";
	outFolder="/home/ikriuchevs/workspace/melted/tests/128320/testsCpp/";
	//		TESTconstNumberOfmodes(100,snapshots_dir,valJobsForlder,outFolder);
	//TESTconstNumberOfmodes(500,snapshots_dir,valJobsForlder,outFolder);

	//		TESTonlinePart(valJobsForlder,outFolder,1);
	//	TESTonlinePart(valJobsForlder,outFolder,2);
	//	TESTonlinePart(valJobsForlder,outFolder,3);
	//	TESTonlinePart(valJobsForlder,outFolder,4);
	//	TESTonlinePart(valJobsForlder,outFolder,5);
	//	TESTonlinePart(valJobsForlder,outFolder,10);
	//		TESTonlinePart(valJobsForlder,outFolder,20);
	//		TESTonlinePart(valJobsForlder,outFolder,40);
	//	TESTonlinePart(valJobsForlder,outFolder,50);
	//	TESTonlinePart(valJobsForlder,outFolder,100);

	snapshots_dir ="/home/ikriuchevs/workspace/melted/tests/1069188/Snapshots/JobsGarbage/";
	valJobsForlder="/home/ikriuchevs/workspace/melted/tests/1069188/Snapshots/JobsGarbage/";
	outFolder="/home/ikriuchevs/workspace/melted/tests/1069188/testsCpp/";
	TESTconstNumberOfmodes(100,snapshots_dir,valJobsForlder,outFolder);

	cout << "!!!End!!!" << endl;
	return 0;
}


#include <odb_API.h>

int main(int argc, char * argv[])
{
	odb_initializeAPI();

	int dim=6;
	int nthreads=8;
	FileManager file_manager;
	NinputData3 ndata;
	intVector sizes(dim-2,0);
	//Small test
//	std::cout<<"Small TEST 1-18: ";
//	//Reverse direction:(
//	sizes[0]=3; //Temperature
//	sizes[1]=3; //Specific heat
//	sizes[2]=2; //Density
//	sizes[3]=1; //Density
//	string folder="/home/ikriuchevs/workspace/melted/tests/Ndim/Jobs/";
//	string outfolder="/home/ikriuchevs/workspace/melted/tests/Ndim/Jobs/118model/";
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
//	CreateModes modesCreator2(ndata);
//	std::cout<< omp_get_wtime() - start <<" sec"<<"\n";
//
//	std::cout<<"Saving the model: \n";
//	file_manager.saveModel(outfolder, modesCreator2.nmodel);



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
		string filename=folder+std::to_string(i)+"/just.odb";
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




