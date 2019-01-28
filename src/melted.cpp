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



void cudaInit(int rows,int cols);
void fitNew(const double& newParam1, cv::Mat &result,const cv::Mat &F1, const cv::Mat &F2, const Matrix &F3, const Vector& param1);
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
	std::clock_t start;
	std::clock_t fin;
	//Reading snapshots
	//-----------------------------------------------------------------------
	FileManager file_manager;
	std::cout<<"Reading files: ";
	start = std::clock();
	file_manager.readFolder(snapshots_dir);
	inputData &inData=file_manager.inData;
	inData.error=0.00001;
	std::cout<<( std::clock() - start ) / (double) CLOCKS_PER_SEC<<" sec"<<"\n";


	//Create model
	//-----------------------------------------------------------------------
	std::cout<<"Creating modes: ";
	start = std::clock();
	CreateModes modesCreator(inData,nmax);
	file_manager.saveModel(outFolderName,modesCreator.model);

	std::cout<<( std::clock() - start ) / (double) CLOCKS_PER_SEC<<" sec"<<"\n";
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
	//		TESTconstNumberOfmodes(100,snapshots_dir,valJobsForlder,outFolder);
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

#include <sys/stat.h>
using namespace std;
/*
 ***************
utility functions
 ***************
 */
bool fileExists(const odb_String  &string);
void rightTrim(odb_String  &string,const char* char_set);
void printExecutionSummary();
/***************************************************************/

using namespace std;
int mainOdb(int argc, char **argv)
{


	//  odb_String elsetName;
	//  bool ifElset = false;
	//  odb_Set myElset;
	//  odb_String region = "over the entire model";
	//  char msg[256];
	//  char *abaCmd = argv[0];
	odb_initializeAPI();
	odb_String odbPath;
	odbPath = "/home/ikriuchevs/workspace/melted/tests/369/Snapshots/JobsGarbage/Conductivity_250/heatTransfer_Conductivity_250.odb";
	FileManager fm;
	Matrix m;
	fm.readODB(odbPath.CStr(),m);
	std::string folder_path = "/home/ikriuchevs/workspace/melted/tests/369/Snapshots/JobsGarbage/";
	fm.readFolder(folder_path);
	fm.write("temp.odb.dat",fm.inData.A[0]);

	//	fm.readFolder(folder_path);
	//	fm.write("temp.txt.dat",fm.inData.A[0]);
	exit(1);
	//	odbPath = "/home/ikriuchevs/workspace/melted/tests/369/Snapshots/TxtOutPutFiles/Conductivity_205.0/viewer_tutorial.odb";
	bool ifOdbName = false;
	if (!fileExists(odbPath))
	{
		cerr << "**ERROR** output database  " << odbPath.CStr()
                    								 << " does not exist\n" << endl;
		exit(1);
	}
	ifOdbName = true;
	cout<<"HEloo"<<endl;
	//  for (int arg = 0; arg<argc; arg++)
	//    {
	//      if (strncmp(argv[arg],"-o**",2) == 0)
	//        {
	//          arg++;
	//          odbPath = "tests/369/Snapshots/TxtOutFiles/conductivity_205.0/heatTransfer_Conductivity_205.0.odb";
	//          rightTrim(odbPath,".odb");
	//          if (!fileExists(odbPath))
	//            {
	//              cerr << "**ERROR** output database  " << odbPath.CStr()
	//                   << " does not exist\n" << endl;
	//              exit(1);
	//            }
	//          ifOdbName = true;
	//        }
	//      else if (strncmp(argv[arg],"-e**",2)== 0)
	//        {
	//          arg++;
	//          elsetName = argv[arg];
	//          ifElset = true;
	//        }
	//      else if (strncmp(argv[arg],"-h**",2)== 0)
	//        {
	//          printExecutionSummary();
	//          exit(0);
	//        }
	//    }
	//  if (ifOdbName)
	//    {
	//      cerr << "**ERROR** output database name is not provided\n";
	//      printExecutionSummary();
	//      exit(1);
	//    }
	// Open the output database

	cout<<"Opening .odb"<<endl;

	odb_Odb& myOdb = openOdb(odbPath, true);

	//	odb_Odb& myOdb = Odb(odbPath);
	cout<<"Accessing rootAssemly"<<endl;
	odb_Assembly& myAssembly = myOdb.rootAssembly();
	odb_InstanceRepositoryIT instIter(myAssembly.instances());
	for (instIter.first(); !instIter.isDone(); instIter.next())
		cout << instIter.currentKey().CStr() << endl;

	cout << "Node set keys:" << endl;
	odb_SetRepositoryIT setIter( myAssembly.nodeSets() );
	for (setIter.first(); !setIter.isDone(); setIter.next())
		cout << setIter.currentKey().CStr() << endl;

	odb_InstanceRepository& iCon =
			myOdb.rootAssembly().instances();
	odb_Instance& instance = iCon["PART-1-1"];

	odb_StepRepositoryIT stepIter( myOdb.steps() );
	for (stepIter.first(); !stepIter.isDone();
			stepIter.next())
		cout << stepIter.currentKey().CStr() << endl;

	odb_Step& step = myOdb.steps()["Step-1"];
	odb_SequenceFrame& allFramesInStep = step.frames();
	int numFrames = allFramesInStep.size();
	odb_Frame& lastFrame = allFramesInStep[numFrames-1];

	odb_FieldOutputRepository& fieldOutputRep =
			lastFrame.fieldOutputs();
	odb_FieldOutputRepositoryIT fieldIter( fieldOutputRep );
	for (fieldIter.first(); !fieldIter.isDone(); fieldIter.next())
		cout << fieldIter.currentKey().CStr() << endl;

	for (fieldIter.first(); !fieldIter.isDone();
			fieldIter.next()) {
		odb_FieldOutput& field =
				fieldOutputRep[fieldIter.currentKey()];
		const odb_SequenceFieldValue& seqVal = field.values();
		const odb_SequenceFieldLocation& seqLoc =
				field.locations();
		cout << field.name().CStr() << " : " << field.description().CStr()
	        				<< endl;
		cout << "    Type: " << field.type() << endl;
		int numLoc = seqLoc.size();
		for (int loc = 0; loc<numLoc; loc++){
			cout << "Position: "<<seqLoc.constGet(loc).position();
		}
		cout << endl;
	}


	//	    const odb_SequenceFieldValue& displacements =
	//	    		lastFrame.fieldOutputs()["NT11"].values();
	int timeID=0;
	const odb_SequenceFieldValue& displacements =
			allFramesInStep[timeID].fieldOutputs()["NT11"].values();
	int numValues = displacements.size();
	int numComp = 0;
	for (int i=0; i<numValues; i++) {
		const odb_FieldValue val = displacements[i];
		cout << "Node = " << val.nodeLabel();
		const float* const NT11 = val.data(numComp);
		cout << ", T = ";
		for (int comp=0;comp<numComp;comp++)
			cout << NT11[comp] << "  ";
		cout << endl;
	}

	myOdb.close();
	odb_finalizeAPI();

	return(0);
}

//Ndimensional matrix test
//searching for F[dimId], the rest F[i] i!=dimId are considered to be known
//R-residuals matrix

void findF(int dimId, myMatrix& R, Vector* F);
void findApprox( myMatrix& R, Vector* F);
NModel HOPGD( const myMatrix& M, double ec);
int main(int argc, char * argv[])
{
//	int dim=3;
//
//	intVector sizes(dim,2);
//	sizes[0]=2;
//	sizes[1]=2;
//	sizes[2]=3;
//
//	int arr[12]={1,2,3,4,    1,4,9,16,  3,12,27,48};
//

	int dim=4;

	intVector sizes(dim,2);
	sizes[0]=2;
	sizes[1]=2;
	sizes[2]=3;
	sizes[3]=2;

	int arr[24]={1,2,3,4,    1,4,9,16,  3,12,27,48,
				 5,1,3,14,  21,4,9, 6,  8,11,17, 8};


	cout<<"Start!"<<endl;


	//Create and fill matrix
	myMatrix M(dim,sizes);
	int i=0;
	while(M.setNext(arr[i]))
	{
		i++;
	}
	//Print matrix M in linear fasion
	M.resetCurrentPosition();
	while(!M.ifEndPosition())
	{
		cout<<M.getCurrentVal()<<" ";
		M.next();
	}
	cout<<"\n";
	//test HOPGD sum routine (findF function)

	double eBc=0.000001;
	double ec=0.001;
	NModel model=HOPGD( M,ec);

	cout<<"\nLets see some F1!"<<endl;
	for(int ii=0; ii< model.F[2].rows();++ii)
	{
		for(int j=0; j< model.F[2].columns();++j)
		{
			cout<<model.F[2](ii,j)<<" ";
		}
		cout<<"\n";
	}
	cout<<"\nFinish!"<<endl;
	return(0);

}

NModel HOPGD( const myMatrix& M, double ec)
{

	const double eBc=0.000001;
	const int c_nmax=500;
	const int vmax=200;
	intVector sizes=M.sizes();
	int dim=M.dimensionality();

	NModel model(dim);
	myMatrix fna(dim, sizes);
	myMatrix R=M;
	Vector *F;
	F=new Vector[dim];


	for (int n=0;n<c_nmax;++n)
	{
		cout<<"    Distance:"<<R.dist(M)<<endl;
		if(R.dist(M) < ec )
		{
			cout<<"!!!!CONGRATULATIONS!!!!!\n";
			break;
		}
		//Initialial values
		Vector bta1=Vector(dim,0);
		Vector br1=Vector(dim,1);
		Vector bt1=Vector(dim,0);
		myMatrix f(dim, sizes);
		for(int d =0;d<dim;d++)
		{
			F[d]=Vector(sizes[d],1); //Must come from M or model
		}

		for (int v=0;v<vmax;++v)
		{
			bool e1=true;
			for(int m=0;m<dim;++m)
			{
				findF(m,R,F);
				bt1[m]=blaze::norm(F[m]);
				if(v==0)
				{
					br1[m]=bt1[m];
				}
				e1=e1&&((bt1[m]-bta1[m])/br1[m]<eBc);
			}
			if(e1)
			{
				findApprox( f, F);

				R=R-f;
				fna=fna+f;
				fna.resetCurrentPosition();
				while(!fna.ifEndPosition())
				{
					cout<<fna.getCurrentVal()<<" ";
					fna.next();
				}

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
	return model;
}
void findApprox( myMatrix& R, Vector* F)
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
		R.setElement(R.getCurrentPosition(),fMult);
		//		cout<<fMult<<" ";
		R.next();
	}
	//	cout<<"\n";

}
void findF(int dimId, myMatrix& R, Vector* F)
{
	R.resetCurrentPosition();
	F[dimId]=Vector(F[dimId].size(),1);
	Vector newF=Vector(F[dimId].size(),0);
	int dim=R.dimensionality();
	intVector idx(dim,0);

	while(!R.ifEndPosition())
	{
		R.getCurrentIdx(idx);
		double fMult=1;
		for(int d=0;d<dim;d++)
			fMult*=F[d][idx[d]];
		newF[idx[dimId]]+=R.getElement(idx)*fMult;
		R.next();
	}


	double normBF=1;
	for(int d=0;d<dim;d++)
	{
		if (d!=dimId)
			normBF*=(trans(F[d])*F[d]);
	}
	//	cout<<"normBf="<<normBF<<endl;
	F[dimId]=newF/normBF;
	//	cout<<"Fi={";
	//
	//	for (int i=0;i<newF.size();i++)
	//		cout<<F[dimId][i]<<", ";
	//
	//	cout<<"}\n";
}

