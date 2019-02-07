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
#include "odb_MaterialTypes.h"
using namespace std;
int main327532(int argc, char **argv)
{


	//  odb_String elsetName;
	//  bool ifElset = false;
	//  odb_Set myElset;
	//  odb_String region = "over the entire model";
	//  char msg[256];
	//  char *abaCmd = argv[0];
	odb_initializeAPI();
	odb_String odbPath;
	odbPath = "/home/ikriuchevs/workspace/melted/tests/Ndim/Jobs/1/just.odb";

	bool ifOdbName = false;
//	if (!fileExists(odbPath))
//	{
//		cerr << "**ERROR** output database  " << odbPath.CStr()
//                    																										 << " does not exist\n" << endl;
//		exit(1);
//	}
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

	odb_Material mymat;
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
//	const odb_SequenceFieldValue& displacements =
//			allFramesInStep[timeID].fieldOutputs()["NT11"].values();
//	int numValues = displacements.size();
//	int numComp = 0;
//	for (int i=0; i<numValues; i++) {
//		const odb_FieldValue val = displacements[i];
//		cout << "Node = " << val.nodeLabel();
//		const float* const NT11 = val.data(numComp);
//		cout << ", T = ";
//		for (int comp=0;comp<numComp;comp++)
//			cout << NT11[comp] << "  ";
//		cout << endl;
//	}

	myOdb.close();
	odb_finalizeAPI();

	return(0);
}

//Ndimensional matrix test
//searching for F[dimId], the rest F[i] i!=dimId are considered to be known
//R-residuals matrix


#include "mymatrix3.h"
void findApprox( myMatrix3& R, Vector* F);
void findF3(const int& dimId, myMatrix3& R, Vector* F);

NModel HOPGD( const myMatrix3& M, double ec)
{
	const int c_nmax=100;
	double start,fin;
//	start = omp_get_wtime();
	const double eBc=0.000001;
	const int vmax=200;
	intVector sizes=M.sizes();
	int dim=M.dimensionality();

	NModel model(dim);
	myMatrix3 fna(dim, sizes);
	myMatrix3 f(dim, sizes);
	myMatrix3 R=M;
	Vector *F;
	F=new Vector[dim];
//	std::cout<< "Prep time: "<<omp_get_wtime()-start<<std::endl;
	double findFtime=0;
	double findFnorm=0;
	cout<<M.m_values[0]<<endl;
	for (int n=0;n<c_nmax;++n)
	{
		double normBF=1;
		double prevNorm=1;

		double start,fin;
		std::cout<<n<<"\n";
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
			bool e1=true;

//			start = omp_get_wtime();
			for(int m=0;m<dim;++m)
			{

				double currNorm=blaze::sqrNorm(F[m]);
				normBF/=currNorm;
				normBF*=prevNorm;
				findF3(m,R,F);
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
				cout<<"=============="<<n<<"=========="<<endl;
				cout<<M.m_values[0]<<endl;

				R=R-f;
				fna=fna+f;
				cout<<fna.m_values[0]<<endl;


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
//	std::cout<<"findFtime: "<<findFtime<<std::endl;
//	std::cout<<"findFnorm: "<<findFnorm<<std::endl;
//	nmodel=model;

}

//---------------------------------------------------------------------------------------------------

void findF(const int& dimId, myMatrix& R, Vector* F)
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

		for(int d=0;d<dim;d++)
			fMult*=F[d][idx[d]];
//
//		for(int d=0;d<dimId;++d)
//			fMult*=F[d][idx[d]];
//		for(int d=dimId+1;d<dim;++d)
//			fMult*=F[d][idx[d]];
		newF[idx[dimId]]+=R.getCurrentVal()*fMult;
		R.next();
	}
	F[dimId]=newF;
}

void findF3(const int& dimId, myMatrix3& R, Vector* F)
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
void findApprox( myMatrix3& R, Vector* F)
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


int main234(int argc, char * argv[])
{

	//initialize and fill test matrix 3x3


	unsigned int dim=4;
	intVector sizes(dim,dim);
	sizes[0]=5;
	sizes[1]=4;
	sizes[2]=3;
	sizes[3]=2;
	Matrix m(sizes[0],sizes[1],0);
	for (int i=0;i<sizes[0];++i)
	{
		for (int j=0;j<sizes[1];++j)
		{
			m(i,j)=i+i*j+1; cout<<m(i,j)<<" ";
		}
		cout<<endl;
	}


	myMatrix3 m3(dim,sizes);
	m3.setNext(m);
	m3.setNext(3*m);
	m3.setNext(5*m);

	m3.setNext(2*m);
	m3.setNext(4*m);
	m3.setNext(6*m);

	NModel nm = HOPGD( m3, 0.01);
	FileManager fmmm;
	fmmm.saveModel("/home/ikriuchevs//workspace/melted/tests/Ndim/testNmodelSave/", nm);
	return(0);


	//	blaze::setNumThreads( 4 );
	//	omp_set_num_threads(4);
	NinputData input;
	FileManager fm;
	myMatrix M;
	//Reading snapshots
	string snapshots_dir ="/home/ikriuchevs/workspace/melted/tests/369/Snapshots/TxtOutPutFiles/";
	string valJobsFolder="/home/ikriuchevs/workspace/melted/tests/369/Snapshots/ValidationJobs/";
	string outFolder="/home/ikriuchevs/workspace/melted/tests/369/testsCpp/ND/";
	string outFolder2="/home/ikriuchevs/workspace/melted/tests/369/testsCpp/ND2/";

	//	snapshots_dir ="/home/ikriuchevs/workspace/melted/tests/128320/Snapshots/TxtOutPutFiles/";
	//	valJobsFolder="/home/ikriuchevs/workspace/melted/tests/128320/Snapshots/ValidationJobs/";
	//	outFolder="/home/ikriuchevs/workspace/melted/tests/128320/testsCpp/ND/";
	//	outFolder2="/home/ikriuchevs/workspace/melted/tests/128320/testsCpp/ND2/";

	//	snapshots_dir ="/home/ikriuchevs/workspace/melted/tests/11256/Snapshots/TxtOutPutFiles/";
	//	valJobsFolder="/home/ikriuchevs/workspace/melted/tests/11256/Snapshots/ValidationJobs/";
	//	outFolder="/home/ikriuchevs/workspace/melted/tests/11256/testsCpp/ND/";


	double start = omp_get_wtime();
	std::cout<<"Reading Files: ";
	fm.readFolder(snapshots_dir, input);

	FileManager file_manager;
	std::cout<<"Reading files: ";
	file_manager.readFolder(snapshots_dir);
	inputData &inData=file_manager.inData;

	std::cout<< omp_get_wtime() - start <<" sec"<<"\n";

	//HOPGD
	input.error=0.00001;
	input.nModesMax=200;

	//Creating modes

	start = omp_get_wtime();
//	std::cout<<"HOPGD: ";
//	CreateModes modesCreator(input);
//	std::cout<< omp_get_wtime() - start <<" sec"<<"\n";
//	std::cout<<"#modes = "<<modesCreator.model.nbrModes<<std::endl;
//	fm.saveModel(outFolder, modesCreator.model);



	start = omp_get_wtime();
	std::cout<<"HOPGD3: ";
	std::cout<<inData.param1DegreOfFreedom<<" "<<inData.spaceDegreOfFreedom<<endl;
	NinputData3 input3(inData);
	input3.error=input.error;
	input3.nModesMax=input.nModesMax;
	CreateModes modesCreator3(input3);
	std::cout<< omp_get_wtime() - start <<" sec"<<"\n";
	std::cout<<"#modes = "<<modesCreator3.model.nbrModes<<std::endl;
	fm.saveModel(outFolder, modesCreator3.model);

	//Validation
	//	validation(valJobsFolder,outFolder,modesCreator, input.nModesMax);
	//	outFolder="/home/ikriuchevs/workspace/melted/tests/369/testsCpp/ND/";
//	TESTconstNumberOfmodes(input.nModesMax,snapshots_dir,valJobsFolder,outFolder2,false);

	cout<<"\nFinish!"<<endl;
	return(0);

}
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
	fprintf (outFile, "#value error time\n");
	blaze::setNumThreads( nthreads );
	for(int i=1;i<=8;++i)
	{
		string filename=folder+std::to_string(i)+"/just.odb";
		cout<<filename<<endl;
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

				std::cout<<temp<<" ";
			}
		}

		inFile.close();
		start = omp_get_wtime();
		modesCreator3.fitNewND2(param,Mfit);
		fin=omp_get_wtime()-start;

		double err=fabs(100.0*(1.0-blaze::norm(M)/blaze::norm(Mfit)));

		cout<<" error: "<<err<<"%; time = "<<fin<<" sec"<<endl;

		fprintf (outFile, "%d %f %lf\n", i, err, fin );


	}
	fclose (outFile);

	sizes[0]=8; //Temperature
	sizes[1]=1; //Specific heat
	sizes[2]=1; //Density
	sizes[3]=1; //Density


	odb_finalizeAPI();

	return(0);
}




