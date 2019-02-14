#include <stdio.h>
#include <stdlib.h> /* strtof */
#include "fileManager.h"
//#include "NDMatrix.h"

#include <boost/filesystem.hpp>
using namespace boost::filesystem;
#include <fstream>

#include <iostream>
#include <sstream>      // std::istringstream
#include <odb_API.h>
//using std::string;
typedef blaze::DynamicVector<string,blaze::rowMajor> SVector;
const int c_tmax=100;

//---------------------------------------------------------------------------------------------------
class sort_indices
{
private:
	Vector mparr;
public:
	sort_indices(Vector parr) : mparr(parr) {}
	bool operator()(int i, int j) const { return mparr[i]<mparr[j]; }
};

void vectorSort2(Vector& p,SVector& s)
{
	int *indices=new int[p.size()];
	for (uint i =0;i<p.size();++i)
		indices[i]=i;
	std::sort(indices,indices+p.size(),sort_indices(p));
	Vector p1=p;
	SVector s1=s;;
	for (uint i =0;i<p.size();++i)
	{
		p[i]=p1[indices[i]];
		s[i]=s1[indices[i]];
	}
	delete[] indices;

}
//---------------------------------------------------------------------------------------------------

double my_stod (std::string const& s) {
    std::istringstream iss (s);
    iss.imbue (std::locale("C"));
    double d;
    iss >> d;
    // insert error checking.
    return d;
}

//---------------------------------------------------------------------------------------------------
void FileManager::read(string filename, Matrix& m)
{

	std::ifstream inFile;
	inFile.open(filename.c_str(), std::ios_base::in);
	double f;
	int xsize=m.columns();
	int ysize=m.rows();
	int k=0;
	if (!inFile) {
			perror ("Error opening file");
		}
	else
	{
		while ( inFile >> m(k/xsize,k-(k/xsize)*xsize) && k!=xsize*ysize )
		{
			++k;
		}
	}
	inFile.close();
}
//---------------------------------------------------------------------------------------------------
void FileManager::write(string filename, const Matrix& m)
{
	std::ofstream outFile;
	outFile.open(filename.c_str(), std::ios_base::out);
	int xsize=m.columns();
	int ysize=m.rows();

	if (!outFile) perror ("Error opening file");
	else
	{
		for (int i=0;i<ysize;++i)
		{
			for (int j=0;j<xsize;++j)
			{
				outFile<<m(i,j)<<" ";
			}
			outFile<<"\n";
		}
	}
	outFile.close();
}

FileManager::~FileManager()
{
	//delete[] inData.A;
}
//---------------------------------------------------------------------------------------------------

void FileManager::saveModel(string dirPath, const NModel& model)
{

	std::ofstream outFile;
	outFile.open((dirPath+"model.dat").c_str(), std::ios_base::out);
	outFile<<model.nbrModes<<" ";
	outFile<<"\n";
	outFile<<model.dim;
	outFile<<"\n";
	for (int d=0;d<model.F.size();++d)
		outFile<<model.F[d].rows()<<" ";
	outFile<<"\n";

	for (int i=0;i<model.params.size();++i)
	{
		for (int j=0;j<model.params[i].size();++j)
			outFile<<model.params[i][j]<<" ";
		outFile<<"\n";
	}
	outFile.close();
	for (int d=0;d<model.F.size();++d)
		write(dirPath+"/F"+std::to_string(d)+".dat", model.F[d]);

}
//---------------------------------------------------------------------------------------------------

bool FileManager::loadModel(string dirPath, NModel& model)
{
	std::ifstream inFile;
	inFile.open((dirPath+"model.dat").c_str(), std::ios_base::in);
	if (!inFile) {
		std::cerr << "Unable to open file "+(dirPath+"model.dat");
		return false;   // call system to stop
	}
	inFile>>model.nbrModes;
	inFile>>model.dim;

	intVector sizes(model.dim);
//	F.resize(model.dim);
	for (int d=0;d<model.dim;++d)
	{
		inFile>>sizes[d];
		Matrix M(sizes[d],model.nbrModes,0);
		model.F.push_back(M);
	}

	for (int d=2;d<model.dim;++d)
	{
		Vector v(sizes[d],0);
		for (int j=0;j<sizes[d];++j)
			inFile>>v[j];
		model.params.push_back(v);
	}
	inFile.close();
	for (int d=0;d<model.F.size();++d)
		read(dirPath+"/F"+std::to_string(d)+".dat", model.F[d]);
	return true;

}

//---------------------------------------------------------------------------------------------------
//Reads parameters
void FileManager::readParams(string foldername, NinputData3& Ndata, int dim, intVector sizes)
{
	std::vector<Vector> allParams;
	for (int k=1;k<Ndata.A.m_paramSize+1;++k)
	{
		string folder=foldername+"/"+std::to_string(k)+"/";
		string fname=folder+c_parameters_filename;
		std::ifstream inFile;
		inFile.open(fname, std::ios_base::in);
		double temp=0;
		Vector param(dim-2,0);

		for(int d=0;d<dim-2;++d)
		{
			inFile>>temp;
			if(d>=0)
			{
				param[d]=temp;
			}
		}
		allParams.push_back(param);
		//		std::cout<<"\n";
		inFile.close();
	}
	//Rearrange parameters (fastest first)
	int jump=1;//Ndata.A.m_matrix_size;
	for(int d=2;d<dim;++d)
	{
		Vector v(Ndata.A.sizes()[d]);
		for(int j=0;j<Ndata.A.sizes()[d];j++)
		{
			v[j]=allParams[j*jump][dim-1-(d)];
		}
		Ndata.params.push_back(v);
		jump*=Ndata.A.sizes()[d];
	}
}
//---------------------------------------------------------------------------------------------------

//Makes sure parameters are sorted and files are read in right order
//Creates EMPTY myMatrix

//---------------------------------------------------------------------------------------------------
void FileManager::readFolder(string foldername, NinputData3& Ndata, int dim, intVector param_sizes)
{
	//!SORT parameters array++++++++!!!
//	vectorSort2(Ndata.params[0],sparam1);
	//+++++++++++++++++++++++++++
	int spaceDegreOfFreedom=0;
	int timeDegreOfFreedom=0;
	string new_filename = foldername + "/1/"+c_odb_filename;
	readODB_SpacexTime(new_filename,spaceDegreOfFreedom,timeDegreOfFreedom);

	intVector sizes(dim,0);
	sizes[0]=spaceDegreOfFreedom;
	sizes[1]=timeDegreOfFreedom;
	for(int d=2;d<dim;++d)
	{
		sizes[d]=param_sizes[d-2];
	}
	Ndata.A.setSize(dim,sizes);
	Ndata.A.resetCurrentPosition();
	Matrix M(spaceDegreOfFreedom,timeDegreOfFreedom);

	std::vector<Vector> allParams;
	for (int k=1;k<Ndata.A.m_paramSize+1;++k)
	{

		string folder=foldername+"/"+std::to_string(k)+"/";
		string fname=folder+c_odb_filename;
		readODB2(fname, M);
		Ndata.A.setNext(M);
		fname=folder+"params.dat";
		std::ifstream inFile;
		inFile.open(fname, std::ios_base::in);
		double temp=0;
		Vector param(dim-2,0);
//		std::cout<<"params: "<<k<<" ";
//		for(int d=0;d<dim-2;++d)
//		{
//			inFile>>temp;
//			param[d]=temp;
//			Ndata.params.push_back(param);
//			std::cout<<temp<<" ";
//		}
		//TEMPORARY?
		for(int d=0;d<dim-2;++d)
		{
			inFile>>temp;
			if(d>=0)
			{
				param[d]=temp;

//				std::cout<<temp<<" ";
			}
		}
		allParams.push_back(param);
//		std::cout<<"\n";
		inFile.close();

	}
//	std::reverse(allParams.begin(), allParams.end());
//	std::cout<<"1\n";
	int jump=1;//Ndata.A.m_matrix_size;
	for(int d=2;d<dim;++d)
	{
		Vector v(Ndata.A.sizes()[d]);
		for(int j=0;j<Ndata.A.sizes()[d];j++)
		{
			v[j]=allParams[j*jump][dim-1-(d)];
		}
		Ndata.params.push_back(v);
		jump*=Ndata.A.sizes()[d];
	}
	//--------------------------------------------

//	std::cout<<"2\n";
}
//---------------------------------------------------------------------------------------------------

void FileManager::readODB2(string filename, Matrix& m)
{

//	odb_initializeAPI();

	odb_Odb& myOdb = openOdb(filename.c_str(), true);

	odb_Step& step = myOdb.steps()["Step-1"];
	odb_SequenceFrame& allFramesInStep = step.frames();

	for(int timeID =0;timeID<m.columns();++timeID)
	{
		const odb_SequenceFieldValue& temp =
				allFramesInStep[timeID].fieldOutputs()["NT11"].values();
		int numComp = 0;
		for (int spaceID=0; spaceID<m.rows(); spaceID++)
		{
			const odb_FieldValue val = temp[spaceID];
			const float* const NT11 = val.data(numComp);
			int comp=0;
//			for (int comp=0;comp<numComp;comp++)
			m(spaceID,timeID)= NT11[comp];
		}
//	cout<<timeID<<endl;

	}
	myOdb.close();
//	odb_finalizeAPI();
}
//---------------------------------------------------------------------------------------------------
void FileManager::readODB_SpacexTime(string filename, int& spaceDegreOfFreedom, int& timeDegreOfFreedom)
{
//	odb_initializeAPI();
//
	odb_Odb& myOdb = openOdb(filename.c_str(), true);
	odb_StepRepository& steps = myOdb.steps();
	odb_String stepName("Step-1",6);
	odb_Step& step = steps[stepName];

//	odb_Step& step = myOdb.steps()["Step-1"];
	odb_SequenceFrame& allFramesInStep = step.frames();
	int numFrames = allFramesInStep.size();
	numFrames=c_tmax;
	timeDegreOfFreedom=numFrames;
	int numValues = allFramesInStep[numFrames-1].fieldOutputs()["NT11"].values().size();
	spaceDegreOfFreedom=numValues;
	myOdb.close();
//	odb_finalizeAPI();
}

