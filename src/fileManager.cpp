#include <stdio.h>
#include <stdlib.h> /* strtof */
#include "fileManager.h"
//#include "NDMatrix.h"
#include <odb_API.h>

#include <boost/filesystem.hpp>
using namespace boost::filesystem;
#include <fstream>
#include <iostream>
#include <sstream>      // std::istringstream
//---------------------------------------------------------------------------------------------------

typedef blaze::DynamicVector<string,blaze::rowMajor> SVector;
const int c_tmax=100;
const string c_field_string="NT11";
const string c_step_name="Step-1";
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
FileManager::FileManager(){	odb_initializeAPI();}
//---------------------------------------------------------------------------------------------------

FileManager::~FileManager(){odb_finalizeAPI();}
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
void FileManager::readParams(string foldername, int& dim, intVector& sizes,std::vector<Vector>& params)
{
	std::ifstream inFile;
	inFile.open(foldername+"/"+c_all_parameters_filename, std::ios_base::in);
	std::string line;
	dim=0;
	while (std::getline(inFile, line))
	{
		dim++;
	    std::istringstream iss(line);
	    double temp;
	    Vector p1;
	    int s=0;
	    while ((iss >> temp))
	    {

	    	s++;
	    	p1.resize(s,true);
	    	p1[s-1]=temp;
//	    	std::cout << temp<<" ";
	    } // error

	    params.push_back(p1);
//	    std::cout <<"\n";
	    // process pair (a,b)
	}
	sizes=intVector(dim);
	for(int d=0;d<dim;++d)
		sizes[d]=params[d].size();
	inFile.close();

}

//---------------------------------------------------------------------------------------------------
void FileManager::readFolder(string foldername, NinputData3& Ndata)
{
	int dim=0;						//dimensionality of the system
	intVector paramSizes;			//the sizes of parameters dimensions
//	std::vector<Vector> parameters;
	readParams(foldername, dim, paramSizes, Ndata.params);

	dim=dim+2;

	readFolder( foldername, Ndata, dim, paramSizes );
}
//---------------------------------------------------------------------------------------------------

void FileManager::readFolder(string foldername, NinputData3& Ndata, int dim, intVector param_sizes)
{
	//SORT parameters array++++++++!!!
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
//		std::cout<<fname<<std::endl;
		std::cout << "\r" << k << "/"<<Ndata.A.m_paramSize<<std::flush;
//		std::cout.flush();
		readODB2(fname, M);
		Ndata.A.setNext(M);
	}
	std::cout<<"\rDone           "<<std::endl;
	//--------------------------------------------

//	std::cout<<"2\n";
}
//---------------------------------------------------------------------------------------------------

void FileManager::readODB2(string filename, Matrix& m)
{
	odb_Odb& myOdb = openOdb(filename.c_str(), true);

//	odb_Step& step = myOdb.steps()[c_step_name];
	odb_String stepName(c_step_name.c_str(),c_step_name.size());
	odb_StepRepository& steps = myOdb.steps();
	odb_Step& step = steps[stepName];
	odb_SequenceFrame& allFramesInStep = step.frames();

	for(int timeID =0;timeID<m.columns();++timeID)
	{
		const odb_SequenceFieldValue& temp =
				allFramesInStep[timeID].fieldOutputs()[c_field_string.c_str()].values();
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
}
//---------------------------------------------------------------------------------------------------
void FileManager::readODB_SpacexTime(string filename, int& spaceDegreOfFreedom, int& timeDegreOfFreedom)
{
//	odb_initializeAPI();
//
	odb_Odb& myOdb = openOdb(filename.c_str(), true);
	odb_StepRepository& steps = myOdb.steps();
	odb_String stepName(c_step_name.c_str(),c_step_name.size());
	odb_Step& step = steps[stepName];

//	odb_Step& step = myOdb.steps()["Step-1"];
	odb_SequenceFrame& allFramesInStep = step.frames();
	int numFrames = allFramesInStep.size();
	numFrames=c_tmax;
	timeDegreOfFreedom=numFrames;
	int numValues = allFramesInStep[numFrames-1].fieldOutputs()[c_field_string.c_str()].values().size();
	spaceDegreOfFreedom=numValues;
	myOdb.close();
//	odb_finalizeAPI();
}

