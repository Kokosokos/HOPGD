#include <stdio.h>
#include <stdlib.h> /* strtof */
#include "fileManager.h"
#include "mymatrix.h"

#include <boost/filesystem.hpp>
using namespace boost::filesystem;
#include <fstream>

#include <iostream>
#include <sstream>      // std::istringstream
#include <odb_API.h>
//using std::string;
typedef blaze::DynamicVector<string,blaze::rowMajor> SVector;
const int c_tmax=100;
//Stupid sorting folder names routine

double my_stod (std::string const& s) {
    std::istringstream iss (s);
    iss.imbue (std::locale("C"));
    double d;
    iss >> d;
    // insert error checking.
    return d;
}

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

//Add Warnings/Errors
//if nFiles is different in folders
//if some of the files are empty

void FileManager::readFolderOld(string foldername)
{
	path p = path(foldername);
	directory_iterator it(p);
	string filenameString="results_";
	int nFiles=0;
	inData.param1DegreOfFreedom=0;
	inData.timeDegreOfFreedom=0;

	inData.param1DegreOfFreedom=std::count_if(directory_iterator(p), directory_iterator(), static_cast<bool(*)(const path&)>(is_directory));
	inData.param1.resize(inData.param1DegreOfFreedom);
	SVector sparam1(inData.param1DegreOfFreedom);
	int k=0;

	string folderBaseName;
	while (it != directory_iterator{})
		if(is_directory(it->path()))
		{
			string folder=it->path().generic_string();
			std::size_t	found=folder.rfind("_");
			folderBaseName=folder.substr(0,found+1);
			string folderParam1=folder.substr(found+1);
			inData.param1[k]=my_stod(folderParam1);
			sparam1[k]=folderParam1;
//			std::cout<<k <<": "+folderParam1<<"/"<< inData.param1[k]<<" "<<std::stof(folderParam1,&sz)<<"\n";
			directory_iterator it_files(it->path());
			nFiles=0;
			while (it_files != directory_iterator{})
			{
				string file=basename(it_files->path());
				//std::cout<<file<<" "<<file.compare(0,filenameString.size(),filenameString)<<'\n';
				it_files++;
				if(!file.compare(0,filenameString.size(),filenameString))
					nFiles++;
			}
			it++;
			k++;
		}
	//!SORT param1 array++++++++
	vectorSort2(inData.param1,sparam1);
	//+++++++++++++++++++++++++++

	inData.timeDegreOfFreedom=nFiles;
	//	std::cout << "\n";
	it=directory_iterator(p);
	//count number of rows and columns in result files
	//--------------------------------------------
	std::ifstream myfile;
	string new_filename = it->path().generic_string() + string("/")+filenameString+string("0.txt");
	myfile.open(new_filename);
	inData.spaceDegreOfFreedom=std::count(std::istreambuf_iterator<char>(myfile),
			std::istreambuf_iterator<char>(), '\n');
	string row1;
	myfile.clear();
	myfile.seekg(0, std::ios::beg);
	std::getline(myfile,row1);
	std::stringstream stream(row1);
	int cols=std::distance(std::istream_iterator<std::string>(stream), std::istream_iterator<std::string>());
	//--------------------------------------------

	//Fill the matrix array
	//--------------------------------------------
	//inData.A=new Matrix[inData.param1DegreOfFreedom];
	inData.init();
	k=0;
	int colNeeded=4;
	//	while (it != directory_iterator{})
	//	{
	for (int k=0;k<inData.param1DegreOfFreedom;++k)
	{
		//		path file=it->path();
		path file(folderBaseName+sparam1[k]);
		Matrix m1p(inData.spaceDegreOfFreedom,inData.timeDegreOfFreedom);
		for (int t =0;t<inData.timeDegreOfFreedom;++t)
		{
			Matrix m1t(inData.spaceDegreOfFreedom,cols);
			new_filename = file.generic_string() + string("/")+filenameString+std::to_string(t)+string(".txt");
			read(new_filename, m1t);
			auto cl=columns(m1t,{colNeeded});

			auto sm=submatrix( m1p, 0, t, m1p.rows(),1 );
			sm = cl;

		}
		inData.A[k]=m1p;

		//		it++;
		//		k++;
	}
	//--------------------------------------------


}
//---------------------------------------------------------------------------------------------------

void FileManager::readFolder(string foldername)
{
	path p = path(foldername);
	directory_iterator it(p);
	string filenameString="heatTransfer_Conductivity_";
	int nFiles=0;
	inData.param1DegreOfFreedom=0;
	inData.timeDegreOfFreedom=0;

	inData.param1DegreOfFreedom=std::count_if(directory_iterator(p), directory_iterator(), static_cast<bool(*)(const path&)>(is_directory));
	inData.param1.resize(inData.param1DegreOfFreedom);
	SVector sparam1(inData.param1DegreOfFreedom);
	int k=0;

	string folderBaseName;
	while (it != directory_iterator{})
		if(is_directory(it->path()))
		{
			string folder=it->path().generic_string();
			std::size_t	found=folder.rfind("_");
			folderBaseName=folder.substr(0,found+1);
			string folderParam1=folder.substr(found+1);
			inData.param1[k]=my_stod(folderParam1);
//			std::cout<<k <<"/"<< inData.param1[k]<<"\n";
			sparam1[k]=folderParam1;

			it++;
			k++;
		}
	//!SORT param1 array++++++++
	vectorSort2(inData.param1,sparam1);
	//+++++++++++++++++++++++++++

//	inData.timeDegreOfFreedom=nFiles;
	//	std::cout << "\n";
	inData.init();
	k=0;
	for (int k=0;k<inData.param1DegreOfFreedom;++k)
	{
		//		path file=it->path();
		path file(folderBaseName+sparam1[k]);
		Matrix m1p;//(inData.spaceDegreOfFreedom,inData.timeDegreOfFreedom);

		string new_filename = file.generic_string() + string("/")+filenameString+sparam1[k]+string(".odb");
//		std::cout<<sparam1[k]<<": "<<new_filename<<std::endl;
		readODB(new_filename, m1p);
		inData.A[k]=m1p;

	}
	//--------------------------------------------


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

void FileManager::saveModel(string dirPath, const Model& model)
{

	std::ofstream outFile;
	outFile.open((dirPath+"model.dat").c_str(), std::ios_base::out);
	outFile<<model.nbrModes<<" "<<model.F1.rows()<<" "<<model.F2.rows()<<" "<<model.F3.rows()<<"\n";

	for (int i=0;i<model.param1.size();++i)
	{
		outFile<<model.param1[i]<<" ";
	}
	outFile.close();

	write(dirPath+"F1.dat", model.F1);
	write(dirPath+"F2.dat", model.F2);
	write(dirPath+"F3.dat", model.F3);
}
//---------------------------------------------------------------------------------------------------

bool FileManager::loadModel(string dirPath, Model& model)
{

	std::ifstream inFile;
	inFile.open((dirPath+"model.dat").c_str(), std::ios_base::in);
	if (!inFile) {
	    std::cerr << "Unable to open file "+(dirPath+"model.dat");
	    return false;   // call system to stop
	}
	int F1rows,F2rows,F3rows;
	inFile>>model.nbrModes;
	inFile>>F1rows;
	inFile>>F2rows;
	inFile>>F3rows;
	printf("%d %d %d %d \n",model.nbrModes, F1rows, F2rows, F3rows);


	model.param1.resize(F3rows);
	printf("%d\n",model.param1.size());
	std::cout<<sizeof(model.param1[0])<<" ";
	for (int i=0;i<F3rows;++i)
	{
		inFile >> model.param1[i];
		std::cout<<model.param1[i]<<" ";
	}
	printf("\n");
	inFile.close();

	model.F1.resize(F1rows,model.nbrModes);
	model.F2.resize(F2rows,model.nbrModes);
	model.F3.resize(F3rows,model.nbrModes);
	read(dirPath+"F1.dat", model.F1);
	read(dirPath+"F2.dat", model.F2);
	read(dirPath+"F3.dat", model.F3);
	model.cF2=trans(model.F2);
	model.cF3=trans(model.F3);

	return true;
}

using namespace std;
void FileManager::readODB(string filename, Matrix& m)
{
	odb_initializeAPI();

	odb_Odb& myOdb = openOdb(filename.c_str(), true);

	odb_Step& step = myOdb.steps()["Step-1"];
	odb_SequenceFrame& allFramesInStep = step.frames();
	int numFrames = allFramesInStep.size();
	numFrames=c_tmax;
	inData.timeDegreOfFreedom=numFrames;
	int numValues = allFramesInStep[numFrames-1].fieldOutputs()["NT11"].values().size();
	inData.spaceDegreOfFreedom=numValues;
	m = Matrix(inData.spaceDegreOfFreedom,inData.timeDegreOfFreedom);

	for(int timeID =0;timeID<numFrames;++timeID)
	{
		const odb_SequenceFieldValue& temp =
				allFramesInStep[timeID].fieldOutputs()["NT11"].values();
		int numComp = 0;
		for (int spaceID=0; spaceID<numValues; spaceID++)
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
	odb_finalizeAPI();
}
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------

//bool FileManager::loadModel(string dirPath, Model& model)
//{
//
//	FILE * pFile;
//	pFile = fopen((dirPath+"model.dat").c_str(), "r");
//	if (pFile == NULL) return false;
//	int F1rows,F2rows,F3rows;
//
//	fscanf(pFile, "%d %d %d %d",&model.nbrModes, &F1rows, &F2rows, &F3rows);
//	printf("%d %d %d %d \n",model.nbrModes, F1rows, F2rows, F3rows);
//	model.param1.resize(F3rows);
//	printf("%d\n",model.param1.size());
//	std::cout<<sizeof(model.param1[0])<<" ";
//	for (int i=0;i<F3rows;++i)
//	{
//		fscanf(pFile,"%lf", &model.param1[i]);
//		std::cout<<model.param1[i]<<" ";
//	}
//	printf("\n");
//	fclose(pFile);
//	model.F1.resize(F1rows,model.nbrModes);
//	model.F2.resize(F2rows,model.nbrModes);
//	model.F3.resize(F3rows,model.nbrModes);
//	read(dirPath+"F1.dat", model.F1);
//	read(dirPath+"F2.dat", model.F2);
//	read(dirPath+"F3.dat", model.F3);
//	model.cF2=trans(model.F2);
//	model.cF3=trans(model.F3);
//	return true;
//}//---------------------------------------------------------------------------------------------------
//void FileManager::read(string filename, Matrix& m)
//{
//	FILE * pFile;
//	pFile = fopen(filename.c_str(), "r");
//	double f;
//	int xsize=m.columns();
//	int ysize=m.rows();
//	int k=0;
//	if (pFile == NULL) perror ("Error opening file");
//
//	else
//	{
//		while ( ! feof (pFile) && k!=xsize*ysize )
//		{
//
//			fscanf (pFile, "%lf", &f);
////			fscanf (pFile, "%f", &f);
////			std::cout<<f<<" ";
//
//			m(k/xsize,k-(k/xsize)*xsize)=f;
//			++k;
//
//		}
//		fclose (pFile);
//	}
//
//}
////---------------------------------------------------------------------------------------------------
//void FileManager::write(string filename, const Matrix& m)
//{
//	FILE * pFile;
//	pFile = fopen(filename.c_str(), "w");
//	int xsize=m.columns();
//	int ysize=m.rows();
//
//	if (pFile == NULL) perror ("Error opening file");
//	else
//	{
////		std::cout<<"xsize="<<xsize<<" ysize="<<ysize<<"\n";
//		for (int i=0;i<ysize;++i)
//		{
//			for (int j=0;j<xsize;++j)
//			{
//
//				fprintf (pFile, "%lf ", m(i,j));
////				std::cout<<m(i,j)<<" ";
//
//			}
//			fprintf (pFile, "\n");
////			std::cout<<std::endl;
//		}
//	}
//	fclose (pFile);
//}
//
//FileManager::~FileManager()
//{
//	//delete[] inData.A;
//}
