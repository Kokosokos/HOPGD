/*
 * fileManager.h
 *
 *  Created on: 6 Nov 2018
 *      Author: ikriuchevs
 */

#ifndef FILEMANAGER_H_
#define FILEMANAGER_H_

#include <string>
#include "mymatrix.h"

//typedef SimpleMatrix2D Matrix;

using std::string;

//Add copy operator!!
struct inputData
{
	//number of result files:
	int   timeDegreOfFreedom;
	//number of rows in each result file
	int  spaceDegreOfFreedom;
	//number of different Temperature(or another param) points = #of folders
	int param1DegreOfFreedom;
	float error;
	Vector param1;

	Matrix* A;
	inputData(){};
	void init(){A=new Matrix[param1DegreOfFreedom];};
	~inputData(){ delete[] A;};
};

struct Model
{
	int nbrModes;
	Matrix F1;
	Matrix F2;
	Matrix F3;
	cMatrix cF2;
	cMatrix cF3;
	Vector param1;
};

class FileManager
{
public:
	//Per element read
	//must be a better way
	~FileManager();
	void read(string filename, Matrix& m);

	//Open snapshot dir,
	//read @param1DegreOfFreedom@, @timeDegreOfFreedom@, @timeDegreOfFreedom@
	void readFolder(string foldername);
	void readFolderOld(string foldername);
	void readODB(string filename, Matrix& m);

	void write(string filename, const Matrix& m);
	void saveModel(string dirPath, const Model& model);
	bool loadModel(string dirPath, Model& model);

	inputData inData;
};


#endif /* FILEMANAGER_H_ */
