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



class FileManager
{
public:

	~FileManager();
	void read(string filename, Matrix& m);

	//Open snapshot dir,
	//read @param1DegreOfFreedom@, @timeDegreOfFreedom@, @timeDegreOfFreedom@
	void readFolder(string foldername);
	void readFolder(string filename, myMatrix& m);
	void readODB(string filename, Matrix& m);


	void write(string filename, const Matrix& m);
	void saveModel(string dirPath, const Model& model);
	bool loadModel(string dirPath, Model& model);

	inputData inData;
//  Old read from txt files
//	void readFolderOld(string foldername);
};


#endif /* FILEMANAGER_H_ */
