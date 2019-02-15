/*
 * fileManager.h
 *
 *  Created on: 6 Nov 2018
 *      Author: ikriuchevs
 */

#ifndef FILEMANAGER_H_
#define FILEMANAGER_H_

#include <string>
#include "NDMatrix.h"
//typedef SimpleMatrix2D Matrix;

using std::string;

const string c_parameters_filename="params.dat";
const string c_all_parameters_filename="params.in";
const string c_odb_filename="just.odb";


class FileManager
{
public:
	FileManager();
	~FileManager();
	void read(string filename, Matrix& m);

	//Open snapshot dir,
	//read @param1DegreOfFreedom@, @timeDegreOfFreedom@, @timeDegreOfFreedom@
	void readFolder(string foldername, NinputData3& Ndata);
	void readFolder(string foldername, NinputData3& Ndata, int dim, intVector sizes);
	void readParams(string foldername, NinputData3& Ndata, int dim, intVector sizes);
	void readParams(string foldername, int& dim, intVector& sizes, std::vector<Vector>& params);
	void readODB_SpacexTime(string filename, int& spaceDegreOfFreedom, int& timeDegreOfFreedom);
	void readODB2(string filename, Matrix& m);

	void write(string filename, const Matrix& m);
	void saveModel(string dirPath, const Model& model);
	void saveModel(string dirPath, const NModel& model);

	bool loadModel(string dirPath, Model& model);
	bool loadModel(string dirPath, NModel& model);

	;
//  Old read from txt files
//	void readFolderOld(string foldername);
};


#endif /* FILEMANAGER_H_ */
