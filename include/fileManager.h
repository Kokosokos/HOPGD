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

using std::string;
//The 1 snapshot parameters file name.
const string c_parameters_filename="params.dat";
//The parameters file name for multiple snapshots.
const string c_all_parameters_filename=""
		"params.in";
//The .odb file name.
const string c_odb_filename="just.odb";

/**
 * @brief The class that manages file input/output routines.
 */
class FileManager
{
public:
	/**
	 * @brief Initialises odbApi.
	 */
	FileManager();
	/**
	 * @brief Finalises odbApi.
	 */
	~FileManager();

	/**
	 * @brief Snapshots reading routine.
	 * @param[in] foldername The name of the folder to read the snapshots from.
	 * The snapshots must be prepared with "prepare_folders.sh" script and "params.in" file must be in the same folder.
	 * @param[out] Ndata The output NinputData3 structure.
	 */
	void readFolder(string foldername, NinputData3& Ndata);
	void readFolder(string foldername, NinputData3& Ndata, int dim, intVector sizes);
	void readParams(string foldername, NinputData3& Ndata, int dim, intVector sizes);
	//old: reads parameters from each folder instead of "params.in"
	void readParams(string foldername, int& dim, intVector& sizes, std::vector<Vector>& params);
	/**
	 * @brief Reads the size of the space and time dimensions from .odb file.
	 * @param[in] filename The name of the file to read the matrix from.
	 * @param[out] spaceDegreOfFreedom The size of the "space" dimension. (number of space points)
	 * @param[out] timeDegreOfFreedom  The size of the "time" dimension.  (number of time points)
	 */
	void readODB_SpacexTime(string filename, int& spaceDegreOfFreedom, int& timeDegreOfFreedom);
	/**
	 * @brief Reads matrix from .odb file.
	 * @param[in] foldername The name of the folder to read the snapshots from.
	 * The snapshots must be prepared with "prepare_folders.sh" script and "params.in" file must be in the same folder.
	 * @param[out] m The output matrix.
	 */
	void readODB2(string filename, Matrix& m);
	/**
	 * @brief Reads matrix from text file.
	 * @param[in] filename The filename to read matrix from.
	 * @param[out] m The output matrix.
	*/
	void read(string filename, Matrix& m);
	/**
	 * @brief Writes the matrix to text file.
	 * @param[in] filename The filename to write matrix to.
	 * @param[in] m The input matrix.
	*/
	void write(string filename, const Matrix& m);

	/**
	 * @brief Saves the ND model to text file (see NModel).
	 *
	 * Creates "model.dat" and "F0.dat", .. "F(N-1).dat" files.<br>
	 *<pre>
	 * File structure of "model.dat":
	 *
	 *  		Number of modes
	 *  		Dimensionality (N)
	 *  		Sizes of each dimension
	 *  		first parameters set
	 *  		...
	 *  		N-th parameters set
	 * </pre>
	 * @param[in] model The ND model.
	 * @param[in] dirPath The directory path to store the model files to.
	 */
	void saveModel(string dirPath, const NModel& model);

	/**
	 * @brief Saves the ND model to text file (see NModel).
	 *
	 * Creates model.dat and {F0.dat, .. F(N-1).dat} files.<br>
	 * See saveModel() for model.dat file structure.
	 */
	bool loadModel(string dirPath, NModel& model);

};


#endif /* FILEMANAGER_H_ */
