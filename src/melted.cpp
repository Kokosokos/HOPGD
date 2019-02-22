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
#include <omp.h>



using namespace boost::filesystem;

using namespace std;
#include "boost/program_options.hpp"
namespace po = boost::program_options;


#include <sys/stat.h> //to check if folder exists
struct stat info;

int main(int argc, char * argv[])
{
	// process arguments
	// Declare the supported options.
	po::options_description desc("Allowed options");
	bool createModelFlag = false;
	bool validationFlag  = false;
	int nmbrModes=20;				//maximum number of modes used
	int nthreads = 8;				//number of threads for openmp blaze-lib parallelization

	string validationFolder;

	desc.add_options()
		    						("help", "produce help message")
		    						("create_model,c", po::bool_switch(&createModelFlag)->default_value(false), "creates model and stores it in \'model\' folder")
		    						("validation,v",po::bool_switch(&validationFlag)->default_value(false), "runs the model validation routine. Validation folder name must be provided with --vfolder")
		    						("vfolder",po::value<string>(&validationFolder), "The folder name of validation snapshots.")
		    						("nmax",po::value<int>(&nmbrModes)->default_value(20), "The maximum number of modes in the model (default = 20).")
		    						("nthreads",po::value<int>(&nthreads)->default_value(8), "The number of cores used by openmp (default = 8).")
		    						;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help") || argc == 1 )
	{
		cout << desc << "\n";
		return 1;
	}
	if (validationFlag && !vm.count("vfolder"))
	{
		cout << "Please provide validation folder (-vfolder)"<< ".\n";
		return 1;
	}

	blaze::setNumThreads( nthreads );
	FileManager file_manager;		//input files reader class
	CreateModes modesCreator3;//(nData);
	string modelFolder =			/*folder to store the model files*/
			"model/";
	//			"/home/ikriuchevs/workspace/melted/tests/Ndim//Jobs_bash_test/model/";
	if(createModelFlag)
	{
		NinputData3 nData;				//input data container

		string snapshotsFolder = 		/*snapshots folder*/
				"./";
		//			"/home/ikriuchevs/workspace/melted/tests/Ndim//Jobs_bash_test/";


		std::cout<<"Reading Snapshots: \n";
		double start = omp_get_wtime();
		file_manager.readFolder(snapshotsFolder, nData );
		std::cout<< omp_get_wtime() - start <<" sec"<<"\n";

		nData.error 	= 0.001;		//Model precision
		nData.nModesMax = nmbrModes;			//Maximum number of modes

		std::cout<<"Creating model: ";
		start = omp_get_wtime();

		modesCreator3.create(nData);
		std::cout<< omp_get_wtime() - start <<" sec"<<"\n";

		std::cout<<"Saving the model: ";
		start = omp_get_wtime();
		file_manager.saveModel(modelFolder, modesCreator3.nmodel);
		std::cout<< "\nDone\n"<<omp_get_wtime() - start <<" sec"<<"\n";
		//
	}
	else
	{
		file_manager.loadModel(modelFolder, modesCreator3.nmodel);
	}

	if(nmbrModes > modesCreator3.nmodel.nbrModes)
	{
		cout<<"Warning: Too large nmax. nmax is set to a number of modes in the model (" <<modesCreator3.nmodel.nbrModes<<")\n";
		nmbrModes = modesCreator3.nmodel.nbrModes;
	}
	if (validationFlag)
	{
		int dim=modesCreator3.nmodel.dim;
//		string validationFolder="/home/ikriuchevs/workspace/melted/tests/Ndim/Validation/";
		Matrix M(modesCreator3.nmodel.F[0].rows(),modesCreator3.nmodel.F[1].rows());
		Matrix Mfit(modesCreator3.nmodel.F[0].rows(),modesCreator3.nmodel.F[1].rows());

		double start = omp_get_wtime();
		double fin;

		FILE * outFile;
		string outFileName=validationFolder+"/nmax"+std::to_string(nmbrModes)+".omp"+std::to_string(nthreads)+".dat";
		outFile = fopen(outFileName.c_str(), "w");
		fprintf (outFile, "#nbr of modes = %d\n",modesCreator3.nmodel.nbrModes);
		fprintf (outFile, "#Folder Error Time CudaError CudaTime\n");

		blaze::setNumThreads( nthreads );

		//CUDA
		modesCreator3.cudaInitND(nmbrModes);
		cv::Mat res;
		res.create(modesCreator3.nmodel.F[0].rows(),modesCreator3.nmodel.F[1].rows(),opencvMatrixType);

		int i=1;
//		for(int i=1;i<=4;++i)
		while(stat((validationFolder+"/"+std::to_string(i)).c_str(), &info) == 0)
		{
			string filename=validationFolder+"/"+std::to_string(i)+"/"+c_odb_filename;//<--odb_filename defined in
			cout<<filename<<endl;
			file_manager.readODB2(filename, M);

			string paramname=validationFolder+"/"+std::to_string(i)+"/params.dat";
			std::ifstream inFile;
			inFile.open(paramname, std::ios_base::in);
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

			inFile.close();
			start = omp_get_wtime();
			cout<<"params: "<<blaze::trans(param);
			if(modesCreator3.fitNewND(param,Mfit, nmbrModes))
			{
				fin=omp_get_wtime()-start;

				double err=fabs(100.0*(1.0-blaze::norm(M)/blaze::norm(Mfit)));

				cout<<"OMP: error = "<<err<<"%; time = "<<fin<<" sec;";
				fprintf (outFile, "%d %f %lf ", i, err, fin );

				start = omp_get_wtime();
				modesCreator3.fitNewNDCuda(param,res,nmbrModes);
				fin = omp_get_wtime()-start;

				opencv2blaze(res,Mfit);

				err=fabs(100.0*(1.0-blaze::norm(M)/blaze::norm(Mfit)));
				std::cout<<" CUDA: error = "<<err<<"%; time = "<<fin  << " sec"<<endl<<endl;
				fprintf (outFile, "%f %lf\n", err, fin );
			}
			else
			{
				fprintf (outFile, "%d Enrichment is needed\n\n",i );
			}
			i++;
		}
		fclose (outFile);

	}
	// VALIDATION


	return(0);
}




