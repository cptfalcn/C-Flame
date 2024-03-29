/*
 * ===========================================================================================
 * This is serial implementation.
 * ===========================================================================================
 */
//This is to be depreciated in the next June release.
#include <stdio.h>
#include <math.h>
#include "Epic.h"
#include <cvode/cvode.h>
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix  */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver  */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */
#include <sunlinsol/sunlinsol_spgmr.h>  /* SPGMR SUNLinearSolver                */
#include <sunlinsol/sunlinsol_spbcgs.h>       /* access to SPBCGS SUNLinearSolver            */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver         */
#include <nvector/nvector_serial.h>
//#include <sundials/sundials_nvector.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TChem_CommandLineParser.hpp"
#include "Epi3SC_KIOPS.h"
#include "EpiP2_KIOPS.h"


#include <chrono>
#include "InitialConditions.h"
#include "Chemistry.h"
#include "Print.h"
#include "CreateIntegrators.h"//This includes the EPI3_KIOPS.cpp/h definitions

//#define TCHEMPB TChem::Impl::IgnitionZeroD_Problem<TChem::KineticModelConstData<Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> >>
//#define TCHEMPB TChem::Impl::IgnitionZeroD_Problem	<TChem::KineticModelConstData	<Kokkos::Device	<Kokkos::Serial, Kokkos::HostSpace>  >	>
//Change to Serial, this can be changed by altering the TChem master build profile to include OPENMP on or off
#define BAR "===================="


//=====================
//Prototypes
//=====================
//Banner
void PrintBanner();

//Misc functions
int CheckStep(realtype, realtype);
int SetInteriorSize(realtype);
int TrackProgress(realtype, realtype, realtype, int);
void ReadData(N_Vector, string);
void PrintProfilingToFile(ofstream &, IntegratorStats* , int , string, string, void*);
void PrintToFile(ofstream &, realtype *, int, realtype, realtype, realtype, realtype, realtype, realtype);
//int CVODEMonitorFunction(N_Vector State);
void postProcess(N_Vector);
//=====================
//Namespaces and globals
//=====================
using namespace std;
//====================
//Main
//====================
int main(int argc, char* argv[])
{
	PrintBanner();//Print One-D, how fancy
	//====================
	// Intial Declarations
	//====================
	static realtype FinalTime	= 0;		//
	static const int NumBands 	= 3;		//Epic stuff, default is 3.
	realtype StepSize			= 0;
	realtype KTime 				= 0;
	int ProgressDots			= 0; 		//From 0 to 100; only care about percentage
	int StepCount 				= 0;
	realtype PercentDone		= 0;
	realtype TNow				= 0;
	realtype TNext				= 0;
	realtype TNextC				= 0;		//CVODE TNext, has a tendency to change TNext
	string MyFile				= "Default.txt";
	string Method				= "EPI3V";	//This is the default method
	static realtype KrylovTol	= 1e-14;
	int SampleNum				= 1;
	int Experiment				= 1; 		//default is hydrogen, set to 0 for kapila
	int Profiling				= 0;		//default to no profiling, need to edit profiling
	int number_of_equations		= 0;		//Varies on problem
	realtype absTol 			= 1e-8;
	realtype relTol 			= 1e-8;
	int startingBasisSizes[]	= {10, 10};	//{3,3} //2855.4227227129413 Ten is the normal number
	int Movie					= 0;
	realtype maxSS				= 1e-3;		//Default value.
	int UseJac					= 1;
	//=====================================================
	//Kokkos Input Parser
	//=====================================================
	std::string chemFile("chem.inp");
 	std::string thermFile("therm.dat");
 	/// parse command line arguments --chemfile=user-chemfile.inp --thermfile=user-thermfile.dat
	TChem::CommandLineParser opts("This example computes reaction rates with a given state vector");
 	opts.set_option<std::string>("chemfile", "Chem file name e.g., chem.inp", &chemFile);
 	opts.set_option<std::string>("thermfile", "Therm file name e.g., therm.dat", &thermFile);
	opts.set_option<realtype>("StepSize", "StepSize desired", &StepSize);//New
	opts.set_option<realtype>("FinalTime", "The final simulation time", &FinalTime);//New
	opts.set_option<std::string>("MyFile","Where we output data", &MyFile);//New
	opts.set_option<realtype>("KrylovTol", "KrylovTolerance", &KrylovTol);
	opts.set_option<int>("UseJac","Will we use the Jacobian", &UseJac);
	opts.set_option<int>("SampleNum", "Sample used", &SampleNum);
	opts.set_option<std::string>("Method", "The time integration method", & Method);
	opts.set_option<int>("Experiment", "Experiment chosen", &Experiment);
	opts.set_option<int>("Profiling", "Output profiling data or not", &Profiling);
    opts.set_option<realtype>("relTol", "Solver Relative Tol", &relTol);
    opts.set_option<realtype>("absTol", "Solver Absolulte Tol", &absTol);
	opts.set_option<realtype>("maxSS", "Solver max step size", &maxSS);
	opts.set_option<int>("Movie", "Generate a data dump set at every step", &Movie);

	const bool r_parse = opts.parse(argc, argv);
	if (r_parse)
		return 0; // print help return
	ofstream myfile(MyFile, std::ios_base::app);
	//=====================================
	//Kokkos Sub-declarations
	//=====================================
	Kokkos::initialize(argc, argv);
	// ALL Kokkos variables are reference counted objects.  They are deallocated within this local scope.
	// Kokkos environments - host device type and multi dimensional arrays
	// note: The 2d view uses row major layout.  Most matrix format uses column major layout.
	{
		using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
		using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
	    // construct TChem's kinect model and read reaction mechanism
    	TChem::KineticModelData kmd(chemFile, thermFile);
    	// construct const (read-only) data and move the data to device.
    	// for host device, it is soft copy (pointer assignmnet).
    	auto kmcd = kmd.createConstData<host_device_type>();
   		// use problem object for ignition zero D problem, interface to source term and jacobian
   		using problem_type = TChem::Impl::IgnitionZeroD_Problem<decltype(kmcd)>;
		//=======================================
		//Set number of equations
		//=======================================
		number_of_equations	= problem_type::getNumberOfEquations(kmcd);
		int vecLength 		= number_of_equations;			//For readability
		int num_eqs 		= number_of_equations;			//To reduce size
		int num_pts 		= 1;							//To reduce size
		//==================================================
		//Prep/Set initial conditions, inherited from Zero-D
		//==================================================
		N_Vector y 			= N_VNew_Serial(vecLength);	//intial conditions/Data
		N_VScale(0.0 , y, y);						//Zero out the y vector.
		realtype *data 		= N_VGetArrayPointer(y);	//Set the y data pointer.
		//N_Vector Test		= N_VClone(y);
		//N_Vector Temp		= N_VClone(y);
		//N_Vector TestJac	= N_VNew_Serial(vecLength*vecLength);
		//N_VScale(0.0, Test, Test);
		//Initial Conditions sub-block
		SetIntCons(Experiment, SampleNum, data);		//Set Initial Conditions
		PrintSuperVector(data, Experiment, 1, BAR);
		//===================================================
		//Set TChem
		//===================================================
    	//TChem does not allocate any workspace internally. You create the work space using NVector,
		//std::vector or real_type_1d_view_type (Kokkos view). Here we use kokkos view.
		//Set up an ignition problem and the workspace. Unused in Experiment 0.
		//=========================================================================
		const int problem_workspace_size = problem_type::getWorkSpaceSize(kmcd);
		real_type_1d_view_type work("workspace", problem_workspace_size);
		myPb2 problem2(num_eqs, work, kmcd, num_pts, y, 0);	//Construct new Problem.
		problem2.SetAdvDiffReacPow(0, 0, 0, 1, 0);//Additional set up.
		problem2.kmd			= kmd;
		void *UserData 			= &problem2;
		int MaxKrylovIters 		= 1000; //max(vecLength, 500);//500 is the base

		//Test Data
		// CHEM_RHS_TCHEM_V2(0, y, Test, UserData);
		// N_VPrint_Serial(Test);
		// CHEM_COMP_JAC_V2(y, UserData);
		// std :: cout << endl;
		// CHEM_JTV_V2(y, Test, 0, y, y, UserData, Temp);
		// //MatrixVectorProduct(number_of_equations, N_VGetArrayPointer(problem2.Jac), y, Temp, N_VGetArrayPointer(Test));
		// N_VPrint_Serial(Test);

		//==================
		//Create integrators
		//==================
		//CVODE
        void * cvode_mem;
		SUNMatrix A					= SUNDenseMatrix(vecLength, vecLength);
		SUNLinearSolver SUPERLS 	= SUNLinSol_SPGMR(y, PREC_NONE, 20);
		SUNNonlinearSolver SUPERNLS = SUNNonlinSol_Newton(y);
		

		//Set EPIxx_KIOPS methods
		IntegratorStats *integratorStats 	= NULL;
		Epi3VChem_KIOPS * EPI3V = new Epi3VChem_KIOPS(CHEM_RHS_TCHEM_V2, CHEM_JTV_V2, 
					SUPER_CHEM_JAC_TCHEM, UserData, MaxKrylovIters, y, vecLength);

		cvode_mem= 	CreateCVODE(CHEM_RHS_TCHEM_V2, CHEM_JTV_V2, CHEM_COMP_JAC_CVODE_V2, UserData, A,
				SUPERLS, SUPERNLS, vecLength, y, relTol, absTol, StepSize, UseJac);
		CVodeSetMaxStep			(cvode_mem, maxSS);
		CVodeSetMaxNumSteps		(cvode_mem, 1e6);
		CVodeSetMaxOrd			(cvode_mem, 2);

		//EPIP2 method, with post processing
		EpiP2_KIOPS * EPIP2 = new EpiP2_KIOPS(CHEM_RHS_TCHEM_V2, CHEM_JTV_V2, UserData, MaxKrylovIters,
				y, vecLength, 3, 0, &postProcess);
		//======================
		//Preparation of rites
		//======================
		if(Movie)
		{
			problem2.Movie	= 1;
			cout << "Enter Jac file location: \n";
			cin >> problem2.dumpJacFile;
		}
		//======================
		//Endure the ritual herein
		//======================
		auto Start=std::chrono::high_resolution_clock::now();//Time integrator
		if(Method == "EPI3V")
        {
			cout << BAR <<"\tEPI3V Selected\t\t" << BAR <<endl;
			integratorStats = EPI3V->Integrate(StepSize, maxSS, absTol, relTol,0.0, FinalTime, NumBands, startingBasisSizes, y);
			problem2.t= FinalTime;
        }
		else if(Method == "CVODEKry")
		{
			cout << BAR <<"\tCVODE Selected\t\t" << BAR <<endl;
			CVode(cvode_mem, FinalTime, y, &TNextC, CV_NORMAL);
			problem2.t = TNextC;
		}
		else if(Method == "EPIP2")
		{
			cout << BAR <<"\tEPIP2 Selected\t\t" << BAR <<endl;
			integratorStats= EPIP2->Integrate(StepSize, 0.0, FinalTime, y, KrylovTol, startingBasisSizes);
			problem2.t = FinalTime;
		}
		auto Stop=std::chrono::high_resolution_clock::now();
		auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
		KTime+=Pass.count()/1e9;
		//=======================================
		//End rites
		//Enscribe upon thyn tablet the results
		//=======================================
		cout << BAR <<"Printing data to "<< MyFile << BAR << endl;
		PrintExpParam(FinalTime, problem2.t, StepSize, StepCount, KrylovTol, absTol, relTol, KTime, BAR);
		cout << setprecision(17);
		PrintSuperVector(data, Experiment, 1, BAR);
		PrintProfiling(integratorStats, Profiling, Method,  BAR, cvode_mem);
		PrintToFile(myfile, data, vecLength, problem2.t, relTol, absTol, KTime, KrylovTol, problem2.ignTime);
		myfile.close();
		cout << "\nTotal integration time: " << KTime << " sec " << endl;
                cout << "\t Chem rhs: " << problem2.rhs_Chem << " sec\n";
                cout << "\t Chem jtv: " << problem2.jtv_Chem << " sec\n";
                cout << "\t jac make: " << problem2.jacMakeTime << " sec\n";
		cout << "Max Step Size: " << problem2.MaxStepTaken << endl;
		cout << "Min Step Size: " << problem2.MinStepTaken << endl;
		//==========================
		//Time to take out the trash
		//==========================
		delete EPI3V;
		N_VDestroy_Serial(y);
		SUNMatDestroy(A);
		SUNLinSolFree(SUPERLS);
		SUNNonlinSolFree(SUPERNLS);
        CVodeFree(&cvode_mem);
		delete EPIP2;
  	}//end local kokkos scope.
  	Kokkos::finalize();	/// Kokkos finalize checks any memory leak that are not properly deallocated.
	cout << BAR << "\tExiting without error\t" << BAR <<endl;
	return 0;
}


/*
//= ___=====================================
// |  _) _ _  ___  __  _-_  _  ___  ___  __
// |  _)| | ||   |/ _)(. .)[_]/ _ \|   |( .)
// |_|  |___||_|_|\__) |_| |_|\___/|_|_|(`_)
//===========================================
*/
//===========================================
//Check the number of steps
//===========================================
int CheckStep(realtype FinalTime, realtype StepSize)
{
	int Steps=0;
	//cout << BAR << "\t Step Check\t\t" << BAR << endl ;
	//cout << std::setprecision(17) << FinalTime/StepSize << " steps proposed...";
	cout << std::setprecision(17) << endl;
        if(floor( FinalTime/StepSize ) == FinalTime/StepSize )
	{
                Steps= FinalTime/StepSize;
                //cout << " accepted...\n";
		return Steps;
        }else if (abs(round(FinalTime/StepSize))-FinalTime/StepSize <1e-6 ){
		Steps= round (FinalTime/StepSize);
		//cout << Steps << " steps approximated...\n";
		return Steps;
	}else{
                cout<<"Cannot perform non-integer number of steps!!\n";
                exit(EXIT_FAILURE);
        }
}


void PrintBanner()
{
	cout << "\n\n\n";
	cout << "\t================================================\n";
	cout << "\t||      <:::>   :>   ::  <::::>      <::::>    ||\n";
	cout << "\t||     <>   <>  |:>  ::  |::         |::  `>   ||\n";
	cout << "\t||    <>     <> |: > ::  |::::> <::> |::   *>  ||\n";
	cout << "\t||     <>   <>  |:  >::  |::         |::  .>   ||\n";
	cout << "\t||      <:::>   |:   >:  <::::>      <::::>    ||\n";
	cout << "\t================================================\n";

}

void PrintToFile(ofstream & myfile, realtype * data, int length, realtype t, realtype rel_tol, realtype abs_tol,
			realtype run_time, realtype kry_tol, realtype ignTime)
{
	for (int i=0; i < length; i++)
        {
                myfile<<setprecision(20)<<fixed <<data[i] <<"\t\t";
        }
	myfile << "\t\t" << t << "\t\t" << run_time << "\t\t" << rel_tol << "\t\t" << abs_tol;
	myfile << "\t\t" << kry_tol << "\t\t" << ignTime << endl;

}
//========================
//Creates the progress bar
//========================
int TrackProgress(realtype FinalTime, realtype TNext, realtype PercentDone, int ProgressDots)
{
	PercentDone=floor(TNext/FinalTime*100);
        for (int k=0; k<PercentDone-ProgressDots;k++)
   	{
      		cout<<":";
		cout.flush();
	}
	return PercentDone;
}


void PrintProfilingToFile(ofstream & myfile, IntegratorStats* integratorStats, int Profiling, string Method,
				string Bar, void* cvode_mem)
{
	if(Profiling == 1 && Method != "CVODEKry")
	{
		myfile << integratorStats->numTimeSteps;//>krylovStats[1].numIterations;
	}
	else if(Profiling == 1 && Method == "CVODEKry")
	{
		long int steps          = 0;
                long int fEval          = 0;
                long int linSetups      = 0;
                long int errTestFails   = 0;
                int lastOrd             = 0;
                int nextOrd             = 0;
                realtype realHInit      = 0;
                realtype hLast          = 0;
                realtype hCurr          = 0;
                realtype tCurr          = 0;
                long int nonLinIters    = 0;
                long int Projections    = 0;
                CVodeGetNumSteps(cvode_mem, &steps);
                CVodeGetIntegratorStats(cvode_mem, &steps, &fEval, &linSetups, &errTestFails, &lastOrd,
                                        &nextOrd, &realHInit, &hLast, &hCurr, &tCurr);
                CVodeGetNumNonlinSolvIters(cvode_mem, &nonLinIters);
                myfile << steps << "\t\t";
	}

}

void postProcess(N_Vector solution)
{
	realtype * data = N_VGetArrayPointer(solution);
	int len 		= N_VGetLength(solution);
	for (int i = 0; i < len; i ++ )
	{
		//data[i] = std::max(1e-10, data[i]);
	}
}
