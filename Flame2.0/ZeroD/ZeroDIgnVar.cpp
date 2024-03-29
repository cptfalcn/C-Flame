/*
 * ===========================================================================================
 * This is allows for openmp in kokkos
 * ===========================================================================================
 */
#include <stdio.h>
#include <math.h>
#include "Epic.h"
#include <nvector/nvector_serial.h>
#include <nvector/nvector_openmp.h>
#include <cvode/cvode.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix  */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver  */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */
#include <sunlinsol/sunlinsol_spgmr.h>  /* SPGMR SUNLinearSolver                */
#include <sunlinsol/sunlinsol_spbcgs.h>       /* access to SPBCGS SUNLinearSolver            */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver         */
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TChem_CommandLineParser.hpp"
#include "TChem_Util.hpp"
#include "TChem_KineticModelData.hpp"
#include "TChem_Impl_IgnitionZeroD_Problem.hpp" // here is where Ignition Zero D problem is implemented
#include <chrono>
#include "InitialConditions.h"
#include "Epi3V.h"
#include "Epi3VCntrl.h"
#include "Print.h"
//#include "Leja.hpp"

using value_type = realtype;
#define TCHEMPB TChem::Impl::IgnitionZeroD_Problem<value_type, Tines::UseThisDevice<TChem::host_exec_space>::type >
//Change to Serial/OMP by altering the TChem master build profile to include OPENMP on or off
#define BAR "===================="

using namespace std;
//=====================
//Prototypes & Classes
//=====================
class myPb : public TCHEMPB{
	public:
	//members
	ordinal_type  			num_equations;
	N_Vector 				Jac;
	realtype 				t;
	realtype 				MaxStepTaken;
	realtype 				MinStepTaken;
	realtype 				ignTime;
	realtype 				KiopsTime;				
	SUNMatrix				Mat;
	int 					Movie;
	std :: string			dumpJacFile;
	int						InternalSteps;
	int						BadErrSteps;
	int						BlowupSteps;
	int						KiopsBlowups;
	std :: string			stepRatioFile;
	realtype 				stepRatio;
	realtype				ProjectTime;
	realtype				OrthogTime;
	//functions
	ordinal_type			get_num_equations(void)
	{
		return this->num_equations;
	}
};


//Banner
void PrintBanner();

//Initial conditions
//See InitialConditions.h

//RHS functions
int RHS_TCHEM(realtype, N_Vector, N_Vector, void *);

//JtV functions
// int Jtv_TCHEM		(N_Vector , N_Vector ,realtype, N_Vector, N_Vector , void* , N_Vector );
int JtV_TCHEM_Fast	(N_Vector , N_Vector ,realtype, N_Vector, N_Vector , void* , N_Vector );


//Misc functions
void PrintFromPtr(realtype *,  int);
void ErrorCheck(ofstream &, N_Vector, realtype *, int, realtype);
void PrintDataToFile(ofstream &, realtype *,int, realtype);
//Main Print function
void PrintToFile(ofstream &, realtype *, int, realtype, realtype, realtype, realtype, realtype, realtype);



//Used in Jtv
int ComputeJac(N_Vector, void*);
int CVodeComputeJacWrapper(realtype t, N_Vector u, N_Vector fy, SUNMatrix Jac,
               void * pb, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
void MatrixVectorProduct(int number_of_equations, realtype * JacD, N_Vector x, 
				N_Vector tmp, realtype * JV);
void MatrixVectorProductFast(int len, N_Vector Jac, N_Vector x, N_Vector tmp, realtype * JV);

//=====================
//Namespaces and globals
//=====================
using namespace std;

int JacCnt=0;
int RHSCnt=0;
int JtvCnt=0;

double JacTime=0;
double RHSTime=0;
double MatVecTime=0;
Kokkos::Timer g_timer;

//====================
//Main
//====================
int main(int argc, char* argv[])
{
	PrintBanner();					//Print Zero-D, how fancy
	sundials::Context sunctx;   	// SUNDIALS context
	//====================
	// Intial Declarations
	//====================
	static realtype FinalTime 	= 0;//1.0e-4;//1e-4 seems to be the max
	static const int NumBands 	= 3;//Epic stuff, default is 3.
	realtype StepSize			= 0;
	realtype KTime				= 0;
	int ProgressDots			= 0; //From 0 to 100; only care about percentage
	long int StepCount 			= 0;
	realtype TNow				= 0;
	realtype TNext				= 0;
	string MyFile				= "Default.txt";
	string Method				= "EPI2";//This is the default method
	static realtype KrylovTol	= 1e-14;
	int UseJac					= 1; //will we use the Jacobian or not
	int SampleNum				= 0;
	int Experiment				= 1; //default is hydrogen, 2 for Gri3.0
	int Profiling				= 0;//default to no profiling
	int number_of_equations		= 0;
	int startingBasisSizes[] 	= {10, 10};
	realtype relTol				= 1e-10;
	realtype absTol				= 1e-10;
	realtype maxSS 				= StepSize;
	int Movie 					= 0;
	string JacFile				= "";
	string inputFile;
	int PressMult = 1;
	

	//=====================================================
	//Kokkos Input Parser
	//=====================================================
	std::string chemFile("chem.inp");
 	std::string thermFile("therm.dat");
 	/// parse command line arguments --chemfile=user-chemfile.inp --thermfile=user-thermfile.dat
	/// --Stepsize=1e-8  --FinalTime=1e-2  --MyFile="Filename"  --KrylovTole=1e-14
	//  --UseJac=0 --SampleNum={1,2,3}, --Method="Integrator name"
	//  --Experiment={0,1,2}
  	/// with --help, the code list the available options.
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
	opts.set_option<string>("Input", "Input TChem file name", &inputFile);
	opts.set_option<realtype>("maxSS", "Maximum StepSize", &maxSS);
	opts.set_option<int>("PressMult", "Pressure Multiplier", &PressMult);
	opts.set_option<int>("Movie", "Jacobian Movie data", &Movie);
	opts.set_option<std::string>("JacFile", "OptionalJacFile", &JacFile);
	const bool r_parse = opts.parse(argc, argv);

	if (r_parse)
		return 0; // print help return

	ofstream myfile(MyFile, std::ios_base::app);
	//============
	//Start Kokkos
	//============
	Kokkos::initialize(argc, argv);
	{//begin local scope
		//======================
		//TChem Setup
		//======================
		using TChem::real_type;
		using TChem::ordinal_type;
		/// Kokkos environments - host device type and multi dimensional arrays
		using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
		using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
		using real_type_2d_view_type = Tines::value_type_2d_view<real_type,host_device_type>;
		//Kokkos ZeroD problem type
		using problem_type = TChem::Impl::IgnitionZeroD_Problem<value_type, host_device_type >;	
		/// construct TChem's kinect model and read reaction mechanism
		TChem::KineticModelData kmd(chemFile, thermFile);
		auto kmcd = TChem::createGasKineticModelConstData<host_device_type>(kmd);
		/// set workspace and #eqs and output information
		real_type_1d_view_type work("workspace", problem_type::getWorkSpaceSize(kmcd));

		//New stuff
		cout << BAR << " Mech and parallel top " << BAR << endl;
		printf("Number of Species %d \n", kmcd.nSpec);
		printf("Number of Reactions %d \n", kmcd.nReac);
		const bool detail = false;
		TChem::exec_space::print_configuration(std::cout, detail);
		//===================================
		//Declare the state variable locally
		//===================================
		number_of_equations		= problem_type::getNumberOfEquations(kmcd);
		N_Vector y 				= N_VNew_Serial(number_of_equations, sunctx); //state
		N_VScale(0.0, y, y);
		realtype *data 			= NV_DATA_S(y);	//set the pointer to the state
		//const int MaxKrylovIters= max(30,10);//use 1000
		const int MaxKrylovIters= max(number_of_equations,10);//use 1000

		//Print mechanism data
		//only For build tests, comment out otherwise.
		// int nBatch= 1;
		// real_type_2d_view_host state_host;
		// const ordinal_type stateVecDim = TChem::Impl::getStateVectorSize(kmcd.nSpec);
		// const auto speciesNames 			= kmcd.speciesNames;
		// const auto SpeciesMolecularWeights 	= kmcd.sMass;
		// TChem::Test::readSample(inputFile, speciesNames, SpeciesMolecularWeights,
		// 		kmcd.nSpec, stateVecDim, state_host, nBatch);
		//======================
		//Init Cons
		//=======================
		
		if(Experiment == 3)
		{
			PressMult = 10;
		}

		if(Experiment == 2 && SampleNum == 10)
		{
			PressMult = 10;
		}
		
		SetIntCons(Experiment, SampleNum, data);		//Will depreciate?
	    //=================
		/// set problem
		//=================
		myPb problem;
		void *pbptr= &problem;
		const real_type pressure(101325);//Constant pressure
		real_type_1d_view_type fac("fac", 2*number_of_equations);
      	problem._p 		= pressure*PressMult; // pressure
		problem._fac	= fac;
      	problem._work 	= work;  // problem workspace array
     	problem._kmcd 	= kmcd;  // kinetic model
		problem.num_equations		= number_of_equations;
		problem.Jac					= N_VNew_Serial(number_of_equations*number_of_equations, sunctx);//Make the Jacobian
		problem.t					= 0;
		problem.stepRatioFile		= "";
		problem.dumpJacFile			= JacFile;
		problem.Movie				= Movie;
		problem.stepRatioFile		= "ZeroDVariableExp"+to_string(Experiment)+"Stepfile.txt";
		problem.KiopsTime			= 0;
		//==============================================
		//Create integrators
		//==============================================
		void * cvode_mem;
		SUNMatrix A						= SUNDenseMatrix(number_of_equations, number_of_equations,sunctx);
		SUNLinearSolver LS 				= SUNLinSol_SPGMR(y, PREC_NONE, number_of_equations,sunctx);
		SUNNonlinearSolver NLS 			= SUNNonlinSol_Newton(y, sunctx);

		int retVal = 0;
        //realtype tret=0;
        N_Vector AbsTol= N_VNew_Serial(number_of_equations, sunctx);
        for ( int i = 0 ; i < number_of_equations ; i++)
			NV_Ith_S(AbsTol,i)=absTol;
        cvode_mem = CVodeCreate (CV_BDF, sunctx);
        retVal = CVodeSetUserData(cvode_mem, pbptr);
        retVal = CVodeInit(cvode_mem, RHS_TCHEM, 0, y);
        retVal = CVodeSetInitStep(cvode_mem, StepSize);
        retVal = CVodeSetMaxStep(cvode_mem, maxSS);
		retVal = CVodeSetMinStep(cvode_mem, 1e-20);
		retVal = CVodeSetMaxErrTestFails(cvode_mem, 200);
        retVal = CVodeSVtolerances(cvode_mem, relTol, AbsTol);
        //Set linear solver
        retVal = CVodeSetLinearSolver(cvode_mem, LS, A);//A might be null, try that
        retVal = CVodeSetNonlinearSolver(cvode_mem, NLS);
		retVal = CVodeSetJacTimes(cvode_mem, NULL, JtV_TCHEM_Fast);
		retVal = CVodeSetLSetupFrequency(cvode_mem, 1); //Remove because also stupid.
		retVal = CVodeSetJacEvalFrequency(cvode_mem, 1);//Remove becuase this is stupid.
		retVal = CVodeSetJacFn(cvode_mem, CVodeComputeJacWrapper);//Error if removed .
		retVal = CVodeSetMaxNumSteps(cvode_mem, 1e6);
		retVal = CVodeSetMaxOrd(cvode_mem, 3);
		//Set Epi3V
		IntegratorStats *integratorStats = NULL;
		Epi3VChem_KIOPS * EPI3V = new Epi3VChem_KIOPS(RHS_TCHEM, JtV_TCHEM_Fast, 
			CVodeComputeJacWrapper, pbptr, MaxKrylovIters, y, number_of_equations);
		Epi3VCntrl * EPI3VCntrl = new Epi3VCntrl(RHS_TCHEM, JtV_TCHEM_Fast, 
			CVodeComputeJacWrapper, pbptr, MaxKrylovIters, y, number_of_equations);

		//================
		//Pre-run printing
		//================
		cout << "Pressure: " << pressure*PressMult << endl;
		//PrintPreRun(StepSize, 0.0, 0.0, KrylovTol, absTol, relTol, Method, number_of_equations, BAR);
		PrintSuperVector(data, Experiment, 1, BAR);
		//=================================
		//Select and run integrator
		//=================================
		cout << BAR << Method << " selected " << BAR << endl;
		auto Start=std::chrono::high_resolution_clock::now();//Time integrator

		if(Method == "CVODEKry")
			CVode(cvode_mem, FinalTime, y, &TNext, CV_NORMAL);//CV_NORMAL/CV_ONE_STEP
		else if(Method == "EPI3V")
		{
			integratorStats = EPI3V->NewIntegrate(StepSize, maxSS, absTol, relTol, 
								0.0, FinalTime, startingBasisSizes, y);
		}
		else if(Method =="EPI3VCntrl")
		{
			integratorStats = EPI3VCntrl->NewIntegrate(StepSize, maxSS, absTol, relTol, 
								0.0, FinalTime, startingBasisSizes, y);
		}
		else
		{
			cout << "Invalid integrator selected\n";
			exit(EXIT_FAILURE);
		}
		auto Stop=std::chrono::high_resolution_clock::now();
		auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
		KTime+=Pass.count()/1e9;

		cout << BAR << "\tIntegration complete\t" << BAR <<endl;
		//=======================================
		//Console Output
		//=======================================
		//Profiling output
		realtype EffRating 		=	1;
		cout << "Mass Fraction error: "<<abs( N_VL1NormLocal(y)-data[0]-1.0)<<endl;
		if (Profiling ==1){//Invalid for experiment 0
			cout << BAR << "    Profiling   " << BAR << endl;
			ofstream ProFile("ZeroDVariableExp"+to_string(Experiment)+ "Method"+ Method + "Sample" + to_string(SampleNum) + "Pressure" + to_string(PressMult)+"Profile.txt", std::ios_base::app);
			//ofstream ProFile("Profiling.txt", std::ios_base::app);//Profiling  File

			if(Method == "EPI3V" || Method == "EPI3VCntrl")
			{
				integratorStats->PrintStats();
				EffRating 	= 100*integratorStats->numTimeSteps/static_cast<realtype>(problem.InternalSteps);
				cout << BAR << "Overall Efficiency rating" << BAR << endl;
				cout << EffRating <<"%\n";
				StepCount 		= integratorStats->numTimeSteps;
				cout << "Step rejection break-down\n";
				cout << "Poor Error Est: " << problem.BadErrSteps << endl;
				cout << "Bad Error Est : " << problem.BlowupSteps << endl;
				cout << "Kiops Errors  : " << problem.KiopsBlowups<< endl;
				cout << "-----------------------------------------\n";
				cout << "Total Internal Error: " << problem.BadErrSteps + problem.BlowupSteps + problem.KiopsBlowups << endl;
				cout << "Sucessful time steps: " << StepCount << endl;
				cout << "Total Internal steps: " << problem.InternalSteps << endl;

				// ProFile << problem.BadErrSteps << "\t";
				// ProFile << problem.BlowupSteps << "\t";
				// ProFile << problem.KiopsBlowups << "\t";
				ProFile << integratorStats->krylovStats->numIterations << "\t";
				ProFile << integratorStats->krylovStats->numProjections << "\t";
				ProFile << integratorStats->krylovStats->numRejections << "\t";

			}
			else if(Method == "CVODEKry")
			{
				FILE * pFile;
				pFile = fopen ("Profiling.txt" , "a");
				//Want:
				//total steps
				//Non-linear iterations
				long int nliters 	= 0 ;
				long int nniters	= 0 ;
				//Fails
				long int nncfails   = 0 ;
				long int nlcfails	= 0 ;
				//Get total steps
				CVodeGetNumSteps(cvode_mem, &StepCount);

				//Get linear iterats
				CVodeGetNumLinIters(cvode_mem, & nliters);
				//Get nonlinear solve iters and nonlinear fails
				CVodeGetNonlinSolvStats(cvode_mem, & nniters, & nncfails);
				CVodeGetNumLinConvFails(cvode_mem, &nlcfails);
				cout << "Steps: " << StepCount << endl;
				cout << "Linear Iterations: " << nliters << endl;
				cout << "NonLinear Iterations: " << nniters << endl;
				//cout << "Linear Convergence fails: " << nlcfails << endl;
				cout << "Linear Iterations per step: " << nliters/(realtype)StepCount << endl;
				cout << "Linear Iterations per NonLinear iteration: " << nliters/(realtype)nniters<<endl;
				cout << "Nonlinear Iterations per Step: " << nniters/(realtype)StepCount << endl;
				CVodePrintAllStats(cvode_mem, pFile, SUN_OUTPUTFORMAT_TABLE);
				fclose(pFile);
				ProFile << nliters << "\t" << nniters << "\t" << nlcfails << "\t";
			}

			ProFile << StepCount << "\t\t";
			cout << "General \t " << BAR << endl;
			cout << "Steps: " << StepCount << endl;
			cout << BAR << "Jacobian\t" << BAR << endl;
			cout << "Jacobian calls: " << JacCnt << endl;
			cout << "Jacobian time:  " << JacTime << endl;
			ProFile << JacCnt << "\t" << JacTime << "\t\t";

			cout << BAR << "RHS\t\t" << BAR << endl;
			cout << "RHS calls: " << RHSCnt << "\n";
			cout << "RHS time:  " << RHSTime << endl;
			ProFile << RHSCnt << "\t" << RHSTime << "\t\t";

			cout << BAR << "Jtv\t\t" << BAR << endl;
			cout << "Jtv calls: " << JtvCnt << endl;
			cout << "Jtv time:  " << MatVecTime << endl;

			cout << "Kiops Time: " << problem.KiopsTime << endl;
			ProFile << JtvCnt << "\t" << MatVecTime << "\t\t" << KTime << "\t" << problem.KiopsTime <<"\t\t"; 

			
			if (Method == "EPI3V")
			{
				cout << "Kiops internal time" << BAR << endl;
				cout << "Orthog time: " << EPI3V->NewKrylov->OrthogTime << endl;
				cout << "Projection time: " << EPI3V->NewKrylov->ProjectTime << endl;
				ProFile << EPI3V->NewKrylov->OrthogTime << "\t";
				ProFile << EPI3V->NewKrylov->ProjectTime << "\t";
				ProFile << EPI3V->NewKrylov->AdaptTime << "\t";

			}
			else
			{
				ProFile << 0.0 << "\t";
				ProFile << 0.0 << "\t";
				ProFile << 0.0 << "\t";

			}

			ProFile << "\t " << Experiment << "\t" << SampleNum << "\t";
			ProFile << FinalTime << "\t" << StepSize << "\t";
			ProFile << absTol <<"\t" << relTol << endl;
			ProFile.close();


		}//End Profiling
		cout << BAR <<"Printing data to "<< MyFile << BAR << endl;
		PrintExpParam(FinalTime, FinalTime, StepSize, StepCount, KrylovTol, absTol, relTol, KTime, BAR);
		PrintSuperVector(data, Experiment, 1, BAR);
		PrintToFile(myfile, data, number_of_equations, problem.t, relTol, absTol, KTime, StepCount, EffRating);
		//N_VPrint_Serial(y);
		//Take out the trash
        myfile.close();
		N_VDestroy_Serial(y);
		N_VDestroy_Serial(AbsTol);
		SUNMatDestroy(A);
		SUNLinSolFree(LS);
		SUNNonlinSolFree(NLS);
		CVodeFree(&cvode_mem);

  	}//end local kokkos scope.
 	/// Kokkos finalize checks any memory leak that are not properly deallocated.
  	Kokkos::finalize();
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

// ====================================================
//  ___    _   _   ___
// |   \  | | | | /   \
// | |) ) | |_| | \ \\/
// |   <  |  _  |  \ \
// | |\	\ | | | | /\\ \
// |_| \_\|_| |_| \___/
// ====================================================
/*
 * ===========================================================================================
 *
 * Function RHS
 *
 * If y' = f(t, y) then RHS is function f.
 * Function RHS is written is such way that it can handle arrays or arbitrarly length.
 * It is important that all data is in form of an array.
 *
 * Inputs:
 * t          time
 * u          input vector
 * udot       result
 * userData
 *
 * ===========================================================================================
//============================================================
//RHS that uses the TCHEM rhs functions
//============================================================
*/
int RHS_TCHEM(realtype t, N_Vector u, N_Vector udot, void * pb)
{
	//TCHEMPB *pbPtr{ static_cast<TCHEMPB*>(pb)} ;//recast the type here.
	myPb * pbPtr{static_cast<myPb *> (pb)};//Recast
	ordinal_type number_of_equations=pbPtr->get_num_equations();
	realtype *y= NV_DATA_S(u);
	realtype *dy=NV_DATA_S(udot);
	/// we use std vector mimicing users' interface
	using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
	using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
    // output rhs vector and x via a kokkos wrapper.
	real_type_1d_view_type x(y,	number_of_equations);
	real_type_1d_view_type f(dy,	number_of_equations);
	// =================================
	// Compute right hand side vector
	// =================================
	auto member =  Tines::HostSerialTeamMember();
	g_timer.reset();
	pbPtr->computeFunction(member, x ,f);
	RHSTime+=g_timer.seconds();
	RHSCnt++;
	return 0;
}
//===============================================================
//   _____   ___     ____   ___    _____   _____    ___   __   _
//  |__ __| /   \   / ___) / _ \  |     \ |__ __|  /   \ |  \ | |
//  _ | |  | (x) | | /    / / \ \ |  x  /   | |   | (x) ||   \| |
// / (| |  |  n  | | \___ \ \_/ / |  x  \  _| |   |  n  || |\   |
// \____/  |_| |_|  \____) \___/  |_____/ |_____| |_| |_||_| \__|
//===============================================================
int ComputeJac(N_Vector u, void* pb)
{
	//problem_type problem;
	myPb * pbPtr{static_cast<myPb *> (pb)};//Recast
	ordinal_type number_of_equations=pbPtr->num_equations;//Get number of equations
	//======================================
	//Set the necessary vectors and pointers
	//======================================
	realtype *y= NV_DATA_S(u);
	realtype *JacD= NV_DATA_S(pbPtr->Jac);
	//=============================================
	//We use std vector mimicing users' interface
	//=============================================
	using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
	using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
	using real_type_2d_view_type = Tines::value_type_2d_view<real_type,host_device_type>;

	//============================================
	// output rhs vector and J matrix wrappers.
	//============================================

	real_type_1d_view_type x(y   ,   number_of_equations);
	real_type_2d_view_type J(JacD,   number_of_equations, number_of_equations);

	//===============================
	/// Compute Jacobian
	//===============================
	auto member =  Tines::HostSerialTeamMember();
	g_timer.reset();
	pbPtr->computeJacobian(member, x ,J);
	JacTime+=g_timer.seconds();
	JacCnt++;
	return 0;
}

//==============
//Cvode wrapper
//==============
int CVodeComputeJacWrapper(realtype t, N_Vector u, N_Vector fy, SUNMatrix Jac,
               void * pb, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	ComputeJac(u, pb);
	return 0;
}


//==============================//
//	  	  | |    _				//
//	      | |  _| |_  __    __	//
//     _  | | (_   _) \ \  / /	//
//    ( (_| |   | |    \ \/ /	//
//     \___/    |_|     \__/	//
//==============================//
/*
 * ===========================================================================================
 *
 * Function Jtv
 *
 * This function computes Jacobian matrix times some vector v. In Epirk Jacobian is never
 * stored and it is computed every time we need it and every time we compute it we acctualy
 * compute it we compute its product with some vector v.
 * Function Jtv is written is such way that it can handle arrays or arbitrarly length.
 * It is important that all data is in form of an array.
 * It is used in JTimesv class.
 *
 * Inputs:
 * v          vector being multiplied with Jacobian Matrix
 * Jv         result
 * t          time
 * u          vector used ot compute Jacobian Matrix
 * fu         f(u), i.e., RHS of u
 * userData
 * tmp
 * ===========================================================================================
//=============================================================
||Jtv using TCHEM
||N_Vector v: 	the v in Jv
||N_Vector Jv:  Jv
||realtype t:  	time, needed by KIOPS
||N_Vector u:	the state vector
||N_Vector fu:	rhs evaluated at u
||void * pb:	pointer to my modified problem class
||N_Vector tmp:	temporary vector used here to grab rows of the Jacobian
\\=============================================================
*/
int JtV_TCHEM_Fast(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void* pb, N_Vector tmp)
{
	myPb * pbPtr{static_cast<myPb *> (pb)};	//Recast
	N_VConst(0.0, Jv);						//Clear data to not poison results
	g_timer.reset();						//Time... lets GO!
	MatrixVectorProductFast(pbPtr->num_equations, pbPtr->Jac, v, tmp, NV_DATA_S(Jv));//Run
	MatVecTime += g_timer.seconds();		//End timing
	JtvCnt ++;
	return 0;
}




//====== ____===========================//
// 		|  _ \     _				//
//  	| (_) )   (_)         _		//
//  	|  __/`._  _  _ __  _| |_	//
//  	| |  |  _)| || `  \(_   _)	//
//		| |  | |  | || |\ |  | |	//
//		|_|  |_|  |_||_||_|  |_|	//
//========================================
//===================================================
//Print data to output file
//myfile:	pointer to the filestream
//data:		pointer to the data
//number_of_equations
//t:		time/Stepsize
//===================================================
void PrintDataToFile(ofstream & myfile, realtype * data, int number_of_equations, realtype t)
{
	for (int i=0; i<number_of_equations; i++)
	{
		myfile<<setprecision(20)<<fixed <<data[i] <<"\t\t";
	}
	myfile << "\t\t " << t;
	myfile.flush();

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


//===============================
//  ___  ___  ___   ___   ___
// |  _)| _ \| _ \ /   \ | _ \
// |  _)|   /|   /|  X  ||   /
// |___)|_\_\|_\_\ \___/ |_\_\
//===============================
//RHS Error Check
//realtype* y 		y Data ptr
//realtype* rhs		rhs Data ptr
//int num_eqs		number of equations in the system
//========================================================
void ErrorCheck(ofstream & myfile, N_Vector y, realtype * data, int number_of_equations, realtype TNext)
{
	realtype MassError=abs( N_VL1NormLocal(y)-data[0]-1.0);
	if(MassError>.1 || abs(data[0])>1e5 )
	{
		PrintDataToFile(myfile, data, number_of_equations, TNext);
		cout<<"\n===============!!!Critical error!!!================\n";
		if(MassError>.1)
		{
			cout<<"Severe mass loss detected. Check output file for details.\n";
			myfile<<"Mass fraction error : " << abs( N_VL1NormLocal(y)-data[0]-1.0)<< endl;
		}
		if(abs(data[0])>1e5)
			cout<<"Temperature instability detected. Check output file for details.\n";
		cout<<"\nFailed computing during time: "<< TNext<<endl;
		exit(EXIT_FAILURE);
	}

}
/*
//=============================================================
//  __ __    ___   _____      __   __ ___   ___
// |  V  | ./ _ \.|_   _| ___ \ \ / /| __) / __)
// | .V. | | (_) |  | |  |___| \ v / | __)( (__
// |_| |_| |_/ \_|  |_|         \_/  |___) \___)
//============================================================
//We cover our bases by doing MatVec directly for EPI methods.
//By doing this, we don't have to transpose the data.
//Does Jac*v = JV
*/
void MatrixVectorProductFast(int len, N_Vector Jac, N_Vector x, N_Vector tmp, realtype * JV)
{
	realtype * X 			= NV_DATA_S(x);
	realtype * JacD			= NV_DATA_S(Jac);
	for (int i=0; i< len; i++)//for each state
	{
		for(int j=0; j<len; j++)//marches across the column
		{
			JV[i] += X[j] * JacD[j+ i * len];
		}
	}
}

void PrintBanner()
{
	cout << "\n\n\n\n\n\n";
	cout << "\t==========================================================\n";
	cout << "\t||  <:::::>  <::::>  <:::::>     <:::>         <::::>    ||\n";
	cout << "\t||     .*/   |::     |::  ::>   <>   <>        |::  `>   ||\n";
	cout << "\t||    .*/    |::::>  |:::::<   <>     <>  <::> |::   *>  ||\n";
	cout << "\t||   .*/     |::     |:: ':>    <>   <>        |::  .>   ||\n";
	cout << "\t||  <:::::>  <::::>  |::   :>    <:::>         <::::>    ||\n";
	cout << "\t==========================================================\n";

}
