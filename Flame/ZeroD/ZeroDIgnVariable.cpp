/*
 * ===========================================================================================
 * This is serial implementation.
 * ===========================================================================================
 */
//This should be the new version of the ZeroDIgnVar, with faster timing
//This supports Epi3VChem_KIOPS and CVODE integrators
#include <stdio.h>
#include <math.h>
#include "Epic.h"
#include <nvector/nvector_serial.h>
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
#include "TChem_Impl_NewtonSolver.hpp"
#include "TChem_Impl_TrBDF2.hpp"
#include <chrono>
#include "InitialConditions.h"
#include "Chemistry.h"
#include "CreateIntegrators.h"
#include <cvode/cvode.h>
#include "Print.h"

//#define TCHEMPB TChem::Impl::IgnitionZeroD_Problem<TChem::KineticModelConstData<Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> >>
//#define TCHEMPB TChem::Impl::IgnitionZeroD_Problem	<TChem::KineticModelConstData	<Kokkos::Device	<Kokkos::Serial, Kokkos::HostSpace>  >	>
//Change to Serial, this can be changed by altering the TChem master build profile to include OPENMP on or off
#define BAR "===================="


//=====================
//Prototypes & Classes
//=====================
// class myPb : public TCHEMPB{
// 	public:
// 	//members
// 	ordinal_type  num_equations;
// 	N_Vector Jac;
// 	//functions
// 	ordinal_type	get_num_equations(void)
// 		return this->num_equations;
// };
//See Chemistry.h


//Banner
void PrintBanner();

//Initial conditions
//See InitialConditions.h

//RHS functions
int RHS_TCHEM(realtype, N_Vector, N_Vector, void *);

//JtV functions
int Jtv_TCHEM(N_Vector , N_Vector ,realtype, N_Vector, N_Vector , void* , N_Vector );

//Misc functions
int CheckStep(realtype, realtype);
void PrintFromPtr(realtype *,  int);
void ErrorCheck(ofstream &, N_Vector, realtype *, int, realtype);
void PrintDataToFile(ofstream &, realtype *,int, realtype);
void PrintToFile(ofstream &, realtype *, int, realtype, realtype, realtype, realtype, realtype, realtype);


//Used in Jtv
int ComputeJac(N_Vector, void*);
int CVodeComputeJacWrapper(realtype t, N_Vector u, N_Vector fy, SUNMatrix Jac,
               void * pb, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

//=====================
//Namespaces and globals
//=====================
using namespace std;

int JacCnt=0;
int RHSCnt=0;

double JacTime=0;
double RHSTime=0;
double MatVecTime=0;
Kokkos::Impl::Timer g_timer;

//====================
//Main
//====================
int main(int argc, char* argv[])
{
	PrintBanner();//Print Zero-D, how fancy
	//====================
	// Intial Declarations
	//====================
	static realtype FinalTime 	= 0;//1.0e-4;//1e-4 seems to be the max
	static const int NumBands 	= 3;//Epic stuff, default is 3.
	realtype StepSize			= 0;
	realtype KTime				= 0;
	realtype PrintTime			= 0;
	int ProgressDots			= 0; //From 0 to 100; only care about percentage
	int StepCount 				= 0;
	realtype PercentDone		= 0;
	realtype TNow				= 0;
	realtype TNext				= 0;
	string MyFile				= "Default.txt";
	string Method				= "EPI2";//This is the default method
	static realtype KrylovTol	= 1e-14;
	int UseJac					= 1; //will we use the Jacobian or not
	int SampleNum				= 0;
	int Experiment				= 1; //default is hydrogen, set to 0 for kapila
	int Profiling				= 0;//default to no profiling
	int number_of_equations		= 0;
	int startingBasisSizes[] 	= {10, 10};
	realtype relTol				= 1e-10;
	realtype absTol				= 1e-10;
	realtype maxSS 				= StepSize;
	int Movie 					= 0;

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
	opts.set_option<realtype>("maxSS", "Solver max step size", &maxSS);
	opts.set_option<int>("Movie", "Generate a data dump set at every step", &Movie);

	const bool r_parse = opts.parse(argc, argv);
	if (r_parse)
		return 0; // print help return

	ofstream myfile(MyFile, std::ios_base::app);
	//==========================
	//Start Kokkos
	//==========================
	/// all Kokkos varialbes are reference counted objects. 
	//  they are deallocated within this local scope.
	Kokkos::initialize(argc, argv);
	{//begin local scope
		//=====================================
		//Kokkos Sub-declarations
		//=====================================
		/// Kokkos environments - host device type and multi dimensional arrays
		/// note that the 2d view use row major layout while most matrix format uses column major layout.
		using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
		using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
		
		/// construct TChem's kinect model and read reaction mechanism
		TChem::KineticModelData kmd(chemFile, thermFile);

		/// construct const (read-only) data and move the data to device.
		/// for host device, it is soft copy (pointer assignmnet).
		auto kmcd = kmd.createConstData<host_device_type>();

		/// use problem object for ignition zero D problem, interface to source term and jacobian
		using problem_type = TChem::Impl::IgnitionZeroD_Problem<decltype(kmcd)>;

		//Set number of equations
		number_of_equations=problem_type::getNumberOfEquations(kmcd);

		//Declare the state variable locally
		N_Vector y = N_VNew_Serial(number_of_equations); //The data y
		N_VScale(0.0 , y, y);						//Zero out the y vector.
		//set the pointer to the state
		realtype *data = NV_DATA_S(y);
		SetIntCons(Experiment, SampleNum, data);
		/// TChem does not allocate any workspace internally.workspace should be explicitly given from users.
		/// you can create the work space using NVector, std::vector or real_type_1d_view_type (Kokkos view)
		/// here we use kokkos view
		//Set up an ignition problem and the workspace.
		const ordinal_type problem_workspace_size = problem_type::getWorkSpaceSize(kmcd);
		real_type_1d_view_type work("workspace", problem_workspace_size);
		myPb problem;
		void *pbptr= &problem;

		//========================
      	// initialize problem
		//=========================
		const real_type pressure(101325);//Constant pressure
		real_type_1d_view_type fac("fac", 2*number_of_equations);
      	problem._p 		= pressure; // pressure
		problem._fac	= fac;
      	problem._work 	= work;  // problem workspace array
     	problem._kmcd 	= kmcd;  // kinetic model

		problem.num_equations		= number_of_equations;
		problem.Jac					= N_VNew_Serial(number_of_equations*number_of_equations);//Make the Jacobian
		const int MaxKrylovIters 	= 1000;//500
		problem.MaxStepTaken		= 0;
		problem.MinStepTaken		= 0;
		problem.ignTime				= 0;
		problem.Movie				= 0;

		//==============================================
		//Create integrators
		//==============================================
		void * cvode_mem;
		SUNMatrix A						= SUNDenseMatrix(number_of_equations, number_of_equations);
		SUNLinearSolver LS 				= SUNLinSol_SPGMR(y, PREC_NONE, 20);
		SUNNonlinearSolver NLS 			= SUNNonlinSol_Newton(y);
		
		IntegratorStats *integratorStats = NULL;
		Epi3VChem_KIOPS * EPI3V = new Epi3VChem_KIOPS(RHS_TCHEM, Jtv_TCHEM, 
			CVodeComputeJacWrapper, pbptr, MaxKrylovIters, y, number_of_equations);

		int retVal = 0;
        N_Vector AbsTol= N_VNew_Serial(number_of_equations);
        for ( int i = 0 ; i < number_of_equations ; i++)
                NV_Ith_S(AbsTol,i)=absTol;
        cvode_mem = CVodeCreate (CV_BDF);				//Set BDF method
        retVal = CVodeSetUserData(cvode_mem, pbptr);	//Give CVODE the user problem
        retVal = CVodeInit(cvode_mem, RHS_TCHEM, 0, y); //Set intial state
        retVal = CVodeSetInitStep(cvode_mem, StepSize); //Set intial stepsize guess
        retVal = CVodeSetMaxStep(cvode_mem, maxSS);		//Set max stepsize
        retVal = CVodeSVtolerances(cvode_mem, relTol, AbsTol);	//Set tolerances
        //Set linear solver
        retVal = CVodeSetLinearSolver(cvode_mem, LS, A);//A might be null, try that
        retVal = CVodeSetNonlinearSolver(cvode_mem, NLS);
		retVal = CVodeSetJacTimes(cvode_mem, NULL, Jtv_TCHEM);
		retVal = CVodeSetLSetupFrequency(cvode_mem, 1); //Remove because also stupid.
		retVal = CVodeSetJacEvalFrequency(cvode_mem, 1);//Remove becuase this is stupid.
		retVal = CVodeSetJacFn(cvode_mem, CVodeComputeJacWrapper);		//Error if removed .
		retVal = CVodeSetMaxNumSteps(cvode_mem, 1e6);
		retVal = CVodeSetMaxOrd(cvode_mem,2);
        
		//===============================
		//Select and run integrator
		//===============================
		auto Start=std::chrono::high_resolution_clock::now();//Time integrator
		if(Method == "EPI3V")
		{
			cout << BAR <<"\tEPI3V Selected\t\t" << BAR <<endl;
			integratorStats = EPI3V->Integrate(StepSize, maxSS, absTol, relTol, 
								0.0, FinalTime, NumBands, startingBasisSizes, y);
			TNext = problem.t;
		}
		if(Method =="CVODEKry")
		{
			cout << BAR <<"\tCVODE Selected\t\t" << BAR <<endl;
			CVode(cvode_mem, FinalTime, y, &TNext, CV_NORMAL);
		}
		auto Stop=std::chrono::high_resolution_clock::now();
		auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
		KTime+=Pass.count()/1e9;

		//=======================================
		//Console Output
		//=======================================
		//Profiling output
		cout << BAR <<"Printing data to "<< MyFile << BAR << endl;
		PrintExpParam(FinalTime, TNext, StepSize, StepCount, KrylovTol, absTol, relTol, KTime, BAR);
		PrintSuperVector(data, Experiment, 1, BAR);
		cout << "Mass Fraction error: "<<abs( N_VL1NormLocal(y)-data[0]-1.0)<<endl;
		PrintProfiling(integratorStats, Profiling, Method,  BAR, cvode_mem);
		PrintToFile(myfile, data, number_of_equations, problem.t, relTol, absTol, KTime, KrylovTol, problem.ignTime);
		if (Profiling ==1){//Invalid for experiment 0
			ofstream ProFile("Profiling.txt", std::ios_base::app);//Profiling  File
			cout << "Integrator CPU time: "<<KTime<<" seconds\n";
	
			cout << BAR << "Jacobian\t" << BAR << endl;
			cout << "Jacobian Calls: " << JacCnt;
			ProFile << JacCnt << "\t" << JacTime << "\t";
			cout << "\t\tJacobian Time: " << JacTime << " seconds\n";
			cout << JacTime*100/KTime<<"%\n";

			//cout << "\t\tAveraged Jac Time cost: " <<JacTime/JacCnt << " seconds\n";
			cout << "MatVec Time: " << MatVecTime << " seconds\n";
			cout << MatVecTime*100/KTime <<"%\n";

			ProFile << RHSCnt << "\t" << RHSTime;
			cout << BAR << "RHS\t\t" << BAR << endl;
			cout << "RHS Function calls: " << RHSCnt << "\t\tRHS time: ";
			cout << RHSTime <<" seconds\n";
			cout << RHSTime/KTime*100<< "%\n";
			//cout << "\t\tAveraged RHS Time cost: " << RHSTime/RHSCnt << " seconds\n\n";
			//cout << "\tRelative cost ratio Jac/RHS: ";
			//cout << (JacTime/JacCnt)/(RHSTime/RHSCnt) << endl;

			//integratorStats->PrintStats();
			ProFile << "\t " << Experiment << "\t" << SampleNum << "\t";
			ProFile << FinalTime << "\t" << StepSize << endl;
			ProFile.close();

		}//End Profiling
		//==========================
		//Time to take out the trash
		//==========================
		myfile.close();
		delete EPI3V;
		N_VDestroy_Serial(AbsTol);
		N_VDestroy_Serial(y);
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

//===========================================
//Check the number of steps
//FinalTime
//StepSize
//===========================================
int CheckStep(realtype FinalTime, realtype StepSize)
{
	cout << BAR << "\t Step Check\t\t" << BAR << endl ;
	cout<<"Checking the proposed number of steps...\n";
	cout<<std::setprecision(17);
	cout<<FinalTime/StepSize << " steps proposed...";
	if(floor( FinalTime/StepSize ) == FinalTime/StepSize )
	{
		static const int Steps= FinalTime/StepSize;
		cout << " accepted...\n";
		return Steps;
	}
	else if (abs(round(FinalTime/StepSize))-FinalTime/StepSize <1e-6 )
	{
		static const int Steps= round (FinalTime/StepSize);
		cout << Steps << " steps approximated...\n";
		return Steps;
	}
	else
	{
		cout<<"Cannot perform non-integer number of steps!!\n";
		exit(EXIT_FAILURE);
	}
}



	//====================================================
	//  ___    _   _   ___
	// |   \  | | | | /   \
	// | |) ) | |_| | \ \\/
	// |   <  |  _  |  \ \
	// | |\	\ | | | | /\\ \
	// |_| \_\|_| |_| \___/
	//====================================================

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
int Jtv_TCHEM(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void* pb, N_Vector tmp)
{
	//problem_type problem;
	myPb * pbPtr{static_cast<myPb *> (pb)};//Recast
	ordinal_type number_of_equations=pbPtr->num_equations;//Get number of equations
	//======================================
	//Set the necessary vectors and pointers
	//======================================
	realtype *y= NV_DATA_S(u);
	realtype *JacD= NV_DATA_S(pbPtr->Jac);
	realtype *JV=NV_DATA_S(Jv);
	realtype *TMP=NV_DATA_S(tmp);
	//===================
	//Compute JV
	//===================
	g_timer.reset();
	MatrixVectorProduct(number_of_equations, JacD, v, tmp, JV);
	MatVecTime+=g_timer.seconds();
	return 0;
}


//====== ____===========================//
// 		|  _ \     _			//
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

//=======================================================
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
