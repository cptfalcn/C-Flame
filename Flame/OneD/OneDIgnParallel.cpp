/*
 * ===========================================================================================
 * This is a Parallel implementation.
 * ===========================================================================================
 */

#include <stdio.h>
#include <math.h>
#include <mpi.h>
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
#include <sundials/sundials_nvector.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TChem_CommandLineParser.hpp"

//These two will be needed in the future
#include "TChem_Impl_NewtonSolver.hpp"
//#include "TChem_KineticModelData.hpp"
//#include "TChem_Impl_TrBDF2.hpp"
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
void ErrorCheck(ofstream &, N_Vector, realtype *, int, realtype);
static int check_flag(void *, const string, int);
void PrintCVODEStats(void *, long int, int, long int, realtype, realtype, realtype, realtype);

//=====================
//Namespaces and globals
//=====================
using namespace std;

//Kokkos::Impl::Timer g_timer;
//====================
//Main
//====================
int main(int argc, char* argv[])
{
	PrintBanner();//Print One-D, how fancy
	//====================
	// Intial Declarations
	//====================
	long int NumStepsTaken = 0;
	int LastOrder = 0;
	long int OrderReductions = 0 ;
	realtype InitStep = 0;
	realtype LastStep = 0;
	realtype CurTime = 0;
	static realtype FinalTime = 0;//1.0e-4;//1e-4 seems to be the max
	static const int NumBands = 3;//Epic stuff, default is 3.
	realtype StepSize=0;
	realtype KTime=0;
	int ProgressDots=0; //From 0 to 100; only care about percentage
	int StepCount = 0;
	realtype PercentDone=0;//void *userData=nullptr; //Empty problem pointer
	realtype TNow=0;
	realtype TNext=0;
	string MyFile="Default.txt";
	string Method="EPI2";//This is the default method
	static realtype KrylovTol=1e-14;
	int UseJac=1; //will we use the Jacobian or not
	int SampleNum=0;
	int Experiment=1; //default is hydrogen, set to 0 for kapila
	int Profiling=0;//default to no profiling, need to edit profiling
	int number_of_equations=0;
	int NumScalarPoints = 0;
	realtype Delx = 5e-1;
	realtype absTol = 1e-8;
        realtype relTol = 1e-8;
	int startingBasisSizes[] = {3, 3};

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
	opts.set_option<realtype>("Delx", "Grid Delx", &Delx);
        opts.set_option<realtype>("relTol", "Solver Relative Tol", &relTol);
        opts.set_option<realtype>("absTol", "Solver Absolulte Tol", &absTol);

	const bool r_parse = opts.parse(argc, argv);
	if (r_parse)
		return 0; // print help return

	//======================================================
	//Check inputs and output to console for review.
	//======================================================
	ofstream myfile(MyFile, std::ios_base::app);
	int Steps=CheckStep(FinalTime, StepSize);//Checks the number of steps
	realtype StopTime = FinalTime;
	//=====================================================
	//Kokkos block: working, pruning old code
	//=====================================================
	Kokkos::initialize(argc, argv);
	//ALL Kokkos varialbes are reference counted objects.
	//They are deallocated within this local scope.
	{//begin local scope
		//=====================================
		//Kokkos Sub-declarations
		//=====================================
		// Define scalar type and ordinal type
		using real_type = double;
		using ordinal_type = int;

		// Kokkos environments - host device type and multi dimensional arrays
		// note: The 2d view uses row major layout.  Most matrix format uses column major layout.
		// to make the 2d view compatible with other codes, we need to transpose it.
		using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
		using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
		//using real_type_2d_view_type = Tines::value_type_2d_view<real_type,host_device_type>;

	    	// construct TChem's kinect model and read reaction mechanism
    		TChem::KineticModelData kmd(chemFile, thermFile);

	    	// construct const (read-only) data and move the data to device.
	    	// for host device, it is soft copy (pointer assignmnet).
	    	auto kmcd = kmd.createConstData<host_device_type>();

    		// use problem object for ignition zero D problem, interface to source term and jacobian
    		// other available problem objects in TChem: PFR, CSTR, and CV ignition
    		using problem_type = TChem::Impl::IgnitionZeroD_Problem<decltype(kmcd)>;

	    	//State vector - Temperature, Y_0, Y_1, ... Y_{n-1}), n is # of species.
		//End Kokkos sub declarations
		//=======================================
		//Set number of equations and grid points
		//=======================================
		if(Experiment==0)
			number_of_equations=3;
		else
			number_of_equations=problem_type::getNumberOfEquations(kmcd);

		if(Delx==0)
			NumScalarPoints = 1 ;
		else
			NumScalarPoints= round(1/Delx) + 2 ;//+3 for Velocity

		//==================================================
		//Prep/Set initial conditions, inherited from Zero-D
		//==================================================
		N_Vector State = N_VNew_Serial(number_of_equations*NumScalarPoints);
		N_Vector StateDot= N_VNew_Serial(number_of_equations*NumScalarPoints);
		N_Vector y = N_VNew_Serial(number_of_equations);//The data vector y
		N_VScale(0.0 , y, y);//Zero out the vector.
		N_Vector yDot = N_VNew_Serial(number_of_equations);
		realtype *data = NV_DATA_S(y);//Set the State data pointer
		realtype *StateData = NV_DATA_S(State);

		//Set Initial conditions based on experiment# & output to terminal
		switch(Experiment)
		{
			case 0:
				IntConKappa(data, 1e-2);
				break;
			case 1:
				IntConHydro(data, SampleNum);
        			break;
			case 2:
				IntConGri30(data, SampleNum);
				break;
		}

		//Set state based off intial conditions
		for( int i = 0 ; i < number_of_equations ; i ++ )
		{
			for (int j = 0 ; j < NumScalarPoints; j ++)
				StateData[j + i * NumScalarPoints] = data[i];
		}
		//End Initial conditions

		//=========================================================================
		//Set TChem
		//=========================================================================
    		//TChem does not allocate any workspace internally.  Workspace should be explicitly given
		//from users.  You can create the work space using NVector,  std::vector or
		//real_type_1d_view_type (Kokkos view). Here we use kokkos view. Set up an ignition
		//problem and the workspace. Unused in Experiment 0.
		//=========================================================================
		const ordinal_type problem_workspace_size = problem_type::getWorkSpaceSize(kmcd);
		real_type_1d_view_type work("workspace", problem_workspace_size);
	    	//set problem
		myPb problem;//Defined in Chemistry.h
		myPb2 problem2(number_of_equations, work, kmcd,NumScalarPoints);//Construct new Problem.
		problem2.PrintGuts();//Try a print
		void *pbptr= &problem;
		void *UserData = &problem2;

      		// initialize problem
		const real_type pressure(101325);//Constant pressure
		real_type_1d_view_type fac("fac", number_of_equations);//originally 2*
      		problem._p = pressure;	// pressure
		problem._fac= fac;	// Not sure what this is.
      		problem._work = work;  	// problem workspace array
     		problem._kmcd = kmcd;  	// kinetic model

		problem.num_equations=number_of_equations;

		problem.Jac=N_VNew_Serial(number_of_equations*number_of_equations);//Make the Jacobian
		int MaxKrylovIters = max(number_of_equations*number_of_equations,500);//500 is the base
		//End set Tchem
		//==============================================
		//Create integrators
		//==============================================
		//CVODE
		//==============
                void * cvode_mem;
		void * cvode_Test;
		SUNMatrix A 			= SUNDenseMatrix(number_of_equations, number_of_equations);
		SUNLinearSolver LS2 		= SUNLinSol_SPGMR(y, PREC_NONE,0);
                SUNNonlinearSolver NLS 		= SUNNonlinSol_Newton(y);
		SUNMatrix B			= SUNDenseMatrix(number_of_equations * NumScalarPoints,
							number_of_equations * NumScalarPoints);
		SUNLinearSolver SUPERLS 	= SUNLinSol_SPGMR(State, PREC_NONE, 0);
		SUNNonlinearSolver SUPERNLS 	= SUNNonlinSol_Newton(State);
                int retVal = 0;

		//================================
		//Set other integrators
		//================================
		Epi2_KIOPS 	*integrator 		= NULL;
		Epi3_KIOPS 	*integrator2		= NULL;
		EpiRK4SC_KIOPS 	*integrator3 		= NULL;
		void 		*GeneralIntegrator	= NULL;
		IntegratorStats *integratorStats 	= NULL;
		//==============================================================
		//Parse the experiment cases
		//Change this loop, the inputs to state. Also change_RHS and JTV
		//==============================================================
		if(Experiment!=0)//If using TChem problems
		{
			//Debugging gets put here for now.
			cout << BAR << "Debugging tests" << BAR << endl;
			GeneralIntegrator = SelectIntegrator(SUPER_CHEM_RHS_TCHEM, SUPER_CHEM_JTV,
					SUPER_CHEM_JAC_TCHEM,
					UserData, GeneralIntegrator, MaxKrylovIters, State,
					number_of_equations, UseJac,
					Method, B, SUPERLS, SUPERNLS, absTol, relTol, StepSize);
			CompareJacobians(pbptr, UserData, y, State, yDot, StateDot, problem2.Mat);
			IntegrateWrapper(UseJac, StepSize, 0, StepSize, NumBands, State, KrylovTol,
					startingBasisSizes, UserData, GeneralIntegrator, GeneralIntegrator,
					Method, integratorStats, SUPER_CHEM_JAC_TCHEM);
			PrintSuperVector(StateData, Experiment, NumScalarPoints);
			cout << BAR << "End debugging tests" << BAR << endl;
			//End debugging
			if(Method=="EPI2")
			{
				integrator = CreateEPI2Integrator(CHEM_RHS_TCHEM, CHEM_JTV, pbptr,
					MaxKrylovIters, y,number_of_equations, UseJac);
			}else if(Method=="EPI3")
			{
				integrator2 = CreateEPI3Integrator(CHEM_RHS_TCHEM, CHEM_JTV, pbptr,
					MaxKrylovIters, y, number_of_equations, UseJac);
			}else if(Method=="EPIRK4")
			{
				integrator3 = CreateEPIRK4SCIntegrator(CHEM_RHS_TCHEM, CHEM_JTV, pbptr,
					 MaxKrylovIters, y, number_of_equations, UseJac);
			}else if(Method=="CVODEKry")
                        {
				cvode_mem=CreateCVODEKrylov(CHEM_RHS_TCHEM, CHEM_JTV, pbptr, A, LS2, NLS,
					number_of_equations, y, relTol, absTol, StepSize, UseJac);
                                if(UseJac==0)
                                cout << BAR << "\tCVODE Kry W/o Jac\t" << BAR << endl;
                                else
                                cout << BAR << "\tCVODE Kry W Jac  \t" << BAR << endl;
                        }else if(Method == "CVODEOS")
			{
				cvode_mem=CreateCVODEKrylov(CHEM_RHS_TCHEM, CHEM_JTV, pbptr, A, LS2, NLS,
					number_of_equations, y, relTol, absTol, FinalTime, UseJac);
			}
		}else if(Experiment==0){
                        if(Method=="EPI2")
                        	integrator = CreateEPI2Integrator(RHS_KAPPA, Jtv_KAPPA, pbptr,
                                        MaxKrylovIters, y,number_of_equations, UseJac);
                        else if(Method=="EPI3")
                                integrator2 = CreateEPI3Integrator(RHS_KAPPA, Jtv_KAPPA, pbptr,
                                        MaxKrylovIters, y, number_of_equations, UseJac);
                        else if(Method=="EPIRK4")
                                integrator3 = CreateEPIRK4SCIntegrator(RHS_KAPPA, Jtv_KAPPA, pbptr,
                                         MaxKrylovIters, y, number_of_equations, UseJac);
                }
		//End  create integrators
		//===============
		//
		//Debugging tests
		//
		//===============
		if (Experiment !=0)//Experiment 0 is bugged, crashes if this is run.
		{
			//This is bugged, commenting out the test section causes the test run to hang.
			//If we do run it, the result is slightly off from what is expected.
			//The Jacobian is probably not being set correctly
			//SUPER_CHEM_RHS_TCHEM(0, State, StateDot, UserData);
			//SUPER_CHEM_JAC_TCHEM(0, State, StateDot, B, UserData,y,y,y);
			/*
			cout << BAR << "Test SUPERVECTOR integrator" << BAR << endl;
                        cvode_Test 	= CreateCVODEKrylov(SUPER_CHEM_RHS_TCHEM, SUPER_CHEM_JTV,
                                        UserData, B, SUPERLS, SUPERNLS,
                                        number_of_equations*NumScalarPoints, State, relTol,
                                        absTol, FinalTime, UseJac);
                        retVal = CVodeSetMaxStep(cvode_Test, FinalTime);
                        retVal = CVodeSetStopTime(cvode_Test, FinalTime);
                        retVal = CVodeSetJacTimes(cvode_Test, 0, SUPER_CHEM_JTV);
                        auto Start=std::chrono::high_resolution_clock::now();//Time integrator
			retVal = CVode(cvode_Test, FinalTime, State, &FinalTime, CV_NORMAL);
			auto Stop=std::chrono::high_resolution_clock::now();
                        auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
                        check_flag(&retVal, "CVODE", 1);
			PrintSuperVector(StateData, Experiment, NumScalarPoints);
			PrintCVODEStats(cvode_Test, NumStepsTaken, LastOrder, OrderReductions,
                       			CurTime, InitStep, LastStep, Pass.count()/1e9);
			cout << BAR << "Test Complete" << BAR << endl;
			FinalTime=StopTime;
			CVodeFree(&cvode_Test);
			*/
		}
        	//=================================
        	// Run the time integrator loop
        	//=================================
		if( Method == "CVODEOS")
		{//Skip the loop
			retVal = CVodeSetMaxStep(cvode_mem, FinalTime);
			retVal = CVodeSetStopTime(cvode_mem, FinalTime);
			auto Start=std::chrono::high_resolution_clock::now();//Time integrator
			CVode(cvode_mem, FinalTime, y, &FinalTime, CV_NORMAL);
			auto Stop=std::chrono::high_resolution_clock::now();
                        auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
                        KTime+=Pass.count()/1e9;
			PrintCVODEStats(cvode_mem, NumStepsTaken, LastOrder, OrderReductions,
                                        CurTime, InitStep, LastStep, Pass.count()/1e9);
		}else{//Do the loop
			cout << "Starting integration at t=" << StepCount*StepSize << endl <<"[";
			cout.flush();
			while(StepCount<Steps)//while(TNext<FinalTime)
        		{
                		TNow= StepCount*StepSize;
                		TNext=(StepCount+1)*StepSize;
				auto Start=std::chrono::high_resolution_clock::now();//Time integrator
				//Integrate
				integratorStats = Integrate(UseJac, StepSize, TNow, TNext, NumBands, y,
						KrylovTol, startingBasisSizes, pbptr, cvode_mem, Method,
						integrator, integrator2, integrator3, integratorStats);
				//Clock the time spent in the integrator
				auto Stop=std::chrono::high_resolution_clock::now();
                        	auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
				KTime+=Pass.count()/1e9;
				//Also run the super stuff
				/**/
				IntegrateWrapper(UseJac, StepSize, TNow, TNext, NumBands, State, KrylovTol,
                                         startingBasisSizes, UserData, GeneralIntegrator,GeneralIntegrator,
                                         Method, integratorStats, SUPER_CHEM_JAC_TCHEM);
				/**/
				//=======================================
				//Error checking
				//=======================================
				if(Experiment!=0) //Error checking invalid for Experiment 0.
					ErrorCheck(myfile, y, data, number_of_equations, TNext);
				//clean
                		for (int j=0; j<number_of_equations; j++)
             			{
					if(data[j]<0)
                        	        	data[j]=0;
                		}
				//=======================
				//Track the progress
				//=======================
                		PercentDone=floor(TNext/FinalTime*100);
                		for (int k=0; k<PercentDone-ProgressDots;k++)
                		{
                	        	cout<<":";
                	        	cout.flush();
                		}

				/*//Use if we want to do timeseries plots
				if(ProgressDots!=PercentDone)
					PrintDataToFile(myfile, data, number_of_equations,TNext);*/
				ProgressDots=PercentDone;
                		StepCount++;
        		}//End integration loop
        		TNow=TNext;
        		cout << "]100%\n\n";
		}
		cout << BAR << "\tIntegration complete\t" << BAR <<endl;
	        //=======================================
        	//Console Output
        	//=======================================
		PrintExpData(data, Experiment, N_VL1NormLocal(y));
		cout << BAR << "Printing SuperVector" << BAR << endl;
		PrintSuperVector(StateData, Experiment, NumScalarPoints);
		//General Simulation paramater output to console
		cout << BAR << "\tSim Parameters\t\t" << BAR << endl;
		PrintExpParam(FinalTime, TNow, StepSize, StepCount, KrylovTol, absTol, relTol, KTime);
		if (Profiling ==1){//Invalid for experiment 0
			cout << BAR << "\tPerformance data\t" << BAR << endl;
			integratorStats->PrintStats();
		}//End Profiling

		cout << BAR <<"Printing data to "<< MyFile << BAR << endl;
		//Print data to ouput file and close
                //PrintDataToFile(myfile,data,number_of_equations,StepSize);//Print here for conv studies
		PrintDataToFile(myfile, data, number_of_equations, absTol);//Print this for Tol Studies
		myfile << "\t\t" << KTime <<endl;//Print the integrator time
        	myfile.close();
		//==========================
		//Time to take out the trash
		//==========================
		cout << BAR << "Deleting objects" << BAR << endl;
		delete integrator;
		delete integrator2;
		delete integrator3;
		//delete problem.JacArr;
		cout << BAR << "Deleted Integrators" << endl;
		N_VDestroy_Serial(y);
		N_VDestroy_Serial(yDot);
		N_VDestroy_Serial(State);
		N_VDestroy_Serial(StateDot);
		//N_VDestroy_Serial(problem.Jac);
		cout << BAR << "Deleted N_Vectors" << endl;
		SUNMatDestroy(A);
		SUNMatDestroy(B);
		cout << BAR << "Deleted SUNMAtrices" << endl;
		SUNLinSolFree(LS2);
		SUNLinSolFree(SUPERLS);
		SUNNonlinSolFree(NLS);
		SUNNonlinSolFree(SUPERNLS);
		cout << BAR << "Deleted Solvers" << endl;
                CVodeFree(&cvode_mem);
		cout << BAR << "Deleted CVODE memory" << endl;
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
//===========================================
int CheckStep(realtype FinalTime, realtype StepSize)
{
	int Steps=0;
	cout << BAR << "\t Step Check\t\t" << BAR << endl ;
	cout << "Checking the proposed number of steps...\n";
	cout << std::setprecision(17) << FinalTime/StepSize << " steps proposed...";

        if(floor( FinalTime/StepSize ) == FinalTime/StepSize )
	{
                Steps= FinalTime/StepSize;
                cout << " accepted...\n";
		return Steps;
        }else if (abs(round(FinalTime/StepSize))-FinalTime/StepSize <1e-6 ){
		Steps= round (FinalTime/StepSize);
		cout << Steps << " steps approximated...\n";
		return Steps;
	}else{
                cout<<"Cannot perform non-integer number of steps!!\n";
                exit(EXIT_FAILURE);
        }
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
                   	myfile <<"Mass fraction error : " << abs( N_VL1NormLocal(y)-data[0]-1.0)<< endl;
              	}
             	if(abs(data[0])>1e5)
            		cout<<"Temperature instability detected. Check output file for details.\n";
             	cout<<"\nFailed computing during time: "<< TNext<<endl;
              	exit(EXIT_FAILURE);
	}
}


void PrintBanner()
{
	cout << "\n\n\n\n\n\n";
	cout << "\t================================================\n";
	cout << "\t||      <:::>   :>   ::  <::::>      <::::>    ||\n";
	cout << "\t||     <>   <>  |:>  ::  |::         |::  `>   ||\n";
	cout << "\t||    <>     <> |: > ::  |::::> <::> |::   *>  ||\n";
	cout << "\t||     <>   <>  |:  >::  |::         |::  .>   ||\n";
	cout << "\t||      <:::>   |:   >:  <::::>      <::::>    ||\n";
	cout << "\t================================================\n";

}
static int check_flag(void *flagvalue, const string funcname, int opt)
{
 	 // Check if the function returned a NULL pointer
  	if (opt == 0)
  	{
    		if (flagvalue == NULL)
    		{
      		cerr << endl << "ERROR: " << funcname << " returned NULL pointer" << endl << endl;
      		return 1;
    		}
  	}
  	// Check the function return flag value
  	else if (opt == 1 || opt == 2)
  	{
    		int errflag = *((int *) flagvalue);
    		if  ((opt == 1 && errflag < 0) || (opt == 2 && errflag != 0))
    		{
      			cerr << endl << "ERROR: " << funcname << " returned with flag = "
           		<< errflag << endl << endl;
      			return 1;
    		}
  	}
	else
	{
    		cerr << endl << "ERROR: check_flag called with an invalid option value" << endl;
    		return 1;
  	}
  	return 0;
}

void PrintCVODEStats(void * cvode_mem,  long int NumStepsTaken, int LastOrder, long int OrderReductions,
			realtype CurTime, realtype InitStep, realtype LastStep, realtype IntTime)
{
			int retVal = 0;
                        retVal = CVodeGetActualInitStep(cvode_mem, &InitStep);
                        retVal = CVodeGetCurrentTime(cvode_mem, &CurTime);
                        retVal = CVodeGetNumSteps(cvode_mem, &NumStepsTaken);
                        retVal = CVodeGetLastStep(cvode_mem,&LastStep);
                        retVal = CVodeGetLastOrder(cvode_mem, &LastOrder);
                        retVal = CVodeGetNumStabLimOrderReds(cvode_mem, &OrderReductions);
                        cout << "Last step size: " << LastStep << "\t\t";
                        cout << "Order of method: " << LastOrder << endl;
                        cout << "Initial Step: " << InitStep << "\t\t";
                        cout << "Current time: " << CurTime << endl;
                        cout << "Order Reductions: " << OrderReductions << "\t\t" ;
                        cout << "Number internal steps: " << NumStepsTaken << "\t\tInt Time:";
                        cout << IntTime  << " seconds" << endl;
}

