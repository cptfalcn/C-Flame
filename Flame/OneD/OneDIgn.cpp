/*
 * ===========================================================================================
 * This is serial implementation.
 * ===========================================================================================
 */

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
#include <sundials/sundials_nvector.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TChem_CommandLineParser.hpp"

//These two will be needed in the future
//#include "TChem_Impl_NewtonSolver.hpp"
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
int TrackProgress(realtype, realtype, realtype, int);
void TrackSlowDown(IntegratorStats *, int *,  int, int *, realtype *, realtype, int*);
void TrackCVODEKryIters(SUNLinearSolver, realtype *, realtype, int *, int*, int*);
void ReadData(N_Vector, string);
void CheckNanBlowUp(N_Vector, int);

//=====================
//Namespaces and globals
//=====================
using namespace std;
//====================
//Main
//====================
int main(int argc, char* argv[])
{
	PrintBanner();							//Print One-D, how fancy
	//====================
	// Intial Declarations
	//====================
	int SlowDown 		= 0,	OldProjections	= 0;
	int MaxIters	= 0,	OldIters	= 0,	CurrStepIters	 = 0, TotalIters = 0;
	realtype SlowTime		= 0;
	static realtype FinalTime = 0;			//1.0e-4;//1e-4 seems to be the max
	static const int NumBands = 3;			//Epic stuff, default is 3, but unused.
	realtype StepSize		= 0;
	realtype KTime			= 0;
	int ProgressDots		= 0; 			//From 0 to 100; only care about percentage
	int StepCount 			= 0;
	realtype PercentDone	= 0;
	realtype TNow			= 0;
	realtype TNext			= 0;
	realtype TNextC			= 0;			//CVODE TNext, has a tendency to change TNext
	string MyFile			= "Default.txt";//Default output file
	string Method			= "EPI2";		//This is the default method
	static realtype KrylovTol=1e-14;
	int UseJac				= 1; 			//weither will we use the Jacobian or not
	int SampleNum			= 0;
	int Experiment			= 1; 			//default is hydrogen
	int Profiling			= 0;			//default to no profiling, need to edit profiling
	int number_of_equations	= 0;
	int NumScalarPoints 	= 0;
	realtype Delx 			= 1e-5;			//10 micrometer discretization
	realtype absTol 		= 1e-8;			//tight defaults 
    realtype relTol 		= 1e-8;			//tight default 
	int startingBasisSizes[]= {10, 10};		//{3,3}
	int VelUp 				= 1;			//Do we update Velocity, refactor to be yes always
	int Movie				= 0;
	int Restart				= 0;			//NULL implementation
	string InitialData		= "InitialData.txt";	// ""

	realtype ADV 			= 1.0;			//Fixed defaults for RHS
	realtype DIFF 			= 1.0;
	realtype CHEM 			= 1.0;
	realtype POW			= 0.0;

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
	opts.set_option<realtype>("Delx", "Grid Delx", &Delx);
	opts.set_option<realtype>("relTol", "Solver Relative Tol", &relTol);
	opts.set_option<realtype>("absTol", "Solver Absolulte Tol", &absTol);
	opts.set_option<realtype>("Adv", "Advection coefficient", &ADV);
	opts.set_option<realtype>("Diff", "Diffusion coefficient", &DIFF);
	opts.set_option<realtype>("Chem", "Chemistry coefficient", &CHEM);
	opts.set_option<realtype>("Pow", "Power coefficient", &POW);
	opts.set_option<int>("VelUp", "Do we do velocity up?", &VelUp);
	opts.set_option<int>("NumPts", "Interior points for each grid", &NumScalarPoints);
	opts.set_option<int>("Movie", "Generate a data set at every step", &Movie);
	opts.set_option<int>("Restart", "Restart from text file?", &Restart);

	const bool r_parse = opts.parse(argc, argv);
	if (r_parse)
		return 0; // print help return

	ofstream myfile(MyFile, std::ios_base::app);
	int Steps=CheckStep(FinalTime, StepSize);		//Checks the number of steps
	realtype VelStep = min(1e-5, FinalTime);

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
		//using real_type_2d_view_type = Tines::value_type_2d_view<real_type,host_device_type>;
		// construct TChem's kinect model and read reaction mechanism
		TChem::KineticModelData kmd(chemFile, thermFile);

		// construct const (read-only) data and move the data to device.
		// for host device, it is soft copy (pointer assignmnet).
		auto kmcd = kmd.createConstData<host_device_type>();

		// use problem object for ignition zero D problem, interface to source term and jacobian
		using problem_type = TChem::Impl::IgnitionZeroD_Problem<decltype(kmcd)>;

		//=======================================
		//Set number of equations and grid points
		//=======================================
		number_of_equations=problem_type::getNumberOfEquations(kmcd);

		//========================================
		//Parameter checking block
		//========================================
		if(Delx<=0.0)
		{
			std :: cout << "Invalid Spatial Step! Must be greater than 0.0. Exiting\n";
			exit(EXIT_FAILURE);
		}
		if(NumScalarPoints <= 1)
		{
			std :: cout << "You must have at least 2 grid points! Exiting\n";
			exit(EXIT_FAILURE);
		}

		int vecLength = number_of_equations * NumScalarPoints;	//For readability
		int num_eqs = number_of_equations;			//To reduce size
		int num_pts = NumScalarPoints;				//To reduce size
		//==================================================
		//Prep/Set initial conditions, inherited from Zero-D
		//==================================================
		N_Vector State 		= N_VNew_Serial(vecLength);	//State vector
		N_Vector StateDot	= N_VNew_Serial(vecLength);	//State Derivative RHS vector
		N_Vector y 			= N_VNew_Serial(num_eqs);	//intial conditions.
		N_VScale(0.0 , y, y);							//Zero out the y vector.
		N_VScale(0.0, State, State);					//Zero out the State vector.
		N_VScale(0.0, StateDot, StateDot);				//Zero out the StateDot vector.
		realtype *data 		= NV_DATA_S(y);				//Set the y data pointer.
		realtype *StateData = NV_DATA_S(State);			//Set the State data pointer.

		//Do an if statement read in here.
		if(Restart == 1)
			ReadData(State, InitialData);
		else
		{
			//Initial Conditions sub-block
			SetIntCons(Experiment, SampleNum, data);		//Set Initial Conditions
			SetSuperInitCons(data, StateData, num_eqs, num_pts);//Copy IC to State.
			TestingInitCons(SampleNum, num_pts, num_eqs, vecLength, Delx, data, StateData);//Keep until final version
		}
		//===================================================
		//Set TChem
		//===================================================
		//TChem does not allocate any workspace internally. You create the work space using NVector,
		//std::vector or real_type_1d_view_type (Kokkos view). Here we use kokkos view.
		//Set up an ignition problem and the workspace. Unused in Experiment 0.
		//=========================================================================
		const int problem_workspace_size = problem_type::getWorkSpaceSize(kmcd);
		real_type_1d_view_type work("workspace", problem_workspace_size);

		myPb2 problem2(num_eqs, work, kmcd, num_pts, y, Delx);	//Construct new Problem.
		problem2.SetGhostPVel(y, Experiment, SampleNum, 0);	//Do additional set up.
		problem2.SetAdvDiffReacPow(ADV, DIFF, CHEM, POW, VelUp);//Additional set up.
		problem2.kmd = kmd;

		//read in data
		ReadData(problem2.CPPoly,		"Cp_fT.txt");
		ReadData(problem2.TempTable,	"Tf.txt");
		ReadData(problem2.RhoTable, 	"rho_fT.txt");
		ReadData(problem2.DiffTable,	"Diff_fT.txt");
		//problem2.VerifyTempTable(State);

		
		void *UserData = &problem2;
		int MaxKrylovIters = 500; //max(vecLength, 500);//500 is the base
		//==================
		//Create integrators
		//==================
		//CVODE
		void * cvode_mem;
		SUNMatrix A							= SUNDenseMatrix(vecLength, vecLength);
		SUNLinearSolver SUPERLS 			= SUNLinSol_SPGMR(State, PREC_NONE, 20);
		SUNNonlinearSolver SUPERNLS 		= SUNNonlinSol_Newton(State);
		//Set EPI_KIOPS methods
		Epi2_KIOPS	*Epi2					= NULL;
		Epi3_KIOPS 	*Epi3					= NULL;
		EpiRK4SC_KIOPS	*EpiRK4				= NULL;
		IntegratorStats *integratorStats 	= NULL;
		//===================================================
		//Parse the experiment cases and make the integrators
		//===================================================
		Epi2 = 	new Epi2_KIOPS(SUPER_RHS,SUPER_JTV,UserData,MaxKrylovIters,State,vecLength);

		Epi3 = 	new Epi3_KIOPS(SUPER_RHS,SUPER_JTV,UserData,MaxKrylovIters,State,vecLength);

		EpiRK4=	new EpiRK4SC_KIOPS(SUPER_RHS,SUPER_JTV, UserData, MaxKrylovIters, State,vecLength);

		cvode_mem= CreateCVODE(SUPER_RHS, SUPER_JTV, SUPER_CHEM_JAC_TCHEM, UserData, A,
					SUPERLS, SUPERNLS, vecLength, State, relTol, absTol, StepSize, 1);

		PrintPreRun(StepSize, Delx, Steps, KrylovTol, absTol, relTol, Method, num_pts, BAR);
		//if(problem2.NumGridPts>1)
		//	problem2.RunTests(State);
        //=================================
        // Run the time integrator loop
        //=================================
		cout <<"[";
		cout.flush();
		while(StepCount<Steps)//while(TNext<FinalTime)
        {
            TNow		= StepCount*StepSize;
			TNext		= TNow + StepSize;
			TNextC 		= TNext;
			problem2.t	= TNow;
			//Integrate
			auto Start	=std::chrono::high_resolution_clock::now();//Time integrator
			if(Method != "CVODEKry")//Set the step Jacobian
				SUPER_CHEM_JAC_TCHEM(TNow, State, StateDot, A, UserData, State, State, State);
			if(Method == "EPI2")
				integratorStats =Epi2->Integrate(StepSize, TNow, TNext, NumBands,
					State, KrylovTol, startingBasisSizes);
			else if(Method == "EPI3")
				integratorStats =Epi3->Integrate(StepSize, TNow, TNext, NumBands,
					State, KrylovTol, startingBasisSizes);
			else if(Method == "EPIRK4")
				integratorStats =EpiRK4->Integrate(StepSize, TNow, TNext, NumBands,
					State, KrylovTol, startingBasisSizes);
			else if(Method == "CVODEKry")
				CVode(cvode_mem, TNext, State, &TNextC, CV_NORMAL);//CV_NORMAL/CV_ONE_STEP
			auto Stop=std::chrono::high_resolution_clock::now();
            auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
			KTime+=Pass.count()/1e9;
			//Step checking, factor into Error check
			if( TNext != TNextC)
			{
				std :: cout << "CVODE step issue\n";
				std :: cout << "TNext : " << TNext << "\t\tTNextC: " << TNextC <<endl;
				exit(EXIT_FAILURE);
			}
			//Error check
			CheckNanBlowUp(State, vecLength);

			//Profiling
			if(Profiling == 1)
			{
				if(Method != "CVODEKry")
					TrackSlowDown(integratorStats, &SlowDown, StepCount,
						&OldIters, &SlowTime, TNext, & MaxIters);
				if(Method == "CVODEKry")
					TrackCVODEKryIters(SUPERLS, &SlowTime,TNext, &OldIters,
						&MaxIters, &TotalIters);
			}

			//Clean
			Clean(vecLength, StateData);

			//Vel Update
			if(VelUp==1)
				problem2.UpdateOneDVel(State);

			//Check heating
			problem2.CheckHeating(State, TNow);

			//Track Progress
			ProgressDots=TrackProgress(FinalTime, TNext, PercentDone, ProgressDots);

			if(Movie ==1 ) //Use if we want a time-series plot
			{
				PrintDataToFile(myfile, StateData,vecLength, absTol, BAR, MyFile, TNext);
				PrintDataToFile(myfile, N_VGetArrayPointer(problem2.Vel), num_pts+1, 0, BAR, MyFile, 0);
			}
            StepCount++;
        }//End integration loop
		TNow=TNext;
		cout << "]100%\n" << BAR << "\tIntegration complete\t" << BAR <<endl;
		//=======================================
		//Various testing
		//=======================================
		//problem2.VerifyHeatingExp(State, State0, FinalTime);
		//=======================================
    	//Console Output
    	//=======================================
		cout << BAR <<"Printing data to "<< MyFile << BAR << endl;
		PrintExpParam(FinalTime, TNow, StepSize, StepCount, KrylovTol, absTol, relTol, KTime, BAR);
		PrintSuperVector(StateData, Experiment, NumScalarPoints, BAR);
		PrintProfiling(integratorStats, Profiling, Method,  BAR, cvode_mem);
		PrintDataToFile(myfile, StateData, vecLength, absTol, BAR, MyFile, KTime);//change  4th
		if(Profiling == 1)//Refactor into profiling later.
		{
			if(Method == "EPI2")
				cout << "Max Kry Iters: " << MaxIters << " at " << SlowTime << endl;
			if(Method == "CVODEKry")
			{
				cout << "Max Krylov Iterates: " << MaxIters << " at " << SlowTime <<endl;
				cout << "Total Krylov iterates: " << TotalIters<< endl;
			}
			cout << "\nTotal integration time: " << KTime << " sec " << endl;
			cout << "Total time spent in rhs: " << problem2.rhsTime << " sec" << endl;
			cout << "\t Diff rhs: " <<problem2.rhs_Diff << " sec\n";
			cout << "\t\t LookupTime: " << problem2.rhsDiffLookTime << " sec\n";
			cout << "\t Adv rhs: " << problem2.rhs_Adv << " sec\n";
			cout << "\t Chem rhs: " << problem2.rhs_Chem << " sec\n";
			cout << "Total time spent in jtv: " << problem2.jacTime << " sec" << endl;
			cout << "\t Diff jtv: " <<problem2.jtv_Diff << " sec\n";
			cout << "\t\t LookupTime: " << problem2.jtvDiffLookTime << " sec\n";
			cout << "\t Adv jtv: " << problem2.jtv_Adv << " sec\n";
        	cout << "\t Chem jtv: " << problem2.jtv_Chem << " sec\n";
			cout << "\t\t Move time: " << problem2.dataMoveTime << " sec\n";
			cout << "Total time spent in jacMake: " << problem2.jacMakeTime << " sec\n";
		}
		//if(problem2.NumGridPts > 1)
		PrintDataToFile(myfile, N_VGetArrayPointer(problem2.Vel), num_pts+1, 0, BAR, MyFile, 0);
        myfile.close();
		//==========================
		//Time to take out the trash
		//==========================
		delete Epi2;
		delete Epi3;
		delete EpiRK4;
		N_VDestroy_Serial(y);
		N_VDestroy_Serial(State);
		N_VDestroy_Serial(StateDot);
		SUNMatDestroy(A);
		SUNLinSolFree(SUPERLS);
		SUNNonlinSolFree(SUPERNLS);
        CVodeFree(&cvode_mem);
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

void TrackSlowDown(IntegratorStats * integratorStats, int * SlowDown, int StepsTaken, int * OldIters,
		realtype *SlowTime, realtype TNow, int * MaxIter)
{
	int CurrIters= (integratorStats->krylovStats->numIterations) - *OldIters;
	if(CurrIters >= integratorStats->krylovStats->maxIterations)
	{//if New iterations is greater than the old average, we are slowing down
		*SlowDown = 1;
		*SlowTime = TNow;			//Mark the slowdown time
		*MaxIter= CurrIters;
	}
	*OldIters = CurrIters;
	*OldIters = integratorStats->krylovStats->numIterations;
}

void TrackCVODEKryIters(SUNLinearSolver LS, realtype * SlowDownTime, realtype TNow, int * OldIters,
			int * MaxIters, int * TotalIters)
{
	if(  *OldIters <= SUNLinSolNumIters(LS)  )
	{
		*SlowDownTime=TNow;
		*MaxIters = SUNLinSolNumIters(LS);
	}
	*OldIters = SUNLinSolNumIters(LS);
	*TotalIters += *OldIters;				//Tracks all iters
}

void CheckNanBlowUp(N_Vector State, int vecLength)
{
	realtype * DATA = N_VGetArrayPointer(State);
	for( int i = 0 ; i < vecLength; i ++)
	{
		if(isnan(DATA[i]) || abs(DATA[i]) >1e5)
		{
			if( isnan(DATA[i]) )
				std :: cout << "Data has NaN: " <<  i << DATA[i] << "\n";
			if(abs(DATA[i]) > 1e5)
				std :: cout << "Data is out of control: "<< i << DATA[i] <<"\n";
			exit(EXIT_FAILURE);

		}
	}


}

//Currently writing, uncompiled.
void ReadData(N_Vector target, string fileName)
{
	realtype * ptr	= N_VGetArrayPointer(target);
	int len		= N_VGetLength(target);
	string line;
	ifstream file (fileName);			//Open the file
	for( int i =0 ; i < len; i ++)
	{
		getline(file, line);
		ptr[i] = std :: stod(line);
	}
	file.close();
}
