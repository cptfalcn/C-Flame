/*
 * ===========================================================================================
 * This is serial implementation.
 * ===========================================================================================
 */
//ToDo:  
//Make a new Initial conditions file that for Hydrogen and NButane(x)
//
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
#include <omp.h>
#include <optional>
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/transport.h"
#include <chrono>
#include "InitialConditions.h"
#include "ChemistryFinal.h"
#include "Epi2RetCode.h"
#include "Print.h"
#include "Epi3VOneD.h"

//#define TCHEMPB TChem::Impl::IgnitionZeroD_Problem<TChem::KineticModelConstData<Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> >>
//#define TCHEMPB TChem::Impl::IgnitionZeroD_Problem	<TChem::KineticModelConstData	<Kokkos::Device	<Kokkos::Serial, Kokkos::HostSpace>  >	>
//Change to Serial, this can be changed by altering the TChem master build profile to include OPENMP on or off
#define BAR "===================="

//=====================
//Prototypes
//=====================
int OneD_RHS			(realtype, N_Vector, N_Vector, void *);
int OneD_RHS_Adv		(realtype, N_Vector, N_Vector, void *);
int OneD_RHS_Chem		(realtype, N_Vector, N_Vector, void *);
int OneD_RHS_Diff		(realtype, N_Vector, N_Vector, void *);
int OneD_RHS_CrossDiff	(realtype, N_Vector, N_Vector, void *);
int OneD_RHS_Heat		(realtype, N_Vector, N_Vector, void *);

//Create Jacobians
int OneD_Jac			(realtype, N_Vector, N_Vector, SUNMatrix, void *, N_Vector,
						N_Vector, N_Vector);
int OneD_JacArray		(realtype, N_Vector, N_Vector, SUNMatrix, void *, N_Vector,
						N_Vector, N_Vector);

int OneD_JtV			(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
int OneD_JtV_Adv		(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
//int OneD_JtV_Diff		(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
int OneD_JtV_DiffFast	(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);

int OneD_JtV_CrossDiff	(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
int OneD_Jtv_Chem		(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);

int OneD_Jtv_ChemArray	(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);

//Finite difference RHS
int OneD_RHS_FD			(realtype, N_Vector, N_Vector, void *);

int OneD_VelDivergence 	(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);


//Banner
void PrintBanner();

//Misc functions
int CheckStep(realtype, realtype);
int TrackProgress(realtype, realtype, realtype, int);
void TrackSlowDown(IntegratorStats *, int *,  int, int *, realtype *, realtype, int*);
void TrackCVODEKryIters(SUNLinearSolver, realtype *, realtype, int *, int*, int*);
void ReadData(N_Vector, string);
void CheckNanBlowUp(N_Vector, int);

int LocalInitCons(realtype *, int, int);
int LocalPrint(realtype *, int, int, string); //Print the data for n-butane and hydrogen

//=======================================================
//SolArray class which interfaces with the problem class.
//=======================================================
int Set_ThermTransData(myPb2*, N_Vector, std::shared_ptr<Cantera::Solution>);
int Set_Derivatives(myPb2*, N_Vector, std::shared_ptr<Cantera::Solution>);
//=====================
//Namespaces and globals
//=====================
using namespace std;

int JacCnt=0;
int RHSCnt=0;
int JtvCnt=0;
int JtvChm=0;
int JtvAdv=0;
int JtvDif=0;
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
	int MaxIters		= 0,	OldIters	= 0,	CurrStepIters	 = 0, TotalIters = 0;
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
	int PressMult			= 1;			//1 atm
	string InitialData		= "InitialData.txt";	// ""
	string Mechanism 		= "gri3.0/gri30.yaml";
	string Mechphase		= "gri30";
	realtype HeatEnd		= 1e-4;

	realtype ADV 			= 1.0;			//Fixed defaults for RHS
	realtype DIFF 			= 1.0;
	realtype CHEM 			= 1.0;
	realtype POW			= 0.0;

	int	StepRate			= 1;

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
	opts.set_option<std::string>("Mechanism", "The path to the yaml mech", & Mechanism);
	opts.set_option<std::string>("Mechphase", "name of the mechanism phase", & Mechphase);
	opts.set_option<int>("Experiment", "Experiment chosen", &Experiment);
	opts.set_option<int>("Profiling", "Output profiling data or not", &Profiling);
	opts.set_option<realtype>("Delx", "Grid Delx", &Delx);
	opts.set_option<realtype>("relTol", "Solver Relative Tol", &relTol);
	opts.set_option<realtype>("absTol", "Solver Absolulte Tol", &absTol);
	opts.set_option<realtype>("Adv", "Advection coefficient", &ADV);
	opts.set_option<realtype>("Diff", "Diffusion coefficient", &DIFF);
	opts.set_option<realtype>("Chem", "Chemistry coefficient", &CHEM);
	opts.set_option<realtype>("Pow", "Power coefficient", &POW);
	opts.set_option<int>("VelUp", "Do we do velocity up?", &VelUp);			//In final version remove this to make it fixed
	opts.set_option<int>("NumPts", "Interior points for each grid", &NumScalarPoints);
	opts.set_option<int>("Movie", "Generate a data set at every step", &Movie);
	opts.set_option<realtype>("HeatEnd", "When to stop heating", &HeatEnd);
	opts.set_option<int>("StepRate", "How many steps to skip when making a video", &StepRate);
	const bool r_parse = opts.parse(argc, argv);
	if (r_parse)
		return 0; // print help return

	ofstream myfile(MyFile, std::ios_base::app);
	int Steps=CheckStep(FinalTime, StepSize);		//Checks the number of steps
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
		// N_VScale(0.0, LittleJac, LittleJac);			//Zero out LittleJac.

		realtype *data 		= NV_DATA_S(y);				//Set the y data pointer.
		realtype *StateData = NV_DATA_S(State);			//Set the State data pointer.

		//Initial Conditions sub-block
		SetIntCons(Experiment, SampleNum, data);		//Set Initial Conditions
		LocalInitCons(data, Experiment, SampleNum);		//Custom for this project

		SetSuperInitCons(data, StateData, num_eqs, num_pts);//Copy IC to State.
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
		problem2.Set_RHS(OneD_RHS_Chem, OneD_RHS_Adv, OneD_RHS_Diff, OneD_RHS_Heat);
		problem2.Set_Jtv(OneD_JtV_Adv, OneD_Jtv_ChemArray, OneD_JtV_CrossDiff, OneD_JtV_DiffFast);
		problem2.RHS_CrossDiff=OneD_RHS_CrossDiff;//Attach optional cross diffusion for testing.
		problem2.JtV_CrossDiff=OneD_JtV_CrossDiff;//Attach optional cross diffusion jtv for testing.
		problem2.RHS_Full=OneD_RHS;
		problem2.pb._p= 4*101325;
		std :: cout << "Pressure (atm): " << problem2.pb._p << " (" << problem2.pb._p/101325  <<") " << std :: endl;
		//===================
		//Create Cantera
		//===================
		//Set a solution per the gri standard
		std :: cout << "Using Mechanism: " << Mechanism << " and phase " << Mechphase << std :: endl;
		std::shared_ptr<Cantera::Solution> sol = Cantera::newSolution(Mechanism, Mechphase);
		Set_ThermTransData(&problem2, State, sol);
		//std::cout << sol->thermo()->report() << std::endl; //Remove the report in the final version

		//Set the pointer to the object into the problem
		problem2.sol = sol;
		//End Cantera
		void *UserData = &problem2;
		int MaxKrylovIters = 500;//max(vecLength, 500);//500 is the base
		//==================
		//Create integrators
		//==================
		//CVODE
		void * cvode_mem;
		int retVal=0;
		SUNMatrix A							= SUNDenseMatrix(vecLength, vecLength);
		SUNLinearSolver SUPERLS 			= SUNLinSol_SPGMR(State, PREC_NONE, 20);
		SUNNonlinearSolver SUPERNLS 		= SUNNonlinSol_Newton(State);

		
		//Set EPI_KIOPS methods
		Epi2_KIOPS		*Epi2					= NULL;
		//Epi2_KIOPS		*Epi2_FD				= NULL;
		Epi3VChem_KIOPS *Epi3V 					= NULL;
		Epi2RetCode  	*Epi2_FD				= NULL;
		IntegratorStats *integratorStats 	= NULL;
		//===================================================
		//Parse the experiment cases and make the integrators
		//===================================================
		Epi2 	= new Epi2_KIOPS(OneD_RHS,OneD_JtV,UserData,MaxKrylovIters,State,vecLength);
		Epi2_FD = new Epi2RetCode(OneD_RHS_FD, UserData, MaxKrylovIters, State, vecLength);
		Epi3V	= new Epi3VChem_KIOPS(OneD_RHS_FD, UserData, MaxKrylovIters, State, vecLength);


		

		PrintPreRun(StepSize, Delx, Steps, KrylovTol, absTol, relTol, Method, num_pts, BAR);
		if(Method == "EPI3V")
			std :: cout << BAR << BAR <<  " EPI3V Selected " << BAR << BAR << std :: endl;
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
			//Set data
			Set_ThermTransData(&problem2, State, sol);//Cantera Change
			problem2.Set_ScalarGradient(State);//Cantera Change
			problem2.Set_TransportGradient(State);//Cantera Change 
			

			if(Method == "EPI2" && UseJac==1)
			{
				OneD_JacArray(TNow, State, StateDot, A, UserData, State, State, State);
				integratorStats =Epi2->Integrate(StepSize, TNow, TNext, NumBands,
					State, KrylovTol, startingBasisSizes);
			}
			else if(Method == "EPI2" && UseJac==0)
			{
				integratorStats =Epi2_FD->Integrate(StepSize, TNow, TNext, NumBands,
					State, KrylovTol, startingBasisSizes);
			}
			else if(Method == "EPI3V")
			{
				//std :: cout << "Attempting EPI3V" << std :: endl;
				integratorStats = Epi3V->NewIntegrate(StepSize, StepSize, KrylovTol, absTol, relTol, TNow, TNext,
					startingBasisSizes, State);
				//std :: cout << "Completed stepping to: " << TNext << std :: endl << std :: endl << std :: endl;
			}
			
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
			}

			//Clean
			Clean(vecLength, StateData);
			Set_ThermTransData(&problem2, State, sol);//Cantera Update
			//Vel Update
			if(VelUp==1)
			{
				problem2.UpdateOneDVelCrossDiff(State);
				problem2.SetVelAve();
			}
			//Check heating
			problem2.CheckHeatingFinal(State, TNext, HeatEnd);

			//Track Progress
			ProgressDots=TrackProgress(FinalTime, TNext, PercentDone, ProgressDots);
			
			if(Movie ==1 ) //Use if we want a time-series plot
			{
				if((StepCount+1)%StepRate==0)
				{
				PrintDataToFile(myfile, StateData,vecLength, absTol, BAR, MyFile, TNext);
				PrintDataToFile(myfile, N_VGetArrayPointer(problem2.Vel), num_pts+1, 0, BAR, MyFile, 0);
				}
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
		LocalPrint(StateData, Experiment, NumScalarPoints, BAR);//Local, in this file.
		PrintProfiling(integratorStats, Profiling, Method,  BAR, cvode_mem);
		PrintDataToFile(myfile, StateData, vecLength, absTol, BAR, MyFile, KTime);//change  4th
		if(Profiling == 1)//Refactor into profiling later.
		{
			cout << "Max Krylov Iterates: " << MaxIters << " at " << SlowTime <<endl;
			cout << "\nTotal integration time: " << KTime << " sec " << endl;
			cout << "RHS calls: " << RHSCnt << "\t";
			cout << "Total time: " << problem2.rhsTime << " sec" << endl;
		}
		PrintDataToFile(myfile, N_VGetArrayPointer(problem2.Vel), num_pts+1, 0, BAR, MyFile, 0);
        myfile.close();
		//==========================
		//Time to take out the trash
		//==========================
		delete Epi2;

		N_VDestroy_Serial(y);
		N_VDestroy_Serial(State);
		N_VDestroy_Serial(StateDot);
		SUNMatDestroy(A);
        //CVodeFree(&cvode_mem);
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

//=================================================
//||      <:::>   :>   ::  <::::>      <::::>    ||
//||     <>   <>  |:>  ::  |::         |::  `>   ||
//||    <>     <> |: > ::  |::::> <::> |::   *>  ||
//||     <>   <>  |:  >::  |::         |::  .>   ||
//||      <:::>   |:   >:  <::::>      <::::>    ||
//=================================================
//=================================================
//  ___    _   _   ___
// |   \  | | | | /   \
// | |) ) | |_| | \ \\/
// |   <  |  _  |  \ \
// | |\	\ | | | | /\\ \
// |_| \_\|_| |_| \___/
//=================================================
//========================================================
//These use the second version of the problem class, myPb2.
//========================================================
//ToDo: Cleanup timing protocols
//Cross Diffusion turned off.
int OneD_RHS(realtype t, N_Vector State, N_Vector StateDot, void * UserData)
{
	//Top Declarations==================================
	myPb2 * pb{static_cast<myPb2 *> (UserData)};//Recast
	int Length 			= N_VGetLength(State);
	realtype * TmpData 	= NV_DATA_S(pb->Tmp);
	realtype * SDD 		= NV_DATA_S(StateDot);
	realtype * GhostD	= NV_DATA_S(pb->Ghost);
	N_VScale(0.0, StateDot, StateDot);
	N_VScale(0.0, pb->Tmp, pb->Tmp);
	//===================================================
	//ChemRHS============================================
	auto StartChem=std::chrono::high_resolution_clock::now();	//Clock
	pb->RHS_Chem(t, State, StateDot, UserData);
	//SUPER_CHEM_RHS_TCHEM(t, State, StateDot, UserData);			//Call RHS_TCHEM onto StateDot
	N_VScale(pb->React, StateDot, StateDot);					//Move Reaction to StateDot
	auto StopChem =std::chrono::high_resolution_clock::now();
	auto PassChem = std::chrono::duration_cast<std::chrono::nanoseconds>(StopChem-StartChem);
	pb->rhs_Chem+=PassChem.count()/1e9;
	pb->rhsTime+=PassChem.count()/1e9;

	//HeatingRHS===========================================
	pb->RHS_Heat(t, State, pb->Tmp, UserData);
	//SUPER_RHS_HEATING(t, State, pb->Tmp, UserData); 				//Heating
	N_VLinearSum(pb->Power, pb->Tmp, 1.0, StateDot, StateDot);		//Add the heating to stateDot

	//DiffRHS================================================
	auto StartDiff=std::chrono::high_resolution_clock::now();	//Start Timing Diff
	pb->RHS_Diff(t, State, pb->Tmp, UserData);
	N_VLinearSum(pb->Diff, pb->Tmp, 1.0, StateDot, StateDot);
	OneD_RHS_CrossDiff(t, State, pb->Tmp, UserData);
	N_VLinearSum(pb->Diff, pb->Tmp, 1.0, StateDot, StateDot);	
	//SUPER_RHS_DIFF_CP(t, State, pb->Tmp, UserData);				//Diff call	
	auto StopDiff=std::chrono::high_resolution_clock::now();
	auto PassDiff = std::chrono::duration_cast<std::chrono::nanoseconds>(StopDiff-StartDiff);
	pb->rhs_Diff+=PassDiff.count()/1e9;							//Finish timing Diff
	pb->rhsTime+=PassDiff.count()/1e9;							//Finish timing Diff
	//N_VLinearSum(pb->Diff, pb->Tmp, 1.0, StateDot, StateDot);	//Add Diff to soln

	//AdvRHS================================================
	auto StartAdv=std::chrono::high_resolution_clock::now();	//Start timing Adv
	pb->RHS_Adv(t,State, pb->Tmp, UserData);
	//SUPER_RHS_ADV_UPW(t,State,pb->Tmp,UserData);				//Centered Adv call into tmp
	auto StopAdv=std::chrono::high_resolution_clock::now();
	auto PassAdv = std::chrono::duration_cast<std::chrono::nanoseconds>(StopAdv-StartAdv);
	pb->rhs_Adv+=PassAdv.count()/1e9;							//Finish timing Adv
	pb->rhsTime+=PassAdv.count()/1e9;							//Finish timing Adv
	N_VLinearSum(pb->Adv, pb->Tmp, 1.0, StateDot, StateDot);	//Add Adv (tmp) to StateDot


	//Finalize==============================================
	N_VScale(0.0, pb->Tmp, pb->Tmp);							//Clean Tmp
	RHSCnt ++;
	return 0;													//Return to caller
}



//==========================================
//  ______  _____  __  __  _____  __     __
// |_    _||   __||  ||  ||  ___||  \   /  |
//   |  |  |  |__ |      ||  ___||   \ /   |
//   |__|  |_____||__||__||_____||__|\_/|__|
//==========================================
int OneD_RHS_Chem(realtype t, N_Vector State, N_Vector StateDot, void * pb)
{
	myPb2 * pbPtr{static_cast<myPb2 *> (pb)};//Recast
	int Length = N_VGetLength(State);
	int num_eqs = pbPtr->num_equations;
	int grid_sz = pbPtr->NumGridPts;
	N_Vector yTemp = N_VNew_Serial(num_eqs);
	N_Vector fTemp = N_VNew_Serial(num_eqs);
	realtype * DATA = NV_DATA_S(yTemp);
	realtype * STATEDATA = NV_DATA_S(State);
	realtype * FDATA = NV_DATA_S(fTemp);
	realtype * STATEDOTDATA = NV_DATA_S(StateDot);

	for(int i = 0 ; i < grid_sz ; i ++ )
	{//March over copies/grid points and grab a data set.
		SUPER_2_VEC(i, DATA, STATEDATA, num_eqs, grid_sz);//Appears to be working
		using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
		using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
		// output rhs vector and x via a kokkos wrapper.
		real_type_1d_view_type x( DATA,     num_eqs);
		real_type_1d_view_type f(FDATA,	    num_eqs);
		//===============================
		// Compute right hand side vector
		//===============================
		auto member =  Tines::HostSerialTeamMember();
		pbPtr->pb.computeFunction(member, x ,f);
		VEC_2_SUPER(i, FDATA, STATEDOTDATA, num_eqs, grid_sz);//Set FDATA into STATEDOTDATA
		//std :: cout << "Step " << i << " successful\n";
	}
	N_VDestroy_Serial(yTemp);
	N_VDestroy_Serial(fTemp);
	return 0;
}

//==================================
//Super Adv
//==================================
// __    __  ___   _
// \ \  / / |  _| | |
//  \ \/ /  |  _| | |__
//   \__/   |___| |____|
//==================================
int OneD_RHS_Adv(realtype t, N_Vector State, N_Vector StateDot, void * userData)
{
	myPb2 * problem{static_cast<myPb2 *> (userData)};//Recast
	realtype divisor 		= -1.0 / (problem->delx);
	realtype * uData   		= N_VGetArrayPointer(State);
	realtype * resultData   = N_VGetArrayPointer(StateDot);
	realtype * Ghost 		= NV_DATA_S(problem->Ghost);
	realtype * VAD			= N_VGetArrayPointer(problem->VelAve);
	int numPts 				= problem->NumGridPts;
	int Grid 				= 0;
	int FI 					= 0;
	double ap 				= 0;				//a^+ for upwind
	double am 				= 0;				//a^- for upwind
	double vp 				= 0;				//v^+ for the scheme
	double vm 				= 0;				//v^- for the scheme
	for (int i = 0; i < problem->vecLength ; i++)
	{
		Grid 	= floor(i/problem->NumGridPts);
		FI  	= i % numPts;						//First grid index
		if(VAD[FI] > 0.0)
			ap = VAD[FI];
		else
			ap = 0.0;
		if(VAD[FI] < 0.0)
			am = VAD[FI];
		else
			am = 0.0;
		//Second order method
		if( FI == 0 ) 							//left boundary
		{
			vp = (-3 * uData[i] + 4 * uData[i+1] - uData[i+2])/2.0;
			vm = ( 3 * uData[i] - 3 * Ghost[Grid] )/2.0;		//uData[i-2] = uData[i-1] = Ghost
		}
		else if( FI == 1)							//Next to left
		{
			vp = (-3 * uData[i] + 4 * uData[i+1] - uData[ i +2 ])/2.0;
			vm = ( 3 * uData[i] - 4 * uData[i-1] + Ghost[ Grid ])/2.0;		//uData[i-2] = Ghost
		}
		else if( FI == (numPts - 2) ) 				//next to right boundary
		{
			vp = (-3 * uData[i] + 3 * uData[i+1])/2.0;					//uData[i+2] = uData[i+1]
			vm = ( 3 * uData[i] - 4 * uData[i-1] + uData[i-2])/2.0;
		}
		else if( FI == (numPts - 1) )				//right boundary
		{
			vp = 0;
			vm = (3*uData[i] -4 * uData[i-1] + uData[i-2])/2.0;
		}
		else
		{
			vp = (-3 * uData[i] + 4 * uData[i+1] - uData[i+2])/2.0;
			vm = (3 * uData[i] -4 * uData[i-1] + uData[i-2])/2.0;
		}
		resultData[i] = divisor * (ap * vm + am * vp);
	}
	return 0;
}

//===============================
//  ____    __   _____   _____
// |  _ \  |  | |   __| |   __|
// | |_) | |  | |   __| |   __|
// |____/  |__| |__|    |__|
//===============================
//Current version uses lookup table
int OneD_RHS_Diff(realtype t, N_Vector u, N_Vector uDot, void * userData)
{
    myPb2 * problem{static_cast<myPb2 *> (userData)};//Recast
    realtype * uData        = NV_DATA_S(u);//Stuff comes from here.
    realtype * resultData   = NV_DATA_S(uDot);//stuff goes in here.
    int numPts 				= problem->NumGridPts;
    realtype divisor        = 1.0/(problem->delx * problem->delx);
    int grid 				= 0;
	int TI					= 0;
    realtype T  			= 1;
    realtype * Ghost 		= NV_DATA_S(problem->Ghost);
    int vecLength 			= problem->num_equations * problem->NumGridPts;
	realtype * 	DiffGridPtr = NV_DATA_S(problem->DiffGrid);
	realtype *  RhoGridPtr	= NV_DATA_S(problem->RhoGrid);
	realtype *  CpGridPtr	= NV_DATA_S(problem->CpGrid); 
	realtype *	MassPtr		= NV_DATA_S(problem->MolarWeights);
	//Start main loop
	for (int i = 0; i < vecLength; i++)
	{
		TI		= i % numPts;
		grid 	= floor(i/problem->NumGridPts);

		if(i < numPts)//If we are looking at temp
		{
			T = DiffGridPtr[TI]/(CpGridPtr[TI]*RhoGridPtr[TI]); 
		}
		else
		{
			T = DiffGridPtr[grid * problem->NumGridPts + TI];
		}
		
		//Set the result
        if(i% numPts == 0)//left
        	resultData[i] = T * divisor * (Ghost[grid] - 2*uData[i] + uData[i+1]);
        else if (i % numPts == (numPts - 1) )//right 0 neumann
            resultData[i] = T * divisor * ( uData[i-1] - uData[i] );
        else
            resultData[i] = T * divisor * (uData[i-1] - 2*uData[i] + uData[i+1]);
    }
    return 0;
}





//=============================
// _______
//|__   __| _
//   | |  _| |_  __    __
// _ | | |_   _| \ \  / /
//\ \| |   | |    \ \/ /
// \___|   |_|     \__/
//=============================
// ______  ______ __   __
//|_    _||_    _|\ \ / /
// _|  |    |  |   \   /
//|____|    |__|    \_/
//=============================
//Cross Diffusion turned off.
int OneD_JtV(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void * pb, N_Vector tmp)
{
	myPb2 * pbPtr{static_cast<myPb2 *> (pb)};			//Recast
	N_VScale(0.0, Jv, Jv);
	N_VScale(0.0, tmp, tmp);
	if(pbPtr->React>0){//Get Chem
		//Chem
		//OneD_Jtv_ChemArray(v, Jv, t, u, fu, pb, tmp);
		pbPtr->JtV_Chem(v, Jv, t, u, fu, pb, tmp);
		N_VScale(pbPtr->React, Jv, Jv);
	}


	auto StartAdv=std::chrono::high_resolution_clock::now();
	
	pbPtr->JtV_Adv(v, tmp, t, u, fu, pb, tmp);
	//OneD_JtV_Adv(v, tmp, t, u, fu, pb, tmp);
	auto StopAdv=std::chrono::high_resolution_clock::now();
	auto PassAdv = std::chrono::duration_cast<std::chrono::nanoseconds>(StopAdv-StartAdv);
	pbPtr->jtv_Adv+=PassAdv.count()/1e9;
	pbPtr->jacTime+=PassAdv.count()/1e9;
	N_VLinearSum(pbPtr->Adv, tmp , 1.0, Jv, Jv);

	// Diff
	auto StartDiff=std::chrono::high_resolution_clock::now();
	pbPtr->JtV_Diff(v,tmp,t,u, fu, pb, tmp);
	N_VLinearSum(pbPtr->Diff, tmp, 1.0, Jv, Jv);
	pbPtr->JtV_CrossDiff(v,tmp, t, u , fu, pb, tmp);
	N_VLinearSum(pbPtr->Diff, tmp, 1.0, Jv, Jv);

	auto StopDiff=std::chrono::high_resolution_clock::now();
	auto PassDiff = std::chrono::duration_cast<std::chrono::nanoseconds>(StopDiff-StartDiff);
	pbPtr->jtv_Diff+=PassDiff.count()/1e9;
	pbPtr->jacTime+=PassDiff.count()/1e9;
	N_VLinearSum(pbPtr->Diff, tmp, 1.0, Jv, Jv);

	//Finalize
	
    return 0;
}



//=================================
//    ___  ____  ___  _____  __
//   |   \/    \|	\|_   _|/  \ 
//   |  _/| () ||	/  | |  \ \/
//   | |  | __ ||	\  | |	/\ \
//   |_|  |_||_||_|\_\ |_|  \__/
//==================================
// __    __  ___   _
// \ \  / / |  _| | |
//  \ \/ /  |  _| | |__
//   \__/   |___| |____|
//==================================
int OneD_JtV_Adv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void* userData, N_Vector tmp)
{
	myPb2 * problem{static_cast<myPb2 *> (userData)};//Recast
	realtype divisor 		= -1.0 / (problem->delx);
	realtype * uData        = NV_DATA_S(u);//Stuff comes from here.
	realtype * vData		= NV_DATA_S(v);//call from here
	realtype * resultData   = NV_DATA_S(Jv);//stuff goes in here.
	realtype * VAD			= N_VGetArrayPointer(problem->VelAve);
	int numPts 				= problem->NumGridPts;
	int Grid 				= 0;
	int FI 					= 0;
	double ap 			= 0;				//a^+ for upwind
	double am 			= 0;				//a^- for upwind
	double vp 			= 0;				//v^+ for the scheme
	double vm 			= 0;				//v^- for the scheme
	for (int i = 0; i < problem->vecLength ; i++)
	{
		Grid 	= floor(i/problem->NumGridPts);
		FI  	= i % numPts;						//First grid index
		//Do max(VAD[i],0)
		if(VAD[FI] > 0.0)
			ap = VAD[FI];
		else
			ap = 0.0;

		//Do min(VAD[i], 0)
		if(VAD[FI] < 0.0)
			am = VAD[FI];
		else
			am = 0.0;
		//Second order.
		if( FI == 0 ) 							//left boundary
		{
			vp = (-3 * vData[i] + 4 * vData[i+1] - vData[i+2])/2.0;
			vm = ( 3 * vData[i] )/2.0;								//vData[i-1] & vData[i-2] drops off
			resultData[i] = divisor * (ap * vm + am * vp);
		}
		else if( FI == 1)							//Next to left
		{
			vp = (-3 * vData[i] + 4 * vData[i+1] - vData[i+2])/2.0;
			vm = ( 3  *vData[i] - 4 * vData[i-1])/2.0;				//vData[i-2] drops off
			resultData[i] = divisor * (ap * vm + am * vp);
		}
		else if( FI == (numPts - 2) ) 				//next to right boundary
		{
			vp = (-3 * vData[i] + 3 * vData[i+1])/2.0;			//vData[i+2] vanishes
			vm = (3*vData[i] -4 * vData[i-1] + vData[i-2])/2.0;
			resultData[i] = divisor * (ap * vm + am * vp);
		}
		else if( FI == (numPts - 1) )				//right boundary
		{
			vp = (3 * vData[i])/2.0;									
			vm = (3 * vData[i] -4 * vData[i-1] + vData[i-2])/2.0;
			resultData[i] = divisor * (ap * vm + am * vp);
		}
		else
		{
			vp = (-3 * vData[i] + 4 * vData[i+1] - vData[i+2])/2.0;
			vm = (3*vData[i] -4 * vData[i-1] + vData[i-2])/2.0;
			resultData[i] = divisor * (ap * vm + am * vp);
		}
	}
	JtvAdv ++;
	return 0;
}

//===============================
//  ____    __   _____   _____
// |  _ \  |  | |   __| |   __|
// | |_) | |  | |   __| |   __|
// |____/  |__| |__|    |__|
//===============================
//=====================
//Diff JTV
//Current version with lookup timing
//=====================
int OneD_JtV_DiffFast(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu,void* userData, N_Vector tmp)
{
	myPb2 * problem{static_cast<myPb2 *> (userData)};//Recast
    realtype * uData        = NV_DATA_S(u);//Stuff comes from here.
    realtype * resultData   = NV_DATA_S(Jv);//stuff goes in here.
    realtype * vData        = NV_DATA_S(v);//Stuff also comes from here
	realtype * CpPtr		= NV_DATA_S(problem->CpGrid);
	realtype * DiffPtr		= NV_DATA_S(problem->DiffGrid);
	realtype * RhoPtr		= NV_DATA_S(problem->RhoGrid);
	realtype * WPtr			= NV_DATA_S(problem->MolarWeights);

    int numPts              = problem->NumGridPts;
    realtype delx           = problem->delx;
    realtype divisor        = 1.0/(delx * delx);
    int TI                  = 0;
    int grid                = 0;
    realtype T              = 0;					//Diff term
    //realtype c              = divisor;
	for (int i = 0; i < problem->vecLength; i++)
	{
		TI   = i%numPts;
		grid = floor(i/numPts);
		//Set Coefficient for step
		if(i < numPts)						//Parse if we are in Temp or not
		{
			//T 		= DiffPtr[ i ];
			T 	= DiffPtr[ i ] / (RhoPtr[ i ] * CpPtr [i]  ); //get lambda/rhoCp for this temp
		}
		else 
		{
			//T			= DiffPtr[TI];///WPtr[grid];
			T 		= DiffPtr[ grid * numPts + TI];								//We need to get DT for species i
		}
		//End parse Thermal lookup
		//Main if
		//std :: cout << T  << std :: endl;
		if(i% numPts == 0)														//left
			resultData[i] = divisor * T * (-2*vData[i] + vData[i+1]);
    	else if (i % numPts == (numPts - 1) )									//right 0 neumann
            resultData[i] = divisor * T * ( vData[i-1] - vData[i] );
        else
        	resultData[i] = divisor * T * (vData[i-1] - 2*vData[i] + vData[i+1]);
    }
	JtvDif++;
    return 0;
}

//==========================================
//  ______  _____  __  __  _____  __     __
// |_    _||   __||  ||  ||  ___||  \   /  |
//   |  |  |  |__ |      ||  ___||   \ /   |
//   |__|  |_____||__||__||_____||__|\_/|__|
//==========================================
//Computes the TChem JtV
//=========================================
//=======================================================
//This Jtv uses the JacArray, an array of small Jacobians
//=======================================================

int OneD_Jtv_ChemArray(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void* pb, N_Vector tmp)
{
	myPb2 * pbPtr{static_cast<myPb2 *> (pb)};//Recast
	auto StartChem=std::chrono::high_resolution_clock::now();
	//Set necessary temp Data
	int num_eqs				= pbPtr->num_equations;
	int num_grid			= pbPtr->NumGridPts;
	//Data that needs local copies
	//realtype * JACDATA		= NV_DATA_S(pbPtr->Jac);
	realtype * JACDATA		= NV_DATA_S(pbPtr->Jacs[0]);
	realtype * VDATA 		= NV_DATA_S(v);
	realtype * JVDATA 		= NV_DATA_S(Jv);
	int Block 				= num_eqs*num_eqs;
	N_Vector SmallJV 		= N_VNew_Serial(num_eqs*num_grid);
	N_VConst(0.0, SmallJV);						//Small fix
	realtype * SmallJVData	= NV_DATA_S(SmallJV);
	N_Vector SmallV			= N_VNew_Serial(num_eqs*num_grid);
	realtype * SmallVData	= NV_DATA_S(SmallV);
	N_Vector SmallTmp 		= N_VNew_Serial(num_eqs*num_grid);

	//#pragma omp parallel for firstprivate(JACDATA, SmallVData, SmallJVData, JVDATA, VDATA, SmallTmp)
	for( int i = 0; i < num_grid; i ++ )
	{//Over each grid
		auto Start=std::chrono::high_resolution_clock::now();
		JACDATA =	NV_DATA_S(pbPtr->Jacs[i]);
		SUPER_2_VEC(i, SmallJVData, JVDATA, num_eqs, num_grid);
		SUPER_2_VEC(i, SmallVData, VDATA, num_eqs, num_grid);
		auto Stop=std::chrono::high_resolution_clock::now();
		auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
		pbPtr->dataMoveTime+=Pass.count()/1e9;

		MatVecProdFast(num_eqs, pbPtr->Jacs[i], SmallV, SmallTmp, SmallJVData);

		//MatrixVectorProduct(num_eqs, JACDATA, SmallV, SmallTmp, SmallJVData);

		auto Start2=std::chrono::high_resolution_clock::now();
		VEC_2_SUPER(i, SmallJVData, JVDATA, num_eqs, num_grid);
		auto Stop2=std::chrono::high_resolution_clock::now();
		auto Pass2 = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop2-Start2);
		pbPtr->dataMoveTime+=Pass2.count()/1e9;
	}
	auto StopChem=std::chrono::high_resolution_clock::now();
	auto PassChem = std::chrono::duration_cast<std::chrono::nanoseconds>(StopChem-StartChem);
	pbPtr->jtv_Chem+=PassChem.count()/1e9;
	pbPtr->jacTime+=PassChem.count()/1e9;
	//N_VDestroy_Serial(SmallJac);
	N_VDestroy_Serial(SmallTmp);
	N_VDestroy_Serial(SmallJV);
	JtvCnt ++;
	JtvChm ++;

	return 0;
}


//===============================
//TCHEM Stuff
//===============================
//Super Jac Via TChem
//===============================
// ______  _____  ______
//|_    _||  _  ||   ___|
// _|  |  |     ||  |___
//|____|  |__|__||______|
//================================
//===========================================
//This method uses an array of smaller Jac's
//No overflow problems as every LittleJac is 54x54
//===========================================
int OneD_JacArray(realtype t, N_Vector State, N_Vector StateDot, SUNMatrix Jac, void * pb, N_Vector tmp1,
				N_Vector tmp2, N_Vector tmp3)
{
	myPb2 * pbPtr{static_cast<myPb2 *> (pb)};//Recast
	auto Start=std::chrono::high_resolution_clock::now();
	int Length 				= N_VGetLength(State);
	int num_eqs 			= pbPtr->num_equations;
	int grid_sz 			= pbPtr->NumGridPts;		//This give # data copies
	N_Vector yTemp 			= N_VNew_Serial(num_eqs);
	realtype * DATA 		= NV_DATA_S(yTemp);
	realtype * STATEDATA	= NV_DATA_S(State);
	realtype * LittleJacPtr	= NULL;
	using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
	using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
	using real_type_2d_view_type = Tines::value_type_2d_view<real_type,host_device_type>;
	int AccessCount	= 0;
	real_type_1d_view_type x;
	real_type_2d_view_type J;
	

	for(int i = 0 ; i < grid_sz ; i ++ )
	{//March over copies/grid points and grab the data needed
		SUPER_2_VEC(i, DATA, STATEDATA, num_eqs, grid_sz);
		LittleJacPtr 	=	NV_DATA_S(pbPtr->Jacs[i]);
		//output Jac vector and x via a kokkos wrapper.
		real_type_1d_view_type x(DATA,     num_eqs);
		real_type_2d_view_type J(LittleJacPtr,     num_eqs, num_eqs);
		//=================================
		// Compute Jacobian
		//=================================
		auto member =  Tines::HostSerialTeamMember();
		pbPtr->pb.computeJacobian(member, x ,J);		//Stores row major
	}
	
	auto Stop=std::chrono::high_resolution_clock::now();
	auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
	pbPtr->jacMakeTime+=Pass.count()/1e9;
	N_VDestroy_Serial(yTemp);
	JacCnt ++;
	return 0;
}

//As of 5/18/23 this function passed the comparison test against NGA.
int OneD_RHS_Heat(realtype t, N_Vector u, N_Vector uDot, void * userData)
{
	myPb2 * problem{static_cast<myPb2 *> (userData)};		//Recast
	realtype * resultData   = NV_DATA_S(uDot);				//stuff goes in here.
	int numPts 				= problem->NumGridPts;
	realtype delx     		= problem->delx;
	realtype * Ghost 		= NV_DATA_S(problem->Ghost);
	int vecLength 			= problem->NumGridPts;			//Only march through temp
	realtype x 				= 0;
	realtype OffSet 		= 0;
	realtype * CpData       = NV_DATA_S(problem->CpGrid);
 realtype * RhoData      = NV_DATA_S(problem->RhoGrid);
 realtype * DiffData     = NV_DATA_S(problem->DiffGrid);
	int TInd				= 0;
	if(problem->HeatingOn==1 )
	{
		OffSet = delx*(round(0.9*problem->NumGridPts) - 0.5);	//Set 1/10 to left
		//OffSet = 0.000895;
        	for (int i = 0; i < numPts; i++) //Only Heat the Temp
		{
			//TInd = problem->TempTableLookUp(uData[i], problem->TempTable);
			x = i*delx + delx/2; 											// x = delx*( i + 0.5)
			// a exp( - 1/(2* rad^2) * (x-center)^2)  !!Need to divide by rho CP at everypoint!!
			resultData[i] = 5e10 * exp( -1e8 * pow(x - OffSet, 2) )/(CpData[i]*RhoData[i]);
		}//Changed below from +delx/2
		problem->HeatingRightGhost=5e10 * exp( -1e8 * pow(x+delx - OffSet, 2) )/(CpData[numPts-1] *RhoData[numPts-1]);
	}
	return 0;
}



//   .oooooo.                                                                                  
//  d8P'  `Y8b                                                                                 
// 888          oooo d8b  .ooooo.   .oooo.o  .oooo.o                                           
// 888          `888""8P d88' `88b d88(  "8 d88(  "8                                           
// 888           888     888   888 `"Y88b.  `"Y88b.                                            
// `88b    ooo   888     888   888 o.  )88b o.  )88b                                           
//  `Y8bood8P'  d888b    `Y8bod8P' 8""888P' 8""888P'                                           
// oooooooooo.    o8o   .o88o.  .o88o.                       o8o                               
// `888'   `Y8b   `"'   888 `"  888 `"                       `"'                               
//  888      888 oooo  o888oo  o888oo  oooo  oooo   .oooo.o oooo   .ooooo.  ooo. .oo.          
//  888      888 `888   888     888    `888  `888  d88(  "8 `888  d88' `88b `888P"Y88b         
//  888      888  888   888     888     888   888  `"Y88b.   888  888   888  888   888         
//  888     d88'  888   888     888     888   888  o.  )88b  888  888   888  888   888         
// o888bood8P'   o888o o888o   o888o    `V88V"V8P' 8""888P' o888o `Y8bod8P' o888o o888o        
// oooooooooooo                                       .    o8o                                 
// `888'     `8                                     .o8    `"'                                 
//  888         oooo  oooo  ooo. .oo.    .ooooo.  .o888oo oooo   .ooooo.  ooo. .oo.    .oooo.o 
//  888oooo8    `888  `888  `888P"Y88b  d88' `"Y8   888   `888  d88' `88b `888P"Y88b  d88(  "8 
//  888    "     888   888   888   888  888         888    888  888   888  888   888  `"Y88b.  
//  888          888   888   888   888  888   .o8   888 .  888  888   888  888   888  o.  )88b 
// o888o         `V88V"V8P' o888o o888o `Y8bod8P'   "888" o888o `Y8bod8P' o888o o888o 8""888P' 
//Jtv Cross Diff
//hold everything but the state constant, so we gradient the v on the second term
//Calculates [1/rho cp * grad(lambda) dot grad(T), 1/rho* grad(rho Di) dot Yi]
int OneD_RHS_CrossDiff(realtype t, N_Vector u, N_Vector uDot, void * userData)
{
    myPb2 * problem{static_cast<myPb2 *> (userData)};	//Recast
    //realtype * uData        = NV_DATA_S(u);//Stuff comes from here.
    realtype * resultData   = NV_DATA_S(uDot);//stuff goes in here.
	//Transport Grids
	realtype * CpPtr		= N_VGetArrayPointer(problem->CpGrid);
	realtype * DiffPtr		= N_VGetArrayPointer(problem->DiffGrid);
	realtype * RhoPtr		= N_VGetArrayPointer(problem->RhoGrid);
	realtype * GradDataPtr	= NV_DATA_S(problem->ScalarGradient);
    int numPts 				= problem->NumGridPts;
    int grid 				= 0;
    //int tempInd 			= 0;
	int TI					= 0;
    realtype T  			= 1;
	int TInd 				= 0;
    realtype * Ghost 		= NV_DATA_S(problem->Ghost);
    int vecLength 			= problem->num_equations * problem->NumGridPts;
	//Transport Gradients
	//realtype * CpGradPtr	= NV_DATA_S(problem->CpGrad);
	realtype * RhoGradPtr	= NV_DATA_S(problem->RhoGrad);
	realtype * DiffGradPtr	= NV_DATA_S(problem->DiffGrad);
	realtype * LambdaGradPtr= NV_DATA_S(problem->lambdaGrad);
	realtype GradData		= 0;
	//Run internal checks.

	for (int i = 0; i < vecLength; i++)
	{

		TI		= i % numPts;
		grid 	= floor(i/problem->NumGridPts);
		GradData= GradDataPtr[i];
		if(i < numPts)//If we are looking at temp
		{
			resultData[i] = 1.0 / (RhoPtr[i] * CpPtr[i] ) * LambdaGradPtr[i] * GradData;
		}
		else
		{	//Use a uniform Dth = DT? No
			resultData[i] = 1.0/ (RhoPtr[TI])	*(RhoGradPtr[TI]*DiffPtr[i] + RhoPtr[TI] * DiffGradPtr[i]) * GradData;
		}
    }
    return 0;
}

int OneD_JtV_CrossDiff(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void* userData, N_Vector tmp)
{
	myPb2 * problem{static_cast<myPb2 *> (userData)};	//Recast
    //realtype * uData        = NV_DATA_S(u);//Stuff comes from here.
	realtype * vData		= NV_DATA_S(v);
    realtype * resultData   = NV_DATA_S(Jv);//stuff goes in here.
	//Transport Grids
	realtype * CpPtr		= NV_DATA_S(problem->CpGrid);
	realtype * DiffPtr		= NV_DATA_S(problem->DiffGrid);
	realtype * RhoPtr		= NV_DATA_S(problem->RhoGrid);
	realtype * GradDataPtr	= NV_DATA_S(problem->ScalarGradient);
    int numPts 				= problem->NumGridPts;
    int grid 				= 0;
    int tempInd 			= 0;
	int TI					= 0;
    //realtype T  			= 1;
	int TInd 				= 0;
    //realtype * Ghost 		= NV_DATA_S(problem->Ghost);
    int vecLength 			= problem->num_equations * problem->NumGridPts;
	//Transport Gradients
	realtype * CpGradPtr	= NV_DATA_S(problem->CpGrad);
	realtype * RhoGradPtr	= NV_DATA_S(problem->RhoGrad);
	realtype * DiffGradPtr	= NV_DATA_S(problem->DiffGrad);		//Contains Grad[lambda,Diff]
	//realtype * LambdaGradPtr= NV_DATA_S(problem->lambdaGrad);
	realtype GradData		= 0;
	realtype denom			= 1.0/(2*problem->delx);
	//Start main loop
	for (int i = 0; i < vecLength; i++)
	{

		TI		= i % numPts;
		grid 	= floor(i/problem->NumGridPts);
		if( TI ==0 ) //Left: fixed
		{
			GradData= vData[i+1]*denom;
		}
		else if (TI == numPts-1)//Right: 0 neumann
		{
			GradData= (vData[i]-vData[i-1])*denom;
		}
		else //mid
		{
			GradData= (vData[i+1]-vData[i-1])*denom;

		}


		if(i < numPts)//If we are looking at temp
		{
			resultData[i] = 1.0 / (RhoPtr[i] * CpPtr[i] ) * DiffGradPtr[i] * GradData;
		}
		else
		{
			resultData[i] = 1.0/ (RhoPtr[TI])	*(RhoGradPtr[TI]*DiffPtr[i] + RhoPtr[TI] * DiffGradPtr[i]) * GradData;			
		}
    }

	return 0;
}

//Initial Conditions
int LocalInitCons(realtype * data, int Experiment, int Sample)
{
	if (Experiment==1) //Hydrogen
	{
		if(Sample==2)
		{
			data[0] 	= 4.4579829e+02;	//Temp
			data[1]		= 0.0;
			data[2] 	= 2.0130322e-02;	//H2
			data[4]		= 2.2822068e-01;	//O2 N2 
			data[7]		= 7.5164900e-01;	//N2
		}
		else if(Sample==1)
		{
			 data[0]=900.0;				//Temp
			data[2]=2.7431213550e-01; 	//H2
			data[4]=1.0-data[2];		//O2
		}
	}
	else if(Experiment==5) //n-butane
	{
		if(Sample==1)
		{
			data[0]         = 1200.0;
			data[11]        = 2.173895119224689421e-01;  //O2
			data[155]       = 7.092943418325631244e-01;  //N2
			data[154]       = 1.256572768198739760e-02;  //AR
			data[58]        = 6.075041856298054460e-02;  //C4H10
		}
		else if(Sample==2)
		{
			data[0]			= 4.4579829e+02; //Temp
			data[11]		= 2.2275981e-01; //O2
			data[155] 		= 7.3366351e-01; //N2
			data[58]		= 4.3576684e-02; //C4H10
		}
	}
	return 0;
}

int LocalPrint(realtype * data, int Experiment, int num_pts, string bar)
{
	cout << bar << "Printing SuperVector" << bar << endl;
	for(int i = 0; i < num_pts; i++)
	{
		if(Experiment==1){//Temp H2 O2 O OH H2O H HO2 H2O2
			std :: cout << "Temp=" << data[i] << "\t\t H2=" << data[2*num_pts+i];
			std :: cout <<"\t\t O2=" << data[4*num_pts + i];
		}else if (Experiment==2){//Gri
			std :: cout << "Temp=" << data[i] << "\t\t CH4=" << data[14*num_pts + i];
			std :: cout <<"\t\t O2=" << data[4*num_pts + i];
		}else if (Experiment==5){//n-butane
			std :: cout << "Temp=" << data[i] << "\t\t C4H10=" << data[58*num_pts + i];
			std :: cout <<"\t\t O2=" << data[11*num_pts + i];
		}else{
			return 1;
		}
		std :: cout << std :: endl;
	}
	return 0;
	
}



//New stuff
int Set_ThermTransData(myPb2* pbPtr, N_Vector State, std::shared_ptr<Cantera::Solution> sol)
{
	//Declare ytemp and pointers
	N_Vector yTemp 		= N_VNew_Serial(pbPtr->num_equations);
	N_Vector Coeff 		= N_VClone(yTemp);
	realtype * yPtr 	= NV_DATA_S(yTemp);
	realtype * DATA		= NV_DATA_S(State);
	realtype * CpPtr	= NV_DATA_S(pbPtr->CpGrid);
	realtype * RhoPtr	= NV_DATA_S(pbPtr->RhoGrid);
	realtype * DiffPtr 	= NV_DATA_S(pbPtr->DiffGrid);  //Stores [Lambda, Diffs]
	realtype * GWPtr	= NV_DATA_S(pbPtr->GasWeight); 
	auto gas 			= sol->thermo();
	auto kin 			= sol->kinetics();
	auto trans 			= sol->transport();
	//Over this loop get the little Y.
	for( int i = 0 ; i < pbPtr->NumGridPts ; i ++)
	{
		SUPER_2_VEC(i, yPtr, DATA, pbPtr->num_equations, pbPtr->NumGridPts);
		gas->setState_TPY(yPtr[0], pbPtr->pb._p, yPtr+1);//Reconfigure the gas.
		kin = sol->kinetics();
		trans = sol->transport();
		trans->getMixDiffCoeffs(NV_DATA_S(Coeff)+1); //Get Diffusion Coefficients
		//Set Cp & Rho into Grids.
		RhoPtr[i] = sol->thermo()->density();
		CpPtr[i] = sol->thermo()->cp_mass();
		//Set GW
		GWPtr[i]= gas->meanMolecularWeight();
		//Set [Lambda,Diff] grid 
		DiffPtr[i] = trans->thermalConductivity();
		for(int j = 1 ; j < pbPtr->num_equations; j++)
		{
			DiffPtr[i + j*pbPtr->NumGridPts] = NV_DATA_S(Coeff)[j];
		}

	}
	//std::cout << gas->report() << std::endl;
	//Set the boundary information.


	N_VDestroy_Serial(yTemp);
	N_VDestroy_Serial(Coeff);
	return 0;
}

int Set_Derivatives(myPb2* pbPtr, N_Vector State, std::shared_ptr<Cantera::Solution> sol)
{
	N_Vector yTemp 		= N_VNew_Serial(pbPtr->num_equations);
	N_Vector Coeff 		= N_VClone(yTemp);
	realtype * yPtr 	= NV_DATA_S(pbPtr->Ghost);
	realtype * DATA		= NV_DATA_S(State);
	realtype * CpPtr	= NV_DATA_S(pbPtr->CpGrad);
	realtype * RhoPtr	= NV_DATA_S(pbPtr->RhoGrad);
	realtype * DiffPtr 	= NV_DATA_S(pbPtr->DiffGrad);  //Stores [Lambda, Diffs]
	realtype  Delx		= pbPtr->delx;
	//Sets with messed up end points
	pbPtr->Set_TransportGradient(State);
	//Manually fix the left end points. RhoGrad, CpGrad, LambdaGrad, DiffGrad.
	//Set a point using inlet information.
	
	auto gas 			= sol->thermo();
	gas->setState_TPY(yPtr[0], pbPtr->pb._p, yPtr+1);//Reconfigure the gas based on the intial point
	auto kin 			= sol->kinetics();
	auto trans 			= sol->transport();
	trans->getMixDiffCoeffs(NV_DATA_S(Coeff)+1); //Get Diffusion Coefficients

	CpPtr[0] 			= 1.0/(2.0*Delx)* (NV_DATA_S(pbPtr->RhoGrid)[1] - sol->thermo()->density() );  
	RhoPtr[0]			= 1.0/(2.0*Delx)* (NV_DATA_S(pbPtr->CpGrid)[1] - sol->thermo()->cp_mass() ); 

	return 0;	
}

int OneD_RHS_FD(realtype t, N_Vector State, N_Vector StateDot, void * UserData)
{
	myPb2 * pbPtr{static_cast<myPb2 *> (UserData)};//Recast	
	Set_ThermTransData(pbPtr, State, pbPtr->sol);
	pbPtr->Set_TransportGradient(State);
	OneD_RHS(t, State, StateDot, UserData);
	return 0;
}