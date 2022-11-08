/*
 * ===========================================================================================
 * This is serial implementation.
 * ===========================================================================================
 */
//ToDo:  Move all the RHS/Jac functions to this file.
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
int OneD_RHS			(realtype, N_Vector, N_Vector, void *);
int OneD_RHS_Adv		(realtype, N_Vector, N_Vector, void *);
int OneD_RHS_Chem		(realtype, N_Vector, N_Vector, void *);
int OneD_RHS_Diff		(realtype, N_Vector, N_Vector, void *);
int OneD_RHS_CrossDiff	(realtype, N_Vector, N_Vector, void *);
int OneD_RHS_Heat		(realtype, N_Vector, N_Vector, void *);
//Split RHS
int OneD_RHS_First		(realtype, N_Vector, N_Vector, void *);	//Adv-Diff
int OneD_RHS_Second		(realtype, N_Vector, N_Vector, void *); //Kinetics-Heating

//Create Jacobians
int OneD_Jac			(realtype, N_Vector, N_Vector, SUNMatrix, void *, N_Vector,
						N_Vector, N_Vector);
int OneD_JacArray		(realtype, N_Vector, N_Vector, SUNMatrix, void *, N_Vector,
						N_Vector, N_Vector);

int OneD_JtV			(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
int OneD_JtV_Adv		(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
int OneD_JtV_Diff		(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
int OneD_JtV_DiffFast	(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);

int OneD_JtV_CrossDiff	(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
int OneD_Jtv_Chem		(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);

int OneD_Jtv_ChemArray	(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
//Split Jac Adv Diff
int OneD_JtV_First		(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);

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
	string InitialData		= "InitialData.txt";	// ""

	realtype ADV 			= 1.0;			//Fixed defaults for RHS
	realtype DIFF 			= 1.0;
	realtype CHEM 			= 1.0;
	realtype POW			= 0.0;

	//omp_set_thread_num(4);
	//std :: cout << "max omp threads: " << omp_get_max_threads() << std :: endl;
	//std :: cout << "omp threads: " << omp_get_num_threads() << std :: endl;

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
	opts.set_option<int>("VelUp", "Do we do velocity up?", &VelUp);			//In final version remove this to make it fixed
	opts.set_option<int>("NumPts", "Interior points for each grid", &NumScalarPoints);
	opts.set_option<int>("Movie", "Generate a data set at every step", &Movie);
	opts.set_option<int>("Restart", "Restart from text file?", &Restart);

	const bool r_parse = opts.parse(argc, argv);
	if (r_parse)
		return 0; // print help return

	ofstream myfile(MyFile, std::ios_base::app);
	ofstream VelFile("VelDiv.txt", std::ios_base::app);
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
		// N_VScale(0.0, LittleJac, LittleJac);			//Zero out LittleJac.

		realtype *data 		= NV_DATA_S(y);				//Set the y data pointer.
		realtype *StateData = NV_DATA_S(State);			//Set the State data pointer.

		//Initial Conditions sub-block
		SetIntCons(Experiment, SampleNum, data);		//Set Initial Conditions
		SetSuperInitCons(data, StateData, num_eqs, num_pts);//Copy IC to State.
		TestingInitCons(SampleNum, num_pts, num_eqs, vecLength, Delx, data, StateData);//Keep until final version
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
	

		//read in data, error if files do not exist
		//These tables extend into a higher temp range, Dth=DT
		// ReadData(problem2.CPPoly,		"BisettiCp.txt");
		// ReadData(problem2.TempTable,	"BisettiTemp.txt");
		// ReadData(problem2.RhoTable, 	"BisettiRho.txt");
		// ReadData(problem2.DiffTable,	"BisettiNewDiff.txt");

		ReadData(problem2.CPPoly,		"Cp_fT.txt");
		ReadData(problem2.TempTable,	"Tf.txt");
		ReadData(problem2.RhoTable, 	"rho_fT.txt");
		ReadData(problem2.DiffTable,	"Diff_fT.txt");
		//Need to modify the first entries of the DiffTable to fix my error.
		realtype * CpPtr 	= NV_DATA_S(problem2.CPPoly);
		realtype * RhoPtr 	= NV_DATA_S(problem2.RhoTable); 
		realtype * DiffPtr	= NV_DATA_S(problem2.DiffTable);
		// for( int i = 0 ; i < 500 ; i++)
		// {
		// 	DiffPtr[i]= DiffPtr[i]*(RhoPtr[i]*CpPtr[i]);//Generates Lambda.
		// }
		//DiffPtr[1:500] are lambdas for fixed temperatures.
		//problem2.VerifyTempTable(State);
		//Create Cantera
		auto sol = Cantera::newSolution("gri3.0/gri30.yaml", "gri30", "None");
		

		//End Cantera
		void *UserData = &problem2;
		int MaxKrylovIters = max(vecLength, 500);//500 is the base
		//==================
		//Create integrators
		//==================
		//CVODE
		void * cvode_mem;
		void * cvode_split;
		cvode_split = CVodeCreate(CV_BDF);
		int retVal=0;
		SUNMatrix A							= SUNDenseMatrix(vecLength, vecLength);
		SUNLinearSolver SUPERLS 			= SUNLinSol_SPGMR(State, PREC_NONE, 20);
		SUNNonlinearSolver SUPERNLS 		= SUNNonlinSol_Newton(State);

		
		//Set EPI_KIOPS methods
		Epi2_KIOPS	*Epi2					= NULL;
		IntegratorStats *integratorStats 	= NULL;
		IntegratorStats *intStatsSplit2		= NULL;

		//Split EPI2
		Epi2_KIOPS	*Epi2_AdvDiff			= NULL;
		Epi2_KIOPS	*Epi2_ChemHeat			= NULL;
		//===================================================
		//Parse the experiment cases and make the integrators
		//===================================================
		Epi2 = 	new Epi2_KIOPS(OneD_RHS,OneD_JtV,UserData,MaxKrylovIters,State,vecLength);
		
		//Segfaults if not done this way.  However, we skip Split if done this way.
		if(Method == "CVODEKry")
		{
			cvode_mem = CreateCVODE(OneD_RHS, OneD_JtV, OneD_JacArray, UserData, A,
					SUPERLS, SUPERNLS, vecLength, State, relTol, absTol, StepSize, 1);
		}
		if(Method == "CVODESplit")
		{
			//cout << "Setting split method\n";
			cvode_mem = CreateCVODE(OneD_RHS_First, OneD_JtV, OneD_JacArray, UserData, A,
				SUPERLS, SUPERNLS, vecLength, State, relTol, absTol, StepSize, 1);
			N_Vector AbsTol= N_VNew_Serial(num_eqs);
			for ( int i = 0 ; i < num_eqs ; i++)
				NV_Ith_S(AbsTol,i)=absTol;

			cvode_mem = CVodeCreate (CV_BDF);
			CVodeSetUserData(cvode_split, UserData);
			CVodeInit(cvode_split, OneD_RHS_First, 0, State);
			CVodeSetInitStep(cvode_split, StepSize);
			CVodeSetMaxStep(cvode_split, StepSize);
			CVodeSVtolerances(cvode_split, relTol, AbsTol);
			CVodeSetLinearSolver(cvode_split, SUPERLS, A);//A might be null, try that
			CVodeSetNonlinearSolver(cvode_split, SUPERNLS);
			if (UseJac==1)
			{
				CVodeSetJacTimes(cvode_split, NULL, OneD_JtV_First);
				CVodeSetLSetupFrequency(cvode_split, 1); //Remove because also stupid.
				CVodeSetJacEvalFrequency(cvode_split, 1);//Remove becuase this is stupid.
				CVodeSetJacFn(cvode_split, OneD_JacArray);		//Error if removed .
				CVodeSetMaxNumSteps(cvode_split, 500);
			}
			N_VDestroy_Serial(AbsTol);
		}
		

		Epi2_AdvDiff = new Epi2_KIOPS(OneD_RHS_First, OneD_JtV_First, UserData,
										MaxKrylovIters,State,vecLength);
		
		Epi2_ChemHeat = new Epi2_KIOPS(OneD_RHS_Second, OneD_Jtv_ChemArray, UserData,
										MaxKrylovIters,State,vecLength);

		PrintPreRun(StepSize, Delx, Steps, KrylovTol, absTol, relTol, Method, num_pts, BAR);

		//====================
		//Testing block
		//====================
		// problem2.Test_ScalarGradient(problem2.Scrap);
		// cout << "Constant 1 test run\n\n";
		// problem2.Test_ScalarGradient(State);
		//problem2.Test_TransportGradient(State);
		//problem2.Test_RHS_CrossDiff(State);
		//problem2.Test_JtV_CrossDiff(State);
		// cout << "State test run\n";
		
		//problem2.Test_GasWeight(State);


		// problem2.Test_VelocityDivergence(State);
		// exit(EXIT_FAILURE);
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
			problem2.Set_ScalarGradient(State);
			problem2.SetTransportGrid(State);
			problem2.Set_TransportGradient(State);
			

			if(Method == "EPI2")
			{
				OneD_JacArray(TNow, State, StateDot, A, UserData, State, State, State);
				integratorStats =Epi2->Integrate(StepSize, TNow, TNext, NumBands,
					State, KrylovTol, startingBasisSizes);
			}
			else if(Method == "EPI2Split")
			{
				OneD_JacArray(TNow, State, StateDot, A, UserData, State, State, State);
				integratorStats =Epi2_AdvDiff->Integrate(StepSize, TNow, TNext, NumBands,
					State, KrylovTol, startingBasisSizes);

				OneD_JacArray(TNow, State, StateDot, A, UserData, State, State, State);
				intStatsSplit2 =Epi2_ChemHeat->Integrate(StepSize, TNow, TNext, NumBands,
					State, KrylovTol, startingBasisSizes);

			}
			else if(Method =="CVODESplit")
			{
				//cout << "Attemping CVODE\n";
				// if(StepCount ==0)
				// 	{
				// 		OneD_JacArray(TNow, State, StateDot, A, UserData, State, State, State);
				// 		integratorStats =Epi2_AdvDiff->Integrate(StepSize, TNow, TNext, NumBands,
				// 			State, KrylovTol, startingBasisSizes);
				// 	}
				// else
					CVode(cvode_split, TNext, State, &TNextC, CV_ONE_STEP);//CV_NORMAL/CV_ONE_STEP
				cout << "CVODE Step done\n";
				OneD_JacArray(TNow, State, StateDot, A, UserData, State, State, State);
				//cout << "Attempting EPI2\n";
				integratorStats =Epi2_ChemHeat->Integrate(StepSize, TNow, TNext, NumBands,
					State, KrylovTol, startingBasisSizes);
			}
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
			{
				problem2.SetTransportGrid(State);
				//problem2.UpdateOneDVel(State);
				problem2.Set_VelocityDivergence(State);
				problem2.Print_MaterialDerivative(problem2.VelAve, VelFile);
				problem2.VelIntegrate(NV_DATA_S(problem2.VelAve), State, problem2.VelAveLeftBnd, problem2.VelAveRightBnd);
				problem2.SetVelAve();
				
				problem2.Set_MaterialDerivative(StepSize, problem2.MatDerivative);
				problem2.Print_MaterialDerivative(problem2.VelScrap, VelFile);								
			}
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
		if(Method == "EPI2Split")
		{
			std :: cout << "\nSplit second stats\n\n";
			intStatsSplit2->PrintStats();
		}
		PrintDataToFile(myfile, StateData, vecLength, absTol, BAR, MyFile, KTime);//change  4th
		if(Profiling == 1)//Refactor into profiling later.
		{
			cout << "Max Krylov Iterates: " << MaxIters << " at " << SlowTime <<endl;
			cout << "\nTotal integration time: " << KTime << " sec " << endl;
			cout << "RHS calls: " << RHSCnt << "\t";
			cout << "Total time: " << problem2.rhsTime << " sec" << endl;
			//cout << "\t Diff rhs: " <<problem2.rhs_Diff << " sec\n";
			//cout << "\t\t LookupTime: " << problem2.rhsDiffLookTime << " sec\n";
			//cout << "\t Adv rhs: " << problem2.rhs_Adv << " sec\n";
			//cout << "\t Chem rhs: " << problem2.rhs_Chem << " sec\n";
			cout << "JtV calls: " << JtvCnt << "\t";
			cout << "Total time: " << problem2.jacTime << " sec" << endl;
			cout << "\t Diff jtv: " <<problem2.jtv_Diff << " sec\n";
			cout << "\t num calls:" << JtvDif << endl;
			//cout << "\t\t LookupTime: " << problem2.jtvDiffLookTime << " sec\n";
			cout << "\t Adv jtv: " << problem2.jtv_Adv << " sec\n";
			cout << "\t num calls: " << JtvAdv << endl;
        	cout << "\t Chem jtv: " << problem2.jtv_Chem << " sec\n";
			cout << "\t Chem calls: " << JtvChm << endl;
			//cout << "\t\t Move time: " << problem2.dataMoveTime << " sec\n";
			cout << "Jac calls: " << JacCnt << "\t";
			cout << "Total time: " << problem2.jacMakeTime << " sec\n";
		}
		PrintDataToFile(myfile, N_VGetArrayPointer(problem2.Vel), num_pts+1, 0, BAR, MyFile, 0);
        myfile.close();
		//==========================
		//Time to take out the trash
		//==========================
		delete Epi2;

		delete Epi2_AdvDiff;
		delete Epi2_ChemHeat;
		N_VDestroy_Serial(y);
		N_VDestroy_Serial(State);
		N_VDestroy_Serial(StateDot);
		SUNMatDestroy(A);
		SUNLinSolFree(SUPERLS);
		SUNNonlinSolFree(SUPERNLS);
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

	//Given the temperature, set all the grids
	// pb->SetTransportGrid(State);
	// pb->Set_ScalarGradient(State);
	// pb->Set_TransportGradient(State);

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

//=================================
//    __   ___   _     _   _____
//   /  \ |   \ | |   | | |_   _|
//   \ \/ |  _/ | |   | |   | |
//   /\ \ | |   | |_  | |   | |
//   \__/ |_|   |___| |_|   |_|
//================================================
//The first RHS of the split method uses ADV-Diff
//================================================
int OneD_RHS_First		(realtype t, N_Vector State, N_Vector StateDot, void * UserData)
{
	myPb2 * pb{static_cast<myPb2 *> (UserData)};//Recast
	N_VScale(0.0, StateDot, StateDot);
	N_VScale(0.0, pb->Tmp, pb->Tmp);
	//Given the temperature, set all the grids
	auto Start=std::chrono::high_resolution_clock::now();
	auto StartDiff=std::chrono::high_resolution_clock::now();	//Start Timing Diff
	pb->RHS_Diff(t, State, pb->Tmp, UserData);
	N_VLinearSum(pb->Diff, pb->Tmp, 1.0, StateDot, StateDot);
	//OneD_RHS_CrossDiff(t, State, pb->Tmp, UserData);
	//N_VLinearSum(pb->Diff, pb->Tmp, 1.0, StateDot, StateDot);
	//SUPER_RHS_DIFF_CP(t, State, pb->Tmp, UserData);				//Diff call	
	auto StopDiff=std::chrono::high_resolution_clock::now();
	auto PassDiff = std::chrono::duration_cast<std::chrono::nanoseconds>(StopDiff-StartDiff);
	pb->rhs_Diff+=PassDiff.count()/1e9;							//Finish timing Diff
	//N_VLinearSum(pb->Diff, pb->Tmp, 1.0, StateDot, StateDot);	//Add Diff to soln


	auto StartAdv=std::chrono::high_resolution_clock::now();	//Start timing Adv
	pb->RHS_Adv(t,State,pb->Tmp, UserData);
	//SUPER_RHS_ADV_UPW(t,State,pb->Tmp,UserData);				//Centered Adv call into tmp
	auto StopAdv=std::chrono::high_resolution_clock::now();
	auto PassAdv = std::chrono::duration_cast<std::chrono::nanoseconds>(StopAdv-StartAdv);
	pb->rhs_Adv+=PassAdv.count()/1e9;							//Finish timing Adv
	N_VLinearSum(pb->Adv, pb->Tmp, 1.0, StateDot, StateDot);	//Add Adv (tmp) to StateDot

	N_VScale(0.0, pb->Tmp, pb->Tmp);							//Clean Tmp


	auto Stop=std::chrono::high_resolution_clock::now();
    auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
    pb->rhsTime+=Pass.count()/1e9;									//Final timing set
	RHSCnt ++;
	return 0;

}

//==========================
//Kinetics and Heating
//==========================
int OneD_RHS_Second		(realtype t, N_Vector State, N_Vector StateDot, void * UserData)
{
	myPb2 * pb{static_cast<myPb2 *> (UserData)};//Recast
	N_VScale(0.0, StateDot, StateDot);
	N_VScale(0.0, pb->Tmp, pb->Tmp);

	auto Start=std::chrono::high_resolution_clock::now();
	

	auto StartChem=std::chrono::high_resolution_clock::now();	//Clock
	pb->RHS_Chem(t, State, StateDot, UserData);
	//SUPER_CHEM_RHS_TCHEM(t, State, StateDot, UserData);			//Call RHS_TCHEM onto StateDot
	//N_VScale(pb->React, StateDot, pb->Scrap);					//Save the reaction in Scrap
	N_VScale(pb->React, StateDot, StateDot);					//Move Reaction to StateDot
	auto StopChem =std::chrono::high_resolution_clock::now();
	auto PassChem = std::chrono::duration_cast<std::chrono::nanoseconds>(StopChem-StartChem);
	pb->rhs_Chem+=PassChem.count()/1e9;
	pb->RHS_Heat(t,State,pb->Tmp, UserData);
	//SUPER_RHS_HEATING(t, State, pb->Tmp, UserData); 				//Heating
	N_VLinearSum(pb->Power, pb->Tmp, 1.0, StateDot, StateDot);		//Add the heating to stateDot

	N_VScale(0.0, pb->Tmp, pb->Tmp);							//Clean Tmp


	auto Stop=std::chrono::high_resolution_clock::now();
    auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
    pb->rhsTime+=Pass.count()/1e9;									//Final timing set
	RHSCnt ++;
	return 0;														//Return to caller


}


//=================================
//    ___  ____  ___  _____  __
//   |   \/    \|	\|_   _|/  \ 
//   |  _/| () ||	/  | |  \ \/
//   | |  | __ ||	\  | |	/\ \
//   |_|  |_||_||_|\_\ |_|  \__/

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
	//realtype * LookupTemp	= N_VGetArrayPointer(problem->TempTable);
	//realtype * LookupCp		= N_VGetArrayPointer(problem->CPPoly);
	//realtype * LookupDiff	= N_VGetArrayPointer(problem->DiffTable);
	//realtype * LookupRho	= N_VGetArrayPointer(problem->RhoTable);
    int numPts 				= problem->NumGridPts;
    realtype divisor        = 1.0/(problem->delx * problem->delx);
    int grid 				= 0;
    //int tempInd 			= 0;
	int TI					= 0;
    realtype T  			= 1;
	//int TInd 				= 0;
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

		//No clue why I would still be using this.
		//auto Start=std::chrono::high_resolution_clock::now();
		//TInd = problem->TempTableLookUp(uData[TI], problem->TempTable);
		//auto Stop=std::chrono::high_resolution_clock::now();
		//auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
		//problem->rhsDiffLookTime+=Pass.count()/1e9;


		if(i < numPts)//If we are looking at temp
		{
			T = DiffGridPtr[TI];//(CpGridPtr[TI]*RhoGridPtr[TI]); 
			//T = LookupDiff[ TInd ] / (LookupRho[ TInd ] * LookupCp [ TInd]  );
		}
		else
		{
			T = DiffGridPtr[TI];///MassPtr[grid];//Uniform DT=Dth
			//T = DiffGridPtr[grid * problem->NumGridPts + TI];
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
	//pbPtr->SetTransportGrid(u);
	//pbPtr->Set_ScalarGradient(u);
	//pbPtr->Set_TransportGradient(u);


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
	//OneD_JtV_Diff(v,tmp,t,u, fu, pb, tmp);
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
//    __   ___   _     _   _____
//   /  \ |  _\ | |   | | |_   _|
//   \ \/ |  _/ | |   | |   | |
//   /\ \ | |   | |_  | |   | |
//   \__/ |_|   |___| |_|   |_|
//==================================
//Split Part 1: Adv-Diff, then call Chem Jtv Directly
//============================
int OneD_JtV_First(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void * pb, N_Vector tmp)
{
	myPb2 * pbPtr{static_cast<myPb2 *> (pb)};			//Recast
	N_VScale(0.0, Jv, Jv);
	N_VScale(0.0, tmp, tmp);
	pbPtr->SetTransportGrid(u);

	auto Start=std::chrono::high_resolution_clock::now();

	//Adv
	auto StartAdv=std::chrono::high_resolution_clock::now();
	OneD_JtV_Adv(v, tmp, t, u, fu, pb, tmp);
	auto StopAdv=std::chrono::high_resolution_clock::now();
	auto PassAdv = std::chrono::duration_cast<std::chrono::nanoseconds>(StopAdv-StartAdv);
	pbPtr->jtv_Adv+=PassAdv.count()/1e9;
	N_VLinearSum(pbPtr->Adv, tmp , 1.0, Jv, Jv);
	
	// Diff
	auto StartDiff=std::chrono::high_resolution_clock::now();
	pbPtr->JtV_Diff(v,tmp,t,u, fu, pb, tmp);
	//OneD_JtV_Diff(v,tmp,t,u, fu, pb, tmp);
	N_VLinearSum(pbPtr->Diff, tmp, 1.0, Jv, Jv);
	//pbPtr->JtV_CrossDiff(v,tmp, t, u , fu, pb, tmp);
	//N_VLinearSum(pbPtr->Diff, tmp, 1.0, Jv, Jv);
	auto StopDiff=std::chrono::high_resolution_clock::now();
	auto PassDiff = std::chrono::duration_cast<std::chrono::nanoseconds>(StopDiff-StartDiff);
	pbPtr->jtv_Diff+=PassDiff.count()/1e9;
	N_VLinearSum(pbPtr->Diff, tmp, 1.0, Jv, Jv);

	
	//Finish total timing
	auto Stop=std::chrono::high_resolution_clock::now();
    auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
    pbPtr->jacTime+=Pass.count()/1e9;
	JtvDif ++;
	JtvAdv ++;
	JtvCnt ++;
	//cout << "Exiting JtV_First\n";
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
			T 		= DiffPtr[ i ];
			//T 	= DiffPtr[ i ] / (RhoPtr[ i ] * CpPtr [i]  ); //get lambda/rhoCp for this temp
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
	
	//member Kokkos::Impl::HostThreadTeamMember<Kokkos::Serial>;
	//omp_set_dynamic(0);     // Explicitly disable dynamic teams
	//omp_set_num_threads(8); // Use 4 threads for all consecutive parallel regions
	//#pragma omp parallel for private(DATA, STATEDATA, LittleJacPtr, x, J, pbPtr)
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





//==================================
//  __  __  ____  ______  ______  
// |  ||  ||  __||  __  ||_    _|
// |      ||  __||      |  |  |
// |__||__||____||__||__|  |__|
//==================================
//COMPLETED AND DEBUGGED
//==================================
int OneD_RHS_Heat(realtype t, N_Vector u, N_Vector uDot, void * userData)
{
	myPb2 * problem{static_cast<myPb2 *> (userData)};		//Recast
	realtype * uData        = NV_DATA_S(u);					//Stuff comes from here.
	realtype * resultData   = NV_DATA_S(uDot);				//stuff goes in here.
	N_VScale(0.0, uDot, uDot);								//Clean, new
	int numPts 				= problem->NumGridPts;
	realtype delx           = problem->delx;
	realtype * Ghost 		= NV_DATA_S(problem->Ghost);
	int vecLength 			= problem->NumGridPts;			//Only march through temp
	realtype x 				= 0;
	realtype OffSet 		= 0;
	realtype * LookupTemp	= N_VGetArrayPointer(problem->TempTable);
	realtype * LookupCp		= N_VGetArrayPointer(problem->CPPoly);
	realtype * LookupDiff	= N_VGetArrayPointer(problem->DiffTable);
	realtype * LookupRho	= N_VGetArrayPointer(problem->RhoTable);
	int TInd				= 0;

	if(problem->HeatingOn==1 )
	{
		OffSet = delx*(round(0.9*problem->NumGridPts) - 0.5);	//Set 1/10 to left
		//OffSet = 0.000895;
        for (int i = 0; i < numPts; i++) //Only Heat the Temp
		{
			TInd = problem->TempTableLookUp(uData[i], problem->TempTable);
			x = i*delx + delx/2; 											// x = delx*( i + 0.5)
			// a exp( - 1/(2* rad^2) * (x-center)^2)  !!Need to divide by rho CP at everypoint!!
			resultData[i] = 5e10 * exp( -1e8 * pow(x - OffSet, 2) )/(LookupCp[TInd]*LookupRho[TInd]);
		}//Changed below from +delx/2
		problem->HeatingRightGhost=5e10 * exp( -1e8 * pow(x+delx - OffSet, 2) )/(LookupCp[TInd]*LookupRho[TInd]);
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
		{	//Use a uniform Dth = DT
			resultData[i] = 1.0/ (RhoPtr[TI])	*(RhoGradPtr[TI]*DiffPtr[TI] + RhoPtr[TI] * DiffGradPtr[TI]) * GradData;
		}
    }
    return 0;
}

int OneD_JtV_CrossDiff(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void* userData, N_Vector tmp)
{
	myPb2 * problem{static_cast<myPb2 *> (userData)};	//Recast
    realtype * uData        = NV_DATA_S(u);//Stuff comes from here.
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
    realtype T  			= 1;
	int TInd 				= 0;
    realtype * Ghost 		= NV_DATA_S(problem->Ghost);
    int vecLength 			= problem->num_equations * problem->NumGridPts;
	//Transport Gradients
	realtype * CpGradPtr	= NV_DATA_S(problem->CpGrad);
	realtype * RhoGradPtr	= NV_DATA_S(problem->RhoGrad);
	realtype * DiffGradPtr	= NV_DATA_S(problem->DiffGrad);
	realtype * LambdaGradPtr= NV_DATA_S(problem->lambdaGrad);
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
		//std :: cout <<  "index: " << i << "V[i]= " << vData[i] << "V Gradient: " << GradData << std :: endl;
		//GradData= GradDataPtr[i];
		if(i < numPts)//If we are looking at temp
		{
			resultData[i] = 1.0 / (RhoPtr[i] * CpPtr[i] ) * LambdaGradPtr[i] * GradData;
		}
		else
		{
			resultData[i] = 1.0/ (RhoPtr[TI])	*(RhoGradPtr[TI]*DiffPtr[i] + RhoPtr[TI] * DiffGradPtr[i]) * GradData;			
		}
    }

	return 0;
}
