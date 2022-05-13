/*
 * ===========================================================================================
 * This is serial implementation.
 * ===========================================================================================
 */
//This implementation is soon to be depreciated
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
void ErrorCheck(ofstream &, N_Vector, realtype *, int, int, realtype);
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
	PrintBanner();//Print One-D, how fancy
	//====================
	// Intial Declarations
	//====================
	int SlowDown 		= 0,	OldProjections	= 0;
	int MaxIters	= 0,	OldIters	= 0,	CurrStepIters	 = 0, TotalIters = 0;
	realtype SlowTime	= 0;
	static realtype FinalTime = 0;//1.0e-4;//1e-4 seems to be the max
	static const int NumBands = 3;//Epic stuff, default is 3.
	realtype StepSize=0;
	realtype KTime=0;
	int ProgressDots=0; //From 0 to 100; only care about percentage
	int StepCount = 0;
	//int SpeedRatio = 1e-2;
	realtype PercentDone		= 0;
	realtype TNow			= 0;
	realtype TNext			= 0;
	realtype TNextC			= 0;			//CVODE TNext, has a tendency to change TNext
	string MyFile 			= "Default.txt";//Default Data dump file
	string JacFile			= "Defaut.txt";	//Default Jac dump file
	string Method 			= "EPI2";		//Default method signature.
	static realtype KrylovTol 	= 1e-14;
	int UseJac				= 1; 			//Use the Jacobian?  1/0
	int SampleNum			= 0;
	int Experiment			= 1; 			//default is hydrogen
	int Profiling			= 0;			//default to no profiling.
	int number_of_equations	= 0;
	int NumScalarPoints 	= 0;
	realtype absTol 		= 1e-8;
	realtype relTol 		= 1e-8;
	int startingBasisSizes[]= {10, 10};//{3,3}
	int Movie	= 0;
	//=====================================================
	//Kokkos Input Parser
	//=====================================================
	std::string chemFile("chem.inp");
 	std::string thermFile("therm.dat");
 	/// parse command line arguments --chemfile=user-chemfile.inp --thermfile=user-thermfile.dat
	/// --Stepsize=1e-8  --FinalTime=1e-2  --MyFile="Filename"  --KrylovTole=1e-14
	//  --UseJac=0 --SampleNum={1,2,3}, --Method="Integrator name"
	//  --Experiment={0,1,2}, Profiling = {0,1}, Delx = SpatialDisc
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
	opts.set_option<int>("NumPts", "Interior points for each grid", &NumScalarPoints);
	opts.set_option<int>("Movie", "Generate a data set at every step", &Movie);
	opts.set_option<std::string>("JacFile", "The dump file for the Jacobian", &JacFile);

	const bool r_parse = opts.parse(argc, argv);
	if (r_parse)
		return 0; // print help return

	ofstream myfile(MyFile, std::ios_base::app);
	ofstream jacfile(JacFile, std::ios_base::app);
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

		realtype Delx 		= 0;					//Forced to 0
		NumScalarPoints		= 1;					//Fixed for 0-D
		int vecLength 		= number_of_equations;	//For readability
		int num_eqs 		= number_of_equations;	//To reduce size
		int num_pts 		= 1;					//To reduce size
		//==================================================
		//Prep/Set initial conditions
		//==================================================
		N_Vector y 			= N_VNew_Serial(num_eqs);	//intial conditions.
		N_VScale(0.0 , y, y);							//Zero out the y vector.
		realtype *data 		= NV_DATA_S(y);				//Set the y data pointer.
		SetIntCons(Experiment, SampleNum, data);		//Set Initial Conditions
		//===================================================
		//Set TChem
		//===================================================
		//TChem does not allocate any workspace internally. You create the work space using NVector,
		//std::vector or real_type_1d_view_type (Kokkos view). Here we use kokkos view.
		//=========================================================================
		const int problem_workspace_size = problem_type::getWorkSpaceSize(kmcd);
		real_type_1d_view_type work("workspace", problem_workspace_size);

		myPb2 problem2(num_eqs, work, kmcd, num_pts, y, Delx);	//Construct new Problem.
		std :: cout << problem2.NumGridPts << std :: endl;
		problem2.SetAdvDiffReacPow(0, 0, 1, 0, 0);//Additional set up.
		problem2.kmd 		= kmd;
		void *UserData 		= &problem2;
		int MaxKrylovIters 	= 500; //max(vecLength, 500);//500 is the base
		//==================
		//Create integrators
		//==================
		//CVODE
		void * cvode_mem;
		SUNMatrix A							= SUNDenseMatrix(vecLength, vecLength);
		SUNLinearSolver SUPERLS 			= SUNLinSol_SPGMR(y, PREC_NONE, 20);
		SUNNonlinearSolver SUPERNLS 		= SUNNonlinSol_Newton(y);
		
		//Set EPIXX_KIOPS methods
		Epi2_KIOPS	*Epi2					= NULL;
		Epi3_KIOPS 	*Epi3					= NULL;
		EpiRK4SC_KIOPS	*EpiRK4				= NULL;
		IntegratorStats *integratorStats 	= NULL;
		//ToDo: Figure out why the code is acting strange when using "Super" versions
		Epi2 = 	new Epi2_KIOPS(CHEM_RHS_TCHEM, CHEM_JTV_V2, UserData,MaxKrylovIters, y,vecLength);

		Epi3 = 	new Epi3_KIOPS(CHEM_RHS_TCHEM_V2, CHEM_JTV_V2,UserData,MaxKrylovIters, y,vecLength);

		//EpiRK4=	new EpiRK4SC_KIOPS(SUPER_RHS,SUPER_JTV, UserData, MaxKrylovIters, y,vecLength);

		cvode_mem= CreateCVODE(CHEM_RHS_TCHEM_V2, CHEM_JTV_V2, SUPER_CHEM_JAC_TCHEM, UserData, A,
					SUPERLS, SUPERNLS, vecLength, y, relTol, absTol, StepSize, 1);
		CVodeSetMaxNumSteps     (cvode_mem, 1e6);

		PrintPreRun(StepSize, Delx, Steps, KrylovTol, absTol, relTol, Method, num_pts, BAR);
		//=================================
		// Run the time integrator loop
		//=================================
		cout <<"[";
		cout.flush();
		while(StepCount<Steps)//while(TNext<FinalTime)
		{
			TNow= StepCount*StepSize;
			TNext=TNow + StepSize;
			TNextC = TNext;
			//Integrate
			auto Start=std::chrono::high_resolution_clock::now();//Time integrator
			if(Method != "CVODEKry")//Set the step Jacobian
				CHEM_COMP_JAC_V2(y, UserData); //State
			if(Method == "EPI2")
				integratorStats =Epi2->Integrate(StepSize, TNow, TNext, NumBands,
					y, KrylovTol, startingBasisSizes);
			else if(Method == "EPI3")
				integratorStats =Epi3->Integrate(StepSize, TNow, TNext, NumBands,
					y, KrylovTol, startingBasisSizes);
			else if(Method == "EPIRK4")
				integratorStats =EpiRK4->Integrate(StepSize, TNow, TNext, NumBands,
					y, KrylovTol, startingBasisSizes);
			else if(Method == "CVODEKry")//This appears to be bugged for the full problem
				CVode(cvode_mem, TNext, y, &TNextC, CV_NORMAL);//CV_NORMAL/CV_ONE_STEP
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
			if(SampleNum!=10) 		//Error check
				ErrorCheck(myfile, y, data, num_eqs, num_pts, TNext);//Needs editing
			CheckNanBlowUp(y, vecLength);

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
			Clean(vecLength, data);

			//Track Progress
			ProgressDots=TrackProgress(FinalTime, TNext, PercentDone, ProgressDots);
			// if( data[0] >1500)
			// 	std :: cout << "Ign Delay acheived" << std :: endl;

			if(Movie ==1 ) //Use if we want a time-series plot
			{
				PrintDataToFile(myfile, data,vecLength, absTol, BAR, MyFile, TNext);
				PrintDataToFile(myfile, N_VGetArrayPointer(problem2.Vel), num_pts+1,
						0, BAR, MyFile, 0);
				realtype * jacPtr       = N_VGetArrayPointer(problem2.Jac);
                        	for(int i = 0 ; i < N_VGetLength(problem2.Jac); i ++)
								jacfile << jacPtr[i] << "\n";
                        	jacfile << StepSize << "\n";
			}
			StepCount++;
		}//End integration loop
        TNow=TNext;
        cout << "]100%\n" << BAR << "\tIntegration complete\t" << BAR <<endl;
		cout << "Congrats\n";
		
	    //=======================================
        //Console Output
        //=======================================
		cout << BAR <<"Printing data to "<< MyFile << BAR << endl;
		PrintExpParam(FinalTime, TNow, StepSize, StepCount, KrylovTol, absTol, relTol, KTime, BAR);
		PrintSuperVector(data, Experiment, 1, BAR);
		PrintProfiling(integratorStats, Profiling, Method,  BAR, cvode_mem);
		PrintDataToFile(myfile, data, vecLength, absTol, BAR, MyFile, KTime);//change  4th

		//Refactor into PrintProfiling later.
		if(Profiling == 1)
		{
			if(Method == "EPI2")
				cout << "Max Kry Iters: " << MaxIters << " at " << SlowTime << endl;
			if(Method == "CVODEKry")
			{
				cout << "Max Krylov Iterates: " << MaxIters << " at " << SlowTime <<endl;
				cout << "Total Krylov iterates: " << TotalIters<< endl;
			}
		}

        myfile.close();
		jacfile.close();
		//==========================
		//Time to take out the trash
		//==========================
		delete Epi2;
		delete Epi3;
		//delete EpiRK4;
		N_VDestroy_Serial(y);
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


//=======================================================
//RHS Error Check
//realtype* y 		y Data ptr
//realtype* rhs		rhs Data ptr
//int num_eqs		number of equations in the system
//========================================================
void ErrorCheck(ofstream & myfile, N_Vector y, realtype * data, int number_of_equations, int num_pts,
			realtype TNext)
{
	realtype MassError=abs( N_VL1NormLocal(y)/num_pts -data[0] -1.0);
	if(MassError>.1 || abs(data[0])>1e5 )
	{
		//PrintDataToFile(myfile, data, number_of_equations, TNext);
                cout<<"\n===============!!!Critical error!!!================\n";
                if(MassError>.1)//originally .1
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
