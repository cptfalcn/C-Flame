/*
 * ===========================================================================================
 * This is serial implementation.
 * ===========================================================================================
 * Console call is as follows:
 * ./CheckCP.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --MyFile="Scrap.txt"
 * --SampleNum=1 --Experiment=2
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
#include "Epi3SC_KIOPS.h"


#include <chrono>
#include "InitialConditions.h"
#include "Chemistry.h"
#include "Print.h"
#include "CreateIntegrators.h"//This includes the EPI3_KIOPS.cpp/h definitions

#define BAR "===================="


//=====================
//Prototypes
//=====================
//Banner
void PrintBanner();

//Misc functions
void SetState(realtype* data, realtype * stateData, double pressure, int length);

int TrackProgress(realtype, realtype, realtype, int);

void CheckNanBlowUp(N_Vector, int);

void PrintData(N_Vector vec);

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
	int ProgressDots	= 0; //From 0 to 100; only care about percentage

	realtype PercentDone	= 0;

	string MyFile		= "Default.txt";

	int SampleNum		= 0;
	int Experiment		= 1; 		//default is hydrogen, set to 0 for kapila

	int number_of_equations	= 0;
	int NumScalarPoints 	= 1;
	int multi 		= 1;		//Added as an input

	//=====================================================
	//Kokkos Input Parser
	//=====================================================
	std::string chemFile("chem.inp");
 	std::string thermFile("therm.dat");
 	// parse command line arguments --chemfile=user-chemfile.inp --thermfile=user-thermfile.dat

 	TChem::CommandLineParser opts("This example computes reaction rates with a given state vector");

 	opts.set_option<std::string>("chemfile", "Chem file name e.g., chem.inp", &chemFile);

 	opts.set_option<std::string>("thermfile", "Therm file name e.g., therm.dat", &thermFile);

	opts.set_option<std::string>("MyFile","Where we output data", &MyFile);//New

	opts.set_option<int>("SampleNum", "Sample used", &SampleNum);

	opts.set_option<int>("Experiment", "Experiment chosen", &Experiment);

	opts.set_option<int>("Multi", "grid multiplier", &multi);

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

		int vecLength 	= number_of_equations;			//For readability
		int num_eqs 	= number_of_equations;			//To reduce size
		int num_pts 	= 1;					//To reduce size
		//==================================================
		//Prep/Set initial conditions, inherited from Zero-D
		//==================================================
		N_Vector y 		= N_VNew_Serial(vecLength);	//intial conditions/Data
		N_Vector yDot		= N_VClone(y);			//State RHS
		N_VScale(0.0 , y, y);					//Zero out the y vector.
		N_VScale(0.0, yDot, yDot);				//Zero out the StateDot vector.
		realtype *data 		= NV_DATA_S(y);			//Set the y data pointer.
		realtype *dataDot 	= NV_DATA_S(yDot);		//Set the State data pointer.
		N_Vector State		= N_VNew_Serial(vecLength+1);	//Need to add pressure for TChem
		realtype *SDot		= N_VGetArrayPointer(State);	//Pointer to State Data.

		//Set up a double length problem
		/**/
		N_Vector y2		= N_VNew_Serial(vecLength*multi);
		N_Vector y2Dot		= N_VClone(y2);
		N_VScale(0.0, y2, y2);
		N_VScale(0.0, y2Dot, y2Dot);
		realtype *data2		= N_VGetArrayPointer(y2);
		realtype *data2Dot	= N_VGetArrayPointer(y2Dot);
		/**/
		cout << BAR << "Initialization complete" << BAR << endl;
		//Initial Conditions sub-block
		SetIntCons(Experiment, SampleNum, data);		//Set Initial Conditions
		SetSuperInitCons(data, data2, num_eqs, num_pts*multi);	//Copy IC to State.
		//===================================================
		//Set TChem
		//===================================================
    		//TChem does not allocate any workspace internally. You create the work space using NVector,
		//std::vector or real_type_1d_view_type (Kokkos view). Here we use kokkos view.
		//Set up an ignition problem and the workspace. Unused in Experiment 0.
		//=========================================================================
		const int problem_workspace_size = problem_type::getWorkSpaceSize(kmcd);
		real_type_1d_view_type work("workspace", problem_workspace_size);

		myPb2 problem(num_eqs, work, kmcd, num_pts, y, 0);	//Construct new Problem.
		problem.SetGhostPVel(y, Experiment, SampleNum, 0);	//Do additional set up.
		problem.SetAdvDiffReacPow(0, 0, 0, 0, 0);//Additional set up.
		problem.kmd		= kmd;
		void *UserData 		= &problem;

		/**/
		myPb2 problem2(num_eqs, work, kmcd, multi, y, 0);
		problem2.SetGhostPVel(y,Experiment, SampleNum,0);
		problem2.kmd		= kmd;
		void *UserData2		= &problem2;
		/**/

		SetState(data, SDot, problem.pb._p, vecLength);
		GetOmegaCP(0, State, problem.CP, UserData);			//Store cp in problem2.CP


		cout << BAR <<"Printing data to "<< MyFile << BAR << endl;
		cout << "Data:\n";
		PrintData(State);
		//CPMass
		cout << "cp:\n";
		PrintData(problem.CP);
		//CPMass Mixture
		cout << "cp mass mix:\n";
		PrintData(problem.SmallScrap);


		//Likely bugged, remove for now
		/*
		cout << "Manual cp mass mix:\n";
		realtype cpmm = problem.ComputeCpTemp(y);
		cout << cpmm << endl;
		*/


		//Compare against the wrapped function
		//Known now to be bugged
		/*
		GetOmegaCP(0, y, problem.CP, UserData);
		cout << "CP with Chem function call\n";
		PrintData(problem.CP);
		cout << "CP mass mix with function call\n";
		PrintData(problem2.SmallScrap);
		*/

		/*
		problem.FetchTherm(y);
		cout << "\nFetch function version\n";
		cout << "cp mass mix:\n";
		PrintData(problem.SmallScrap);
		*/

		/**/
		problem2.FetchTherm(y2);
		//cout << "Lifted y\n";
		//PrintData(y2);
		cout << "\nFetch function on lifted version\n";
		PrintData(problem2.CP);
		cout << "cp mass mix:\n";
		PrintData(problem2.SmallScrap);
		/**/


		myfile.close();
		N_VDestroy(y);
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


void PrintBanner()
{
	//cout << "\n\n\n";
	//cout << "\t================================================\n";
	//cout << "\t||      <:::>   :>   ::  <::::>      <::::>    ||\n";
	//cout << "\t||     <>   <>  |:>  ::  |::         |::  `>   ||\n";
	//cout << "\t||    <>     <> |: > ::  |::::> <::> |::   *>  ||\n";
	//cout << "\t||     <>   <>  |:  >::  |::         |::  .>   ||\n";
	//cout << "\t||      <:::>   |:   >:  <::::>      <::::>    ||\n";
	//cout << "\t================================================\n";

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

void SetState(realtype * data, realtype * stateData, double pressure, int length)
{

	stateData[0] = data[0];		//Temperature is unmoved
	stateData[1] = pressure;	//pressure goes here
	for ( int i = 2; i < length; i ++)
	{
		stateData[i+1] = data[i];
	}

}

void PrintData(N_Vector vec)
{
	realtype * data 	= N_VGetArrayPointer(vec);
	int length		= N_VGetLength(vec);
	for( int i = 0 ;  i < length; i++)
	{
		cout << data[i] << "\t\t";
		if(i%6==5)
			cout << "\n";
	}
	cout << "\n";

}
