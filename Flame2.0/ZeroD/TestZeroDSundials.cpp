/* =====================================================================================
New Sundials version
===================================================================================== */
// #include "TChem_CommandLineParser.hpp"

// #include "TChem_Util.hpp"
// #include "TChem_KineticModelData.hpp"
// #include "TChem_Impl_IgnitionZeroD_Problem.hpp" // here is where Ignition Zero D problem is implemented
#include <sundials/sundials_config.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_math.h>
#include <nvector/nvector_serial.h>
//CVODE includes
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sunlinsol/sunlinsol_spgmr.h>  /* SPGMR SUNLinearSolver                */
#include <sunlinsol/sunlinsol_spbcgs.h>       /* access to SPBCGS SUNLinearSolver            */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver         */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <chrono>

using namespace std;
#define NEQ 3
#define EPS 1e-2

N_Vector InitialConditions();
int RHS(realtype, N_Vector, N_Vector, void *);
int Jtv(N_Vector, N_Vector, realtype, N_Vector, N_Vector , void *, N_Vector);
int CheckStep(realtype, realtype);
void postProcess(N_Vector);



int main(int argc, char* argv[])
{
  sundials::Context sunctx;   // SUNDIALS context

  static const realtype FinalTime = 1e1;
	static const int NumBands = 3;
	string MyFile = "OutPutFile.txt";
	realtype StepSize=1e-2;
	int ProgressDots=0; //From 0 to 100; only care about percentage
	int StepCount = 0;
	realtype PercentDone=0;
	void *userData=nullptr;
	realtype TNow=0;
	realtype TNext=0;
	N_Vector y = N_VNew_Serial(NEQ, sunctx); //The data y
	N_VScale(1.0, InitialConditions(), y);
	realtype *data = NV_DATA_S(y);
	string opt = "EPI2" ;
	bool LastStepAdjusted = 0;
	realtype KTime = 0;
	realtype PrintStepSize = 0;
	realtype KrylovTol = 1e-14;
  realtype absTol    = 1e-10;
  realtype relTol    = 1e-6;
	//Input Parsing
	if(argc == 1)
	{
		//CheckStep(FinalTime, StepSize);
	}
	else if(argc== 2)
	{
		StepSize=stof(argv[1]);
		std :: cout << "Step-size: " << StepSize << endl;
		//CheckStep(FinalTime, StepSize);
		opt = "EPI2";
	}
	else if(argc == 3)
	{
		StepSize=stof(argv[1]);
		std :: cout << "Step-size: " << StepSize << endl;
		opt = argv[2];
		cout << "Method: " << opt << endl;
	}
	else if(argc == 4 )
	{
		StepSize = stof(argv[1]);
		std :: cout << "Step-size: " << StepSize << endl;
		opt = argv[2];
		cout << "Method: " << opt << endl;
		MyFile = argv[3];
		cout << "Output file: " << MyFile << endl;
	}
	else if(argc == 5)
	{
		StepSize = stof(argv[1]);
		opt = argv[2];
		MyFile = argv[3];
		KrylovTol = stof(argv[4]);
	}
	else
	{
		cout << "!!!Failure:  Incorrect number of args!!!" << endl;
		exit(EXIT_FAILURE);
	}
	PrintStepSize=StepSize;
	//===================================================
	// Create the cvode integrator
  //===================================================
	void * cvode_mem;
  void * UserData;
  SUNMatrix A						    = SUNDenseMatrix(NEQ, NEQ, sunctx);
  SUNLinearSolver LS 				= SUNLinSol_SPGMR(y, PREC_NONE, 20, sunctx);
  SUNNonlinearSolver NLS 		= SUNNonlinSol_Newton(y,sunctx);
	int retVal = 0;
	realtype tret=0;
	N_Vector AbsTol= N_VNew_Serial(NEQ, sunctx);
	for ( int i = 0 ; i < NEQ ; i++)
    NV_Ith_S(AbsTol,i)=absTol;

	cvode_mem = CVodeCreate (CV_BDF, sunctx);
	//Give CVODE the user problem
	retVal = CVodeSetUserData(cvode_mem, UserData);
	//Set intial state
	retVal = CVodeInit(cvode_mem, RHS, 0, y);
	//Set intial stepsize guess
	retVal = CVodeSetInitStep(cvode_mem, StepSize);
	//Set max stepsize
  retVal = CVodeSetMaxStep(cvode_mem, 1);
  //retVal = CVodeSetMaxStep(cvode_mem, StepSize);
	//retVal = CVodeSetMinStep(cvode_mem, StepSize);
	//Set tolerances
	retVal = CVodeSVtolerances(cvode_mem, relTol, AbsTol);
	//Set linear solver
	retVal = CVodeSetLinearSolver(cvode_mem, LS, A);
	retVal = CVodeSetNonlinearSolver(cvode_mem, NLS);
	retVal = CVodeSetJacTimes(cvode_mem, NULL, Jtv);
  retVal = CVodeSetMaxNumSteps(cvode_mem, max(2000,NEQ*NEQ*NEQ));
	N_VDestroy_Serial(AbsTol);
  
  
  
  const int MaxKrylovIters = 500;
  //EpiRK4SV *integrator = new EpiRK4SV(RHS,
	//EpiRK4SC_KIOPS *integrator = new EpiRK4SC_KIOPS(RHS,
	//EpiRK5C *integrator = new EpiRK5C(RHS,
	// IntegratorStats *Stats = NULL;
  //   	Epi2_KIOPS *integrator 	= new Epi2_KIOPS(RHS, Jtv, userData, MaxKrylovIters, y, NEQ);
	// Epi3_KIOPS * Epi3 	= new Epi3_KIOPS(RHS, Jtv, userData, MaxKrylovIters, y, NEQ);
	// Epi3SC_KIOPS * Epi3SC  	= new Epi3SC_KIOPS(RHS, Jtv, userData, MaxKrylovIters, y, NEQ);
	// //EpiP2_KIOPS * EpiP2 	= new EpiP2_KIOPS(RHS,Jtv, userData, MaxKrylovIters, y, NEQ, 3, 0, &postProcess);
	//========================
	// Set integrator parameters
	//========================
	//const realtype KrylovTol = RCONST(1.0e-14);//1e-14
	int startingBasis[] = {3, 3};
	cout<<"=======================Initial Data======================\n";
	cout << setprecision(17);
	cout <<"y=" << data[0] << "\t\t z=" << data[1] <<"\t\t Temp=" << data[2]<<endl;

  //Run the cvode integrator
  auto Start=std::chrono::high_resolution_clock::now();
  CVode(cvode_mem, FinalTime, y, &TNext, CV_NORMAL);
  auto Stop=std::chrono::high_resolution_clock::now();
  auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
  KTime+=Pass.count()/1e9;


	//=================================
    	// Run the time integrators outside of the loop
	//=================================
	/**/
	// if(opt == "EPI3SC")
	// {
	// 	auto Start=std::chrono::high_resolution_clock::now();//Time integrator
	// 	Stats = Epi3SC->Integrate(StepSize, StepSize*100, 1e-5, 1e-8, 0, FinalTime,
	// 					NumBands, startingBasis, y);
	// 	auto Stop=std::chrono::high_resolution_clock::now();
  //               auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
  //               KTime+=Pass.count()/1e9;
	// 	TNow = FinalTime;
	// }
	// else if(opt == "EPIP2")
	// {
	// 	// auto Start=std::chrono::high_resolution_clock::now();
	// 	// Stats= EpiP2->Integrate(StepSize, 0.0, FinalTime, y, KrylovTol, startingBasis);
	// 	// auto Stop=std::chrono::high_resolution_clock::now();
	// 	// auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
	// 	// KTime+=Pass.count()/1e9;
	// 	// TNow = FinalTime;
	// }

	
	
	// TNext = StepSize;
	// if(TNext>FinalTime)
	// 	StepSize=FinalTime;


	// cout<<"[";
	// while(TNow<FinalTime)
	// {
	// 	//Integrate
	// 	if(opt=="EPI2")
	// 	{
	// 		auto Start=std::chrono::high_resolution_clock::now();
	// 		Stats= integrator->Integrate(StepSize,TNow, TNext, NumBands, y,
	// 					KrylovTol, startingBasis);
	// 		auto Stop=std::chrono::high_resolution_clock::now();
	// 					auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
	// 					KTime+=Pass.count()/1e9;

	// 	}
	// 	else if(opt=="EPI3")
	// 	{
	// 		auto Start=std::chrono::high_resolution_clock::now();
	// 		Stats= Epi3->Integrate(StepSize, TNow, TNext, NumBands, y,
	// 					KrylovTol, startingBasis);
	// 		auto Stop=std::chrono::high_resolution_clock::now();
	// 					auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
	// 					KTime+=Pass.count()/1e9;
	// 	}
	// 	else if(opt=="EPI3SC"){
	// 		break;
	// 	}
	// 	else if(opt == "EPIP2")
	// 	{
	// 		break;
	// 	}
	// 	else
	// 	{
	// 		std :: cout <<"!!!Error:  Invalid Integrator selected\n";
	// 		exit(EXIT_FAILURE);
	// 	}

	// 	//Removing this cleaning loop will cause stalls.
	// 	postProcess(y);
		
	// 	// for (int j=0; j<NEQ; j++)
	// 	// 	if(data[j]<0)
	// 	// 		data[j]=0;
		
	// 	//Track the progress
	// 	PercentDone=floor(TNext/FinalTime*100);
	// 	for (int k=0; k<PercentDone-ProgressDots;k++)
	// 	{
	// 		cout<<":";
	// 		cout.flush();
	// 	}
	// 	ProgressDots=PercentDone;

	// 	//Set next iterate information
	// 	TNow 	+= StepSize;
	// 	TNext 	=  TNow+StepSize;
	// 	if(TNext>FinalTime)
	// 	{
	// 		StepSize=FinalTime-TNow;
	// 		LastStepAdjusted=1;
	// 		TNext = TNow+StepSize;
	// 	}
	// 	StepCount ++;
	// 	//End loop
	// }
	// cout << "]100%\n\n";
	/**/

	


	// printf("Run stats:\n");
  //   	Stats->PrintStats();
  //Runs a quick error check
	realtype error_inf = 0;
  realtype error_tmp = 0;
	
	if(EPS == 1e-1)
	{
		error_inf = max( 
			data[0] - 0.000000002988028,
			max(
				data[1] - 0.032702315935000,
				data[2] - 1.967297681076978
			));
      error_tmp = data[2] - 1.967297681076978;
	}
	else if(EPS == 1e-2)
	{
		error_inf = max( 
			data[0] - 0.000000000000000,
			max(
				data[1] - 0.005848132117459,
				data[2] - 1.994151867882557
			));
    error_tmp = data[2] - 1.994151867882557;
	}

	//=======================================
	//Console Output
	//=======================================
	cout<<"===========================Data========================\n";
	cout << setprecision(17);
	cout <<"y=" << data[0] << "\t\t z=" << data[1] <<"\t\t Temp=" << data[2]<<endl;
	if (LastStepAdjusted == 1)
		cout<<"Altered final ";
	cout << "Simulation ended at time: " << setprecision(16) << TNext << endl;
	cout << "Steps taken: " << StepCount <<endl;
	cout << "Integration time: " <<KTime << endl;
	cout << "StepSize: " << StepSize << endl;

	cout << "Temp Error:" << error_tmp << std::endl;
	cout << "Error Inf:" << error_inf << std::endl;
	//cout <<"Sum of y=" <<data[0]+ data [1] + data [2]<<"\t\t";
	//cout <<"Error in sum=" << data[0]+data[1]+data[2]-Y0[0]-Y0[1]-Y0[2]<<endl;
	cout <<"===================Printing data to " <<MyFile <<"==============\n";

	// error_inf = math.max( 
	// 	 data[0] - 0.000000002988028
	// 	 math.max(
	// 		 data[1] - 0.032702315935000,
	// 		 data[2] - 1.967297681076978
	// 	 )
	// );


//    0.000000002988028
//    0.032702315935000
//    1.967297681076978


	//=======================================
	//Data file output
	//=======================================
	ofstream myfile(MyFile, std::ios_base::app);
	myfile << setprecision(17) << fixed << data[0] << "\t\t" << data[1] << "\t\t";
	myfile << data[2] << "\t\t" << PrintStepSize << "\t\t" << KTime;
	myfile << "\t\t" << KrylovTol;
	//myfile << "\t\t" << Stats->numTimeSteps << endl;
	myfile.close();
	cout << "======================Simulation complete=============\n";

  // Clean up the integrator
  // delete integrator;
	// delete Epi3SC;
	// //delete EpiP2;



}



/*
//=========================================================
//  //====  ==   ==  ==     ==    _.==_    ========  _===_
//  ||      ||   ||  ||\\   ||  ./     \\     ||    //   \\
//  ||====  ||   ||  || \\  ||  ||            ||    \\___
//  ||      ||   ||  ||  \\ ||	~\            ||         \\
//  ||      \\___//  ||   \\||   ~=___-//     ||    \\___//
//=========================================================
*/


//===========================================
//Check the number of steps
//FinalTime
//StepSize
//===========================================
int CheckStep(realtype FinalTime, realtype StepSize)
{
	cout<<"======================Steps=========================\n";
	cout<<"Checking the proposed number of steps...\n";
	cout<<std::setprecision(17);
	cout<<FinalTime/StepSize << " steps proposed...";
  if(floor( FinalTime/StepSize ) == FinalTime/StepSize )
	{
    static const int Steps= FinalTime/StepSize;
    cout << " accepted...\n";
		return Steps;
  }else if (abs(round(FinalTime/StepSize))-FinalTime/StepSize <1e-3 )
  {
	  	static const int Steps= round (FinalTime/StepSize);
	  	cout << Steps << " steps has been approximated...\n";
	  	return Steps;
	}else
  {
    cout<<"Cannot perform non-integer number of steps!!\n";
    exit(1);
  }
}

/*
 * ===========================================================================================
 * 
 * Function InitialConditions
 * 
 * Computes initial condition. 
 * 
 * Output:
 * y0   data
 * 
 * ===========================================================================================
 */

N_Vector InitialConditions()
{
  sundials::Context sunctx;   // SUNDIALS context
  N_Vector y0 = N_VNew_Serial(NEQ, sunctx);
  realtype *data = NV_DATA_S(y0);
	data[0]=1.0-.5*EPS;
	data[1]=.5*EPS;
	data[2]=1.0;
  return y0;
}


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
 */

int RHS(realtype t, N_Vector u, N_Vector udot, void *userData)
{
	realtype *y = NV_DATA_S(u);
	realtype *dy = NV_DATA_S(udot);
	double Omega= .5 * y[0] * y[1] * exp ( (y[2]-1.0 ) /  (EPS * y[2]  ) );
	dy[0] = -Omega;
	dy[1] = Omega-y[1];
	dy[2] = y[1];
//	cout << dy[0] <<"\t\t" << dy[1] << "\t\t" << dy[2] <<endl;
        return 0;
}

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
 * 
 * ===========================================================================================
 */
int Jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void *userData, N_Vector tmp)
{
	realtype *JV = NV_DATA_S(Jv);
	realtype *Y = NV_DATA_S(u);
	realtype *V = NV_DATA_S(v);
  double EXP= exp ( (Y[2] - 1.0 ) / (EPS * Y[2]));
  double DEXP=1/(EPS * Y[2]*Y[2]);
	JV[0] = -.5*EXP* (V[0]*Y[1]+V[1]*Y[0]+V[2]*Y[0]*Y[1]*DEXP);
	JV[1] =-JV[0]-V[1];
	JV[2] = V[1];
  return 0;
}



void postProcess(N_Vector solution)
{
	realtype * data = N_VGetArrayPointer(solution);
	int len 		= N_VGetLength(solution);
	for (int i = 0; i < len; i ++ )
	{
		data[i] = std::max(0.0, data[i]);
	}
	//std::cout << "I AM HERE\n";
	
	
	/* 
		Maximum Brilliance: 
		clipping is equivalent to: x = (|x| + x ) / 2
	
		tmp = N_VClone(solution)
		N_VScale(
			0.5,
			N_VLinearSum(N_VAbs(temp, temp), solution)
		)
	*/

}