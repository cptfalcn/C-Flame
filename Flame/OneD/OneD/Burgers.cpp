/*
 * ===========================================================================================
 * 
 * This is implementaion of Burgers equation for Epic,
 * 
 *                                    u_t + u u_x = Mu u_xx
 * 
 * with Dirichlet boundary conditions and initial condition is given by
 * 
 *                               u = (sin(3 pi x))^3 (1-x)^(3/2)
 * 
 * All constants can be found in Burgers1dParallel.h. 
 * This impelentation is written so that it can work with MPI. 
 * For more details about variables and structures being used look at Burgers1dParallel.h.
 * 
 * -------------------------------------------------------------------------------------------
/*
 * ===========================================================================================
 * This is a Parallel implementation.
 * ===========================================================================================
 */

#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "Epic.h"
#include "EpiP2_KIOPS.h"
#include "Epi3SC_KIOPS.h"
#include "Epi3Ros.h"
#include "CreateIntegrators.h"
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
#include <nvector/nvector_parallel.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <cstdlib>
// Number of points in the grid
#define NEQ	1024  //was 1500
#define Nu	1e-3 //was 0.03

using namespace std;

struct UserData
{
    MPI_Comm comm;
    int myProcessor;
    int numProcessors;
    N_Vector scratch1;
    N_Vector scratch2;
	N_Vector V;
	int numPts;
	int numGrids;
	int eqPerProc;
	int localStart;
	realtype h;
};

struct EPICState
{
    N_Vector u;
    UserData *userData;
};

void Initialize(MPI_Comm, EPICState *, int);
void ReadData(N_Vector, string);
void Cleanup(EPICState *);
int RHS(realtype, N_Vector, N_Vector, void *);
int Jtv(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void *, N_Vector);
void Print2File(realtype *, EPICState *, string FileName);
void PrintErrStepTime2File(N_Vector, N_Vector, N_Vector, realtype, realtype, string, string);
void PrintErrStepTime2Console(N_Vector, N_Vector, N_Vector, realtype, realtype, string);
void PrintStats(int, IntegratorStats*);
void Print2Console(realtype*, EPICState *, string);
void PrintMPIReportIn(int, int);


void SortSerialVector(N_Vector, N_Vector);
void InitialConditions(N_Vector, int, EPICState *);
void Lu(N_Vector, N_Vector, UserData *, realtype, realtype);
void Mu(realtype, N_Vector, N_Vector, UserData *, realtype, realtype);
//End Prototypes


//Main
int main(int argc, char *argv[])
{
    // Opening statement
	std :: cout << BAR << "Starting Burgers.x" << BAR << "\n";
    // Initialize MPI.
    	MPI_Init(&argc, &argv);

    	int numProcessors, myProcessor;
    	MPI_Comm comm = MPI_COMM_WORLD;
    	MPI_Comm_size(comm, &numProcessors);
    	MPI_Comm_rank(comm, &myProcessor);

    // Declaration of start and finish time
    	static const realtype InitTime = 0.0;
    	static const realtype FinalTime = 0.1;//.01
    	static const int NumBands = 3;
	int size= 1024;

    	EPICState *EPI2State	= new EPICState;
	EPICState *RefState	= new EPICState;
	EPICState *EPI3VState	= new EPICState;
	EPICState *EPI3State	= new EPICState;
	EPICState *EPIP2State   = new EPICState;
	EPICState *EPI3RosState	= new EPICState;
   //Create integrator stats
    	IntegratorStats * integratorStats 	= NULL;
	IntegratorStats * integratorStats2 	= NULL;
	IntegratorStats * integratorStats5 	= NULL;
   //Initialize States
    	Initialize(comm, EPI2State,	size);//grids and num interior pts, makes delx easy!!
	Initialize(comm, RefState,	size);
	Initialize(comm, EPI3VState,	size);
	Initialize(comm, EPI3State,	size);
	Initialize(comm, EPIP2State,	size);
	Initialize(comm, EPI3RosState,	size);
	N_VScale(1.0, RefState->u, RefState->userData->scratch1);
    //Set data pointers
	realtype * EPI2DATA	= NV_DATA_P(EPI2State->u);
	realtype * REFDATA	= NV_DATA_P(RefState->u);
	realtype * EPI3VDATA	= NV_DATA_P(EPI3VState->u);
	realtype * EPI3DATA	= NV_DATA_P(EPI3State->u);
	//Print2File(UDATA, problemState, "AdvDiffIC.txt");


    // Create the integrator. In this case we use EpiRK5C, for other integrators just use their name.
    	const int MaxKryIts = 500;

    	Epi2_KIOPS * EPI2 	= new Epi2_KIOPS(RHS, Jtv, EPI2State->userData, MaxKryIts, EPI2State->u, NEQ);

	Epi3_KIOPS * REFINT 	= new Epi3_KIOPS(RHS, Jtv, RefState->userData, MaxKryIts, RefState->u, NEQ);

	Epi3SC_KIOPS * EPI3V 	= new Epi3SC_KIOPS(RHS, Jtv, EPI3VState->userData, MaxKryIts, EPI3VState->u, NEQ);

	EpiP2_KIOPS * EPIP2 	= new EpiP2_KIOPS(RHS, Jtv, EPIP2State->userData, MaxKryIts, EPIP2State->u, NEQ);

	//EpiRK4SV * EPIRK4V = new EpiRK4SV(RHS, Jtv, problemState->userData, MaxKryIts, problemState3->u, NEQ);

	Epi3_KIOPS * EPI3	= new Epi3_KIOPS(RHS, Jtv, EPI3State->userData, MaxKryIts, EPI3State->u, NEQ);

	Epi3Ros_KIOPS * Ros3 	= new Epi3Ros_KIOPS(RHS, Jtv, EPI3RosState->userData, MaxKryIts, EPI3RosState->u, NEQ);

    // Set run parameters.
    	const realtype StepSize = RCONST(1e-5);		//Originall 1e-5
	realtype TestStep = 1e-5;
	realtype RefStep =  1e-7;
	realtype ModStepSize = TestStep*1e1;		//StepSize*10;
    	const realtype KrylovTol = RCONST(1.0e-8);	//wtf was it set 1e-17?
    	int startingBasisSizes[] = {10, 10};

	ReadData(RefState->u, "BurgersRef.txt");
    // Run the integrators

	//Base EPI2
	auto Start = std::chrono::high_resolution_clock::now();
    	integratorStats = EPI2->Integrate(ModStepSize,	InitTime,	FinalTime,	NumBands,
						EPI2State->u,	KrylovTol,	startingBasisSizes);

	auto Stop=std::chrono::high_resolution_clock::now();
        auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
	std :: cout << "Completed Epi2 Integration\n";

	//EPIP2
	auto StartEPIP2 = std::chrono::high_resolution_clock::now();
			EPIP2->Integrate(ModStepSize, InitTime, FinalTime, NumBands,
				EPIP2State->u, KrylovTol, startingBasisSizes);
	std :: cout << "Completed EpiP2 integration\n";

	//REFIntegrator
	/*
	auto Start2= std::chrono::high_resolution_clock::now();
			  REFINT->Integrate(RefStep,       InitTime,       FinalTime,      NumBands,
                                                RefState->u,        KrylovTol,      startingBasisSizes);

	auto Stop2=std::chrono::high_resolution_clock::now();
        auto Pass2 = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop2-Start2);
	std :: cout << "Completed Reference integrations.\n";
	*/
	//ReadData(RefState->u, "BurgersRef.txt");
	//EPI3V

	auto Start3 = std::chrono::high_resolution_clock::now();
	integratorStats2 =EPI3V->Integrate(ModStepSize,  FinalTime/10, KrylovTol, KrylovTol, InitTime,
						FinalTime, NumBands, startingBasisSizes,  EPI3VState->u);

	auto Stop3 =std::chrono::high_resolution_clock::now();
        auto Pass3 = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop3-Start3);
	std :: cout << "Completed Epi3SV integration\n";

	//EPI3
	auto Start4 = std::chrono::high_resolution_clock::now();
	integratorStats5 =EPI3->Integrate(ModStepSize, InitTime, FinalTime, NumBands, EPI3State->u,
						KrylovTol, startingBasisSizes);

	auto Stop4 =std::chrono::high_resolution_clock::now();
	auto Pass4 = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop4-Start4);
	std :: cout << "Completed Epi3 integration.\n";

	/*
			  Ros3->Integrate(ModStepSize, InitTime, FinalTime, NumBands, EPI3RosState->u,
						KrylovTol, startingBasisSizes);
	std :: cout << "Completed Rosenbrock3 integration.\n";
	*/

	/**/

	cout << BAR << "Integrations complete" << BAR << endl;
    // Print the statistics
	/*
	std :: cout << "EPI2 Run\n";
	PrintStats(myProcessor, integratorStats);
	std :: cout << "EPI3Var Run\n";
	PrintStats(myProcessor, integratorStats2);
	std :: cout << "EPI3 Run\n";
	PrintStats(myProcessor, integratorStats5);
	*/
    // Clean up the integrator
   	delete EPI2;
	delete REFINT;
	delete EPI3V;
	//delete EPIRK4V;
	delete EPI3;
	delete EPIP2;
    //Print to Console

	std :: cout << "\n" << BAR << "Print to console" << BAR  << "\n";
	//EPI2 vs reference
	PrintErrStepTime2File(EPI2State->u, RefState->u, EPI2State->userData->scratch1,
				Pass.count()/1e9, ModStepSize, "BurgersDebug.txt", "EPI2");
	PrintErrStepTime2Console(EPI2State->u, RefState->u, EPI2State->userData->scratch1,
                                Pass.count()/1e9, ModStepSize, "EPI2");
	//EPI3V vs reference experiment
	/**/
	PrintErrStepTime2File(EPI3VState->u, RefState->u, EPI3VState->userData->scratch1,
                                Pass3.count()/1e9, ModStepSize, "BurgersDebug.txt", "EPI3SV");
	PrintErrStepTime2Console(EPI3VState->u, RefState->u, EPI3VState->userData->scratch1,
                                Pass3.count()/1e9, ModStepSize, "EPI3SV");
	std :: cout << "TimeSteps: " << integratorStats2->numTimeSteps << endl << endl;

	/**/
	//EPI3 Vs reference experiment
	PrintErrStepTime2File(EPI3State->u, RefState->u, EPI3State->userData->scratch1,
                                Pass.count()/1e9, ModStepSize, "BurgersDebug.txt", "EPI3");
	PrintErrStepTime2Console(EPI3State->u, RefState->u, EPI3State->userData->scratch1,
                                Pass4.count()/1e9, ModStepSize, "EPI3");
	std :: cout << endl;


	//EPIP2 Vs Ref
	PrintErrStepTime2Console(EPIP2State->u, RefState->u, EPIP2State->userData->scratch1,
                                Pass.count()/1e9, ModStepSize, "EPIP2");

	//Ros3 vs Ref
	/*
	N_VLinearSum(1.0, EPI3RosState->u, -1.0, RefState->u, EPIP2State->userData->scratch1);
	std :: cout << "EPI3Ros Stepsize: " << ModStepSize << endl;
	std :: cout << "EPI3 L2 error norm: " <<
	sqrt(N_VDotProd(EPIP2State->userData->scratch1, EPIP2State->userData->scratch1))<< std::endl;

	std :: cout << "EPI3 Max error: " << N_VMaxNorm(EPIP2State->userData->scratch1) << endl;
	*/

	//Final print
	FILE * FP = fopen("BurgersOut.txt", "w");
	N_VPrintFile_Parallel(EPI3RosState->u, FP);

    // Cleanup
    	Cleanup(EPI2State);
	Cleanup(RefState);
	Cleanup(EPI3VState);
	Cleanup(EPI3State);
	Cleanup(EPIP2State);
    	MPI_Finalize();
}

/*
 * ===========================================================================================
 * Routine Initialize
 *
 * Allocates the data.
 *
 * Inputs:
 * comm
 * state   pointer to state
 *
 * ===========================================================================================
 */

void Initialize(MPI_Comm comm, EPICState *state, int Eqs)
{
        int numProcessors, myProcessor;
        MPI_Comm_size(comm, &numProcessors);
        MPI_Comm_rank(comm, &myProcessor);
	int neqs = Eqs;// + 2); //+2 for the hanging ghosts.

        // Determine local vector length.
        int eqnsPerProcessor = neqs / numProcessors;  // note integer arithmetic
        int numRemainder = neqs - numProcessors * eqnsPerProcessor;
        int localNumEqns = (myProcessor < numRemainder) ? eqnsPerProcessor + 1 : eqnsPerProcessor;
        int myStartLocation = (myProcessor < numRemainder) ? myProcessor * localNumEqns : (myProcessor * eqnsPerProcessor) + numRemainder;

        // Allocate user data.
        state->userData = new UserData;
        state->userData->comm = comm;
        state->userData->myProcessor = myProcessor;
        state->userData->numProcessors = numProcessors;
        state->userData->scratch1 = N_VNew_Parallel(comm, localNumEqns, neqs);
        state->userData->scratch2 = N_VNew_Parallel(comm, localNumEqns, neqs);
	//state->userData->numGrids= Grids;
	state->userData->numPts=  Eqs;
	state->userData->localStart = myStartLocation;
	state->userData->h = 1.0/( neqs + 1 );

        // Allocate and produce initial conditions.
        state->u = N_VNew_Parallel(comm, localNumEqns, neqs);
        //CalculateInitialConditions(state->u, myStartLocation);
	InitialConditions(state->u, myStartLocation, state);
}

void Cleanup(EPICState *state)
{
        N_VDestroy_Parallel(state->u);
        N_VDestroy_Parallel(state->userData->scratch1);
        N_VDestroy_Parallel(state->userData->scratch2);
        delete state->userData;
}

/*
 * ===========================================================================================
 *
 * Routine  InitialConditions
 *
 * Computes initial condition. This rutine is called inside rutine Initialize.
 *
 * Input:
 * state   pointer to state
 *
 * ===========================================================================================
 *///Needs to be adjusted to have only the interior.
void InitialConditions(N_Vector u, int myStartLocation, EPICState * bin)
{
	realtype *data = NV_DATA_P(u);
        int myLength = NV_LOCLENGTH_P(u);
	realtype x = 0;
	int myIndex;
	//const realtype h = 1.0 / (bin->userData->numPts + 1);
	realtype h = bin->userData->h;

        for (int i = 0; i < myLength; i++)
        {
		myIndex = (myStartLocation + i ) % (bin->userData->numPts);
		//cout << myIndex << "\t\t";
		x = (myIndex ) * h ;//+ h/2;  //The closest to the reference
		//cout << x << endl;
                //realtype temp1 = sin(3 * M_PI * x);
                //realtype temp2 = pow(1-x, 3.0/2.0);
                //data[i] = temp1*temp1*temp1*temp2;
		realtype temp1 = pow( x-0.3, 2);
		realtype temp2 = 2* pow(.05, 2);
		data[i] = exp( -1* temp1/ temp2);
        }
	//Print2Console(data, bin, "Pre-Integration f(x)");

}

void Print2File(realtype* UDATA, EPICState * box, string FileName)
{
	int neq = box->userData->numGrids * (box->userData->numPts + 2); //+2 for the hanging ghosts.
        ofstream myfile( FileName, std::ios_base::app);
        for( int i = 0 ; i < neq ; i ++ )
                myfile << UDATA[i] << "\t\t";
	myfile << endl;
        myfile.close();
}

void PrintStats(int myProcessor, IntegratorStats * integratorStats)
{
    if (myProcessor == 0 && integratorStats != NULL )
    {
        printf("Run stats:\n");
        integratorStats->PrintStats();
        printf("\n\n");
    }
}

void Print2Console(realtype * DATA, EPICState * box, string Output)
{
	cout << "\t\t" << Output << endl;
	int neq = box->userData->numGrids * (box->userData->numPts + 2); //+2 for the hanging ghosts.
        for( int i = 0 ; i < neq ; i ++ )
		cout << DATA[i] << "\t\t";
	cout << endl;
}

void PrintMPIReportIn(int maxPro, int thisPro)
{
	cout << "============================Report======================================\n";
	cout << "Process number " << (thisPro + 1) << "/" << (maxPro) << " reported in.\n";
}

/*
 * ===========================================================================================
 * 
 * Routine SortSerialVector
 * 
 * Utility routine used when serializing a parallel vector out to a serial one.
 * Each processor's data is gathered to some root node, but is for a 2D patch of the global
 * grid.
 * The data must be sorted to be contiguous over the entire 2D region, and this routine does
 * that.
 * 
 * Inputs:
 * source   
 * target
 * 
 * ===========================================================================================
 */

void SortSerialVector(N_Vector source, N_Vector target)
{
        realtype *sourceData = NV_DATA_S(source);
        realtype *targetData = NV_DATA_S(target);

        for (int i = 0; i < NEQ; i++)
        {
                targetData[i] = sourceData[i];
        }
}

/*
 * ===========================================================================================
 * 
 * Routine Lu
 * 
 * This routine computes Laplacian.
 * It is used in functions RHS and Jtv.
 * 
 * Inputs:
 * u          input vector
 * result     Laplacian of u
 * userData   processor informations
 * 
 * ===========================================================================================
 */
/**/
void Lu(N_Vector u, N_Vector result, UserData *userData, realtype leftNeighbor, realtype rightNeighbor)
{
        int myProcessor = userData->myProcessor;
        int numProcessors = userData->numProcessors;
        int maxProcessor = numProcessors - 1;

        realtype *uData = NV_DATA_P(u);
        realtype *resultData = NV_DATA_P(result);
        int localLength = NV_LOCLENGTH_P(u);  // assume local length same for u and result

        //realtype h = 1.0 / (static_cast<realtype>(NEQ) + 1);
	realtype h = userData->h;
        realtype c = Nu / (h * h);
	//c = c/Nu;

        // Left edge node.
        if (myProcessor == 0)
        {
                resultData[0] = c * (-2*uData[0] + uData[1]);
        }
        else
        {
                resultData[0] = c * (leftNeighbor - 2*uData[0] + uData[1]);
        }

        // Interior nodes.
        for (int i = 1; i < localLength - 1; i++)
        {
                resultData[i] = c * (uData[i-1] - 2*uData[i] + uData[i+1]);
        }

        // Right edge node.
        if (myProcessor == maxProcessor)
        {
                resultData[localLength - 1] = c * (uData[localLength - 2] - 2*uData[localLength - 1]);
        }
        else
        {
                resultData[localLength - 1] = c * (uData[localLength - 2] - 2*uData[localLength - 1] + rightNeighbor);
        }
}



/*
 * ===========================================================================================
 * 
 * Routine Mu
 * 
 * This routine discretizeds u u_x term in Burgers equation.
 * Constant c passed in because different uses need different c values.
 *
 * Inputs:
 * c               constant defined as -1/(4 h) where h is spacing between nodes
 * u               input vector
 * result          result
 * userData
 * leftNeighbor    input for processor comunication
 * rightNeighbor   input for processor comunication
 * 
 * ===========================================================================================
 */
/**/
void Mu(realtype c, N_Vector u, N_Vector result, UserData *userData, realtype leftNeighbor, realtype rightNeighbor)
{
        int myProcessor = userData->myProcessor;
        int numProcessors = userData->numProcessors;
        int maxProcessor = numProcessors - 1;

        realtype *uData = NV_DATA_P(u);
        realtype *resultData = NV_DATA_P(result);
        int localLength = NV_LOCLENGTH_P(u);  // assume local length same for u and result

        // Left edge node.
        if (myProcessor == 0)
        {
                resultData[0] = c * uData[1];
        }
        else
        {
                resultData[0] = c * (uData[1] - leftNeighbor);
        }

        // Interior nodes.
        for (int i = 1; i < localLength - 1; i++)
        {
                resultData[i] = c * (uData[i + 1] - uData[i - 1]);
        }

        // Right edge node.
        if (myProcessor == maxProcessor)
        {
                resultData[localLength - 1] = -c * uData[localLength - 2];
        }
        else
        {
                resultData[localLength - 1] = c * (rightNeighbor - uData[localLength - 2]);
        }
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
/**/
int RHS(realtype t, N_Vector u, N_Vector udot, void *userData)
{
        UserData *data = static_cast<UserData *>(userData);

        N_Vector scratch1 = data->scratch1;
        N_Vector scratch2 = data->scratch2;

        realtype *uData = NV_DATA_P(u);
        int localLength = NV_LOCLENGTH_P(u);  // assume local length same for u and result
        MPI_Comm comm = NV_COMM_P(u);

        int myProcessor = data->myProcessor;
        int numProcessors = data->numProcessors;
        int maxProcessor = numProcessors - 1;

        const realtype h = 1.0 / (static_cast<realtype>(NEQ) + 1);
        const realtype c = -1.0 / (4*h);

        // Send neighbor values.
        if (myProcessor != 0)
        {
                MPI_Send(&uData[0], 1, PVEC_REAL_MPI_TYPE, myProcessor - 1, 0, comm);
        }
        if (myProcessor != maxProcessor)
        {
                MPI_Send(&uData[localLength - 1], 1, PVEC_REAL_MPI_TYPE, myProcessor + 1, 0, comm);
        }

        // Receive neighbor values.
        realtype leftNeighbor = 0.0;
        realtype rightNeighbor = 0.0;
        MPI_Status mpiStatus;
        if (myProcessor != 0)
        {
                MPI_Recv(&leftNeighbor, 1, PVEC_REAL_MPI_TYPE, myProcessor - 1, 0, comm, &mpiStatus);
        }
        if (myProcessor != maxProcessor)
        {
                MPI_Recv(&rightNeighbor, 1, PVEC_REAL_MPI_TYPE, myProcessor + 1, 0, comm, &mpiStatus);
        }

        N_VProd(u, u, scratch1);
        realtype leftNeighborSquared = leftNeighbor * leftNeighbor;
        realtype rightNeighborSquared = rightNeighbor * rightNeighbor;
        Mu(c, scratch1, scratch2, data, leftNeighborSquared, rightNeighborSquared);

        Lu(u, udot, data, leftNeighbor, rightNeighbor);
        N_VLinearSum(1.0, udot, 1.0, scratch2, udot);

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
/**/
int Jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void *userData, N_Vector tmp)
{
        UserData *data = static_cast<UserData *>(userData);

        N_Vector scratch1 = data->scratch1;
        N_Vector scratch2 = data->scratch2;

        realtype *vData = NV_DATA_P(v);
        realtype *scratch1Data = NV_DATA_P(scratch1);
        int localLength = NV_LOCLENGTH_P(u);  // assume local length same for u and result
        MPI_Comm comm = NV_COMM_P(u);

        int myProcessor = data->myProcessor;
        int numProcessors = data->numProcessors;
        int maxProcessor = numProcessors - 1;

        const realtype h = 1.0 / (NEQ + 1);
        const realtype c = -1.0 / (2 * h);

        N_VProd(u, v, scratch1);

        // Communicate the neighbor values for both v and scratch1.
        // First send my values.
        realtype myLeft[2], myRight[2];
        myLeft[0] = scratch1Data[0];
        myLeft[1] = vData[0];
        myRight[0] = scratch1Data[localLength - 1];
        myRight[1] = vData[localLength - 1];
        if (myProcessor != 0)
        {
                MPI_Send(myLeft, 2, PVEC_REAL_MPI_TYPE, myProcessor - 1, 0, comm);
        }
        if (myProcessor != maxProcessor)
        {
                MPI_Send(myRight, 2, PVEC_REAL_MPI_TYPE, myProcessor + 1, 0, comm);
        }

        // Now receive the neighbor values for v and scratch1.
        realtype leftNeighbors[2], rightNeighbors[2];
        leftNeighbors[0] = 0.0; leftNeighbors[1] = 0.0;
        rightNeighbors[0] = 0.0; rightNeighbors[1] = 0.0;
        MPI_Status mpiStatus;
        if (myProcessor != 0)
        {
                MPI_Recv(leftNeighbors, 2, PVEC_REAL_MPI_TYPE, myProcessor - 1, 0, comm, &mpiStatus);
        }

        if (myProcessor != maxProcessor)
        {
                MPI_Recv(rightNeighbors, 2, PVEC_REAL_MPI_TYPE, myProcessor + 1, 0, comm, &mpiStatus);
        }

        Mu(c, scratch1, scratch2, data, leftNeighbors[0], rightNeighbors[0]);
        Lu(v, Jv, data, leftNeighbors[1], rightNeighbors[1]);
        N_VLinearSum(1.0, Jv, 1.0, scratch2, Jv);

        return 0;
}


void ReadData(N_Vector Ref, string file)
{
	string line;
	cout << "Getting Reference from: " << file << endl;
	realtype * REF = N_VGetArrayPointer(Ref);
	//realtype Test = 0;
	ifstream inputFile;
	inputFile.open(file);
	cout << "Is the file valid? " << (inputFile.good()? "Yes" : "No" ) << endl;
	int i = 0;
	while(!inputFile.eof())
	{
		getline(inputFile,line);
		//std :: cout << line << endl;
		//Test = stof(line);
		REF[i+1] = stod(line);		//This is the problem line, shifts the solution by one.
		i++;
		inputFile.peek();		//Set the flag to see if the next entry is good.
	}
	//cout << "There are " << i << " lines in " << file << endl;
	inputFile.close();


}

void PrintErrStepTime2File(N_Vector data, N_Vector ref, N_Vector scratch, realtype time, realtype step, string File,
				string method)
{
	ofstream myfile( File, std::ios_base::app);
	N_VLinearSum(1.0, data, -1.0, ref, scratch);
       	myfile << endl << method << " integration time: " << time  << std :: endl;
        myfile << method << " StepSize: " << step << endl;
        myfile << "L2 error norm: " << sqrt( N_VDotProd(scratch, scratch)) << std :: endl;
        myfile << std :: endl;
}

void PrintErrStepTime2Console(N_Vector data, N_Vector ref, N_Vector scratch, realtype time, realtype step, string method)
{
        N_VLinearSum(1.0, data, -1.0, ref, scratch);
        std :: cout << endl << method << " integration time: " << time  << std :: endl;
        std :: cout << method << " StepSize: " << step << endl;
        std :: cout << "L2 error norm: " << sqrt( N_VDotProd(scratch, scratch)) << std :: endl;

}
