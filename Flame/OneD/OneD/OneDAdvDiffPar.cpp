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
#include <nvector/nvector_parallel.h>
#include <iostream>
#include <fstream>
#include <iomanip>
// Number of points in the grid
#define NEQ	100 //was 1500
#define Nu	0.03

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
};

struct EPICState
{
    N_Vector u;
    UserData *userData;
};

void Initialize(MPI_Comm, EPICState *, int, int);
void Cleanup(EPICState *);
int RHS(realtype, N_Vector, N_Vector, void *);
int RHS2(realtype, N_Vector, N_Vector, void *);
int Jtv(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void *, N_Vector);
int JtV2(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void *, N_Vector);
void Print2File(realtype *, EPICState *, string FileName);
void PrintStats(int, IntegratorStats*);
void Print2Console(realtype*, EPICState *, string);
void PrintMPIReportIn(int, int);


void SortSerialVector(N_Vector, N_Vector);
void InitialConditions(N_Vector, int, EPICState *);
//void Lu(N_Vector, N_Vector, UserData *, realtype, realtype);
void Diff(N_Vector, N_Vector, UserData *, realtype, realtype, int);
void Adv(N_Vector,N_Vector, N_Vector, UserData *, realtype, realtype, int);
//void Mu(realtype, N_Vector, N_Vector, UserData *, realtype, realtype);
//End Prototypes


//Main
int main(int argc, char *argv[])
{
    // Initialize MPI.
    	MPI_Init(&argc, &argv);

    	int numProcessors, myProcessor;
    	MPI_Comm comm = MPI_COMM_WORLD;
    	MPI_Comm_size(comm, &numProcessors);
    	MPI_Comm_rank(comm, &myProcessor);

    // Declaration of start and finish time
    	static const realtype InitTime = 0.0;
    	static const realtype FinalTime = .1;//.01
    	static const int NumBands = 3;

    	EPICState *problemState = new EPICState;
    	IntegratorStats * integratorStats = NULL;
    	Initialize(comm, problemState, 1, 100);//grids and num interior pts, makes delx easy!!
	realtype * UDATA= NV_DATA_P(problemState->u);
	//Print2File(UDATA, problemState, "AdvDiffIC.txt");

    // Create the integrator. In this case we use EpiRK5C, for other integrators just use their name.
    	const int MaxKryIts = 500;
    	Epi2_KIOPS * integrator = new Epi2_KIOPS(RHS2,	JtV2,
					problemState->userData, MaxKryIts,
					problemState->u,	NEQ);
    // Set run parameters.
    	const realtype StepSize = RCONST(1e-5);
    	const realtype KrylovTol = RCONST(1.0e-8);//wtf was it set 1e-17?
    	int startingBasisSizes[] = {10, 10};
    // Run the integrator
	/**/
    	integratorStats = integrator->Integrate(StepSize,	InitTime,	FinalTime,	NumBands,
						problemState->u,	KrylovTol,	startingBasisSizes);
	/**/
    // Print the statistics
	//PrintStats(myProcessor, integratorStats);
    // Clean up the integrator
   	delete integrator;
    //Print to file
	PrintMPIReportIn(numProcessors, myProcessor);
	Print2File(UDATA, problemState, "AdvDiffDataDebug.txt");
	Print2Console(UDATA, problemState, "\n\tPost Integration f(x)");
    // Cleanup
    	Cleanup(problemState);
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

void Initialize(MPI_Comm comm, EPICState *state, int Grids, int Eqs)
{
        int numProcessors, myProcessor;
        MPI_Comm_size(comm, &numProcessors);
        MPI_Comm_rank(comm, &myProcessor);
	int neqs = Grids * (Eqs + 2); //+2 for the hanging ghosts.
	int vEqs = Eqs + 1 ;//But these are staggered.

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
	state->userData->numGrids= Grids;
	state->userData->numPts=  Eqs;
	state->userData->localStart = myStartLocation;
	//state->userData->V = N_VNew_Serial( (Eqs+1) * Grids );

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
 */
void InitialConditions(N_Vector u, int myStartLocation, EPICState * bin)
{
	realtype *data = NV_DATA_P(u);
        int myLength = NV_LOCLENGTH_P(u);
	realtype x = 0;
	int myIndex;
	const realtype h = 1.0 / (bin->userData->numPts );

        for (int i = 0; i < myLength; i++)
        {
		myIndex = (myStartLocation + i ) % (bin->userData->numPts + 2 );
		cout << myIndex << "\t\t";
		x = ( myIndex -1 ) * h + h/2;
		//cout << x << endl;
                realtype temp1 = cos(2 * M_PI * x);
                //realtype temp2 = pow(1-x, 3.0/2.0);
                data[i] = temp1;
        }
	Print2Console(data, bin, "Pre-Integration f(x)");

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
/*
void Lu(N_Vector u, N_Vector result, UserData *userData, realtype leftNeighbor, realtype rightNeighbor)
{
        int myProcessor = userData->myProcessor;
        int numProcessors = userData->numProcessors;
        int maxProcessor = numProcessors - 1;

        realtype *uData = NV_DATA_P(u);
        realtype *resultData = NV_DATA_P(result);
        int localLength = NV_LOCLENGTH_P(u);  // assume local length same for u and result

        realtype h = 1.0 / (static_cast<realtype>(NEQ) + 1);
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
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//To be removed




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
/*
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

*/
//------------------------------------------
//
//Discretizes :  -  v u'
//Note V needs to be averaged to fit on the u grid
//
//-------------------------------------------
//Advection
//-------------------------------------------
//For this simple test have a constant Velocity of 1, so no MP
void Adv(N_Vector vAve, N_Vector u, N_Vector result, UserData *userData, realtype leftNeighbor,
			realtype rightNeighbor, int localStart)
{
        int myProcessor = userData->myProcessor;
        int numProcessors = userData->numProcessors;
        int maxProcessor = numProcessors - 1;

        realtype *uData = NV_DATA_P(u);
	//realtype *vData = NV_DATA_P(vAve);
        realtype *resultData = NV_DATA_P(result);
	int numPts = userData->numPts;
        int localLength = NV_LOCLENGTH_P(u);  // assume local length same for u and result
	realtype delx = 1.0 / userData->numPts;
	realtype divisor =	- 1.0 / ( 2 * delx );

        // Left edge node- boundary condition
        if (myProcessor == 0)
		resultData[0] = divisor * uData[1];
        else
		resultData[0] = divisor * (uData[1] - leftNeighbor);

        // Interior nodes.
	// cout << "local start: " << localStart << endl;
        for (int i = 1; i < localLength - 1; i++)
	{
		/**/
		if((localStart+i)% (numPts + 2 ) == 0) //left boundary
		{
			//cout << "index " <<  i << " acted up for left boundary check" << endl;
			resultData[i] = divisor * uData[i + 1];
		}
		else if((localStart+i)% (numPts+ 2) == numPts+1)//right boundary
		{
			//cout << "index " << i << " acted up for right boundary check" << endl;
			resultData[i] = divisor * (uData[i] - uData[i - 1]);
		}
		else
		/**/
                	resultData[i] = divisor* (uData[i + 1] - uData[i - 1]);
	}
        // Right edge node- boundary condition
        if (myProcessor == maxProcessor)
		resultData[localLength - 1] = divisor * (uData[localLength-1] - uData[localLength - 2]);
        else
                resultData[localLength - 1] = divisor * (rightNeighbor - uData[localLength - 2]);
}
//------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------
//Diffusion
//------------------------------------------------------------------------------------------
void Diff(N_Vector u, N_Vector result, UserData *userData, realtype leftNeighbor, realtype rightNeighbor,
		int localStart)
{
        int myProcessor = userData->myProcessor;
        int numProcessors = userData->numProcessors;
        int maxProcessor = numProcessors - 1;

        realtype *uData = NV_DATA_P(u);
        realtype *resultData = NV_DATA_P(result);
        int localLength = NV_LOCLENGTH_P(u);  // assume local length same for u and result
	int numPts = userData->numPts;

        realtype delx = 1.0 / (userData->numPts);
	//cout << endl << "Delx = " << delx << endl;
        realtype c = 1.0 / (delx * delx);
	//c = userData->numPts*userData->numPts;
	c = c * Nu;

        // Left edge node.
        if (myProcessor == 0)
                resultData[0] = c * (-2*uData[0] + uData[1]);
        else
                resultData[0] = c * (leftNeighbor - 2*uData[0] + uData[1]);
        // Interior nodes.
        for (int i = 1; i < localLength - 1; i++)
        {
		if((localStart+i)%(numPts + 2) == 0)//left
			resultData[i] = c * (-2*uData[i] + uData[i+1]);
		else if ( (localStart+i)% (numPts + 2) == (numPts + 1) )//right
			resultData[i] = c * ( uData[i-1] - uData[i] );
		else
                	resultData[i] = c * (uData[i-1] - 2*uData[i] + uData[i+1]);
        }

        // Right edge node.
        if (myProcessor == maxProcessor)
        {//Edit the node for 0 boundary condition
                resultData[localLength - 1] = c * (uData[localLength - 2] - uData[localLength - 1]);
        }
        else
        {
                resultData[localLength-1] = c*(uData[localLength-2] - 2*uData[localLength - 1] + rightNeighbor);
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
/*
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
*/
//============================================
//My Version of the Righthand side
//Remember the boundaries are cooked into this version
//============================================
int RHS2(realtype t, N_Vector u, N_Vector udot, void *userData)
{
        UserData *data = static_cast<UserData *>(userData);

	int  	ADV 	=	1;
	int 	DIFF	=	1;

        N_Vector scratch1 = data->scratch1;
        N_Vector scratch2 = data->scratch2;

        realtype *uData = NV_DATA_P(u);
	//realtype *vData = NV_DATA_P(data->V);
        int localLength = NV_LOCLENGTH_P(u);  // assume local length same for u and result
	int localStart 	= data->localStart;
        MPI_Comm comm 	= NV_COMM_P(u);
	//MPI_Comm commV	= NV_COMM_P(data->V);

        int myProcessor = data->myProcessor;
        int numProcessors = data->numProcessors;
        int maxProcessor = numProcessors - 1;

        // Send neighbor values.
        if (myProcessor != 0)
	{
                MPI_Send(&uData[0], 1, PVEC_REAL_MPI_TYPE, myProcessor - 1, 0, comm);
		//MPI_Send(&vData[0], 1, PVEC_REAL_MPI_TYPE, myProcessor - 1, 0, commV);
	}
        if (myProcessor != maxProcessor)
	{
                MPI_Send(&uData[localLength - 1], 1, PVEC_REAL_MPI_TYPE, myProcessor + 1, 0, comm);
		//MPI_Send(&vData[localLength - 1], 1, PVEC_REAL_MPI_TYPE, myProcessor + 1, 0, commV);
	}

        // Receive neighbor values.
        realtype leftNeighbor = 0.0;
        realtype rightNeighbor = 0.0;
        MPI_Status mpiStatus;
        if (myProcessor != 0)
                MPI_Recv(&leftNeighbor, 1, PVEC_REAL_MPI_TYPE, myProcessor - 1, 0, comm, &mpiStatus);
        if (myProcessor != maxProcessor)
                MPI_Recv(&rightNeighbor, 1, PVEC_REAL_MPI_TYPE, myProcessor + 1, 0, comm, &mpiStatus);

        //N_VProd(u, u, scratch1);
	//N_VProd(u, u, scratch1);
	realtype leftNeighborSquared = leftNeighbor * leftNeighbor;
	realtype rightNeighborSquared = rightNeighbor * rightNeighbor;
	//Adv
	if (ADV ==  1 )
        	Adv(data->V, u, scratch2, data, leftNeighbor, rightNeighbor, localStart);
	//Diff
	if (DIFF == 1 )
		Diff(u, udot, data, leftNeighbor, rightNeighbor, localStart);

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
/*
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
*/
int JtV2(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void *userData, N_Vector tmp)
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
	int localStart	=	data->localStart;

        const realtype h = 1.0 / (data-> numPts);
        const realtype c = -1.0 / (2 * h);

        //N_VProd(u, v, scratch1);
	N_VScale(1.0, v, scratch1);

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

        Adv(scratch2, scratch1, scratch2, data, leftNeighbors[1], rightNeighbors[1], localStart);
        Diff(v, Jv, data, leftNeighbors[1], rightNeighbors[1], localStart);
        N_VLinearSum(1.0, Jv, 1.0, scratch2, Jv);

        return 0;
}
