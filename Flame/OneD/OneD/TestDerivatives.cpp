	//====================//
	//Test the Derivatives//
	//====================//

#include "Derivatives.h"
#include "Flow.h"
#define BAR "===================="

void SetGrid	(realtype *, realtype, int, realtype);
void SetFunction(realtype *, realtype *, realtype, int, int);
void CheckErr 	(realtype *, realtype *, int, int);

void PrintState	(realtype *, int, string);

using namespace std;

int main()
{
	int OPT=		1;
	realtype OFFSET =	0;
	realtype DELTA= 	1e-1;
	int SIZE = 		round(1/DELTA) + 3;

	cout << SIZE << " slots\n";
	N_Vector FX = 		N_VNew_Serial(SIZE);
	realtype *FX_PTR =	NV_DATA_S(FX);
	N_Vector FV = 		N_VNew_Serial(SIZE+1);
	realtype * FV_PTR =	NV_DATA_S(FV);
	//Declare Staggered
	N_Vector FVS =		N_VNew_Serial(SIZE);
	realtype * FVS_PTR =	NV_DATA_S(FVS);

	N_Vector X =		N_VNew_Serial(SIZE);
	N_Vector V =		N_VNew_Serial(SIZE+1);
	realtype *XPTR =	NV_DATA_S(X);
	realtype *VPTR =	NV_DATA_S(V);

	N_Vector RES =		N_VNew_Serial(SIZE-2);//No Ghosts
	realtype *RES_PTR =	NV_DATA_S(RES);

	N_Vector RESV =		N_VNew_Serial(SIZE-1);
	realtype *RESV_PTR =	NV_DATA_S(RESV);

	cout<<"Generating test data...\n";

	SetGrid(XPTR, DELTA, SIZE, OFFSET);
	SetGrid(VPTR, DELTA, SIZE+1, DELTA/2);
	SetFunction(FX_PTR, XPTR, DELTA, SIZE, OPT);

	PrintState(XPTR, SIZE, "X");
	PrintState(FX_PTR, SIZE, "f(x)");


	PrintState(VPTR, SIZE+1 , "V");
	SetFunction(FV_PTR, VPTR, DELTA, SIZE+1, OPT);
	PrintState(FV_PTR, SIZE+1, "f(v)");


	//Test some derivatives
	DDCentered(RES_PTR, FX_PTR, DELTA, SIZE);
	PrintState(RES_PTR, SIZE-2, "Second Derivative-Centered");

	DCentered(RES_PTR, FX_PTR, DELTA, SIZE);
	PrintState(RES_PTR, SIZE-2, "First Derivative-Centered");

	DStaggered(RESV_PTR, FV_PTR, DELTA, SIZE+1);
	PrintState(RESV_PTR, SIZE-1,  "V to State Derivative-Centered");

	DStaggeredDown(RES_PTR, FX_PTR, DELTA, SIZE);
	PrintState(RES_PTR, SIZE-2, "State to Velocity Derivative-Centered");

	//Attempt to stagger
	StaggerV(FV_PTR, FVS_PTR, SIZE);
	PrintState(FVS_PTR, SIZE-2,"Staggered averaged velocity function");


	//Clean up and free memory
	N_VDestroy_Serial(RES);
	N_VDestroy_Serial(V);
	N_VDestroy_Serial(X);
	N_VDestroy_Serial(FX);
	N_VDestroy_Serial(FV);
	N_VDestroy_Serial(RESV);
	N_VDestroy_Serial(FVS);
	return 0;

}


//=========
//Functions
//=========
void SetGrid(realtype * GRID, realtype DELTA, int SIZE, realtype OFFSET)
{
	for( int i = 0; i < SIZE ; i++)
		GRID[i]= (i-1)*DELTA - OFFSET;


}

void SetFunction(realtype * OUT, realtype * X, realtype DELTA, int SIZE, int OPT)
{
	for(int i = 0 ; i< SIZE; i++)
	{
		switch(OPT)
		{
		case 1:
			OUT[i] = 1;
			break;
		case 2:
			OUT[i] = X[i];
			break;
		case 3:
			OUT[i] = pow(X[i],2);
			break;
		case 4:
			OUT[i] = sin(M_PI*X[i]);
			break;
		case 5:
			OUT[i] = pow(X[i], 3);
			break;
		}
	}
}

void CheckErr(realtype *DATA, realtype *X , int SIZE, int OPT)
{
	for(int i = 0 ; i< SIZE; i++)
	{
		switch(OPT)
		{
		case 1:
			//OUT[i] = 1;
			break;
		case 2:
			//OUT[i] = X[i];
			break;
		case 3:
			//OUT[i] = pow(X[i],2);
			break;
		case 4:
			//OUT[i] = sin(M_PI*X[i]);
			break;
		case 5:
			//OUT[i] = pow(X[i], 3);
			break;
		}
	}
}
void PrintState(realtype * PTR, int SIZE, string NAME)
{
	cout << BAR << NAME << BAR <<endl;
	for (int i = 0 ; i < SIZE ; i++)
		cout << PTR[i] << "\t";
	cout << endl;

}
