	//====================//
	//Test the Derivatives//
	//====================//

#include "Derivatives.h"
#include "Flow.h"
#include <iomanip>
#include <fstream>
#define BAR "===================="

void SetGrid	(realtype *, realtype, int, realtype);
void SetFunction(realtype *, realtype *, realtype, int, int);

double CheckDerErr (realtype *, realtype *, int, int, int);
double CheckDerDerErr(realtype *, realtype *, int, int, int);

void PrintState	(realtype *, int, string);

using namespace std;

int main()
{
	int OPT=		4;
	realtype OFFSET =	0;
	realtype DELTA= 	1e-4;
	int SIZE = 		round(1/DELTA) + 3;
	realtype Err = 		0;

	ofstream myfile("DerivativeConvergence.txt", std::ios_base::app);

	cout << SIZE << " slots\n";
	//Declare F(X) and F(V) with their pointers
	N_Vector FX = 		N_VNew_Serial(SIZE);
	realtype * FX_PTR =	NV_DATA_S(FX);
	N_Vector FXSTAG = 	N_VNew_Serial(SIZE-1);
	realtype * FXS_PTR =	NV_DATA_S(FXSTAG);

	//Set Data grid X and staggered supergrid V with pointers
	N_Vector X =		N_VNew_Serial(SIZE);
	N_Vector XSTAG =	N_VNew_Serial(SIZE-1); //We lose a slot when we stagger down
	realtype *XPTR =	NV_DATA_S(X);
	realtype *XSTAG_PTR =	NV_DATA_S(XSTAG);
	//Set the result vector and pointer
	N_Vector RES =		N_VNew_Serial(SIZE-2);//No Ghosts
	realtype *RES_PTR =	NV_DATA_S(RES);
	//Set the velocity sized result vector with pointer
	N_Vector RESSTAG =	N_VNew_Serial(SIZE-1);
	realtype *RESSTAG_PTR = NV_DATA_S(RESSTAG);

	cout<<"Generating test data...\n";
	//Set grids to be X=0:Delta:1, V= 0-Delta/2 : Delta : 1+Delta/2
	SetGrid(XPTR, DELTA, SIZE, OFFSET);
	SetGrid(XSTAG_PTR, DELTA, SIZE-1, -DELTA/2);//Shift to the right Delta/2
	SetFunction(FX_PTR, XPTR, DELTA, SIZE, OPT);
	//Print the Various states
	PrintState(XPTR, SIZE, "X");
	PrintState(FX_PTR, SIZE, "f(X)");
	PrintState(XSTAG_PTR, SIZE-1 , "X_Stag");
	SetFunction(FXS_PTR, XSTAG_PTR, DELTA, SIZE-1, OPT);
	PrintState(FXS_PTR, SIZE-1, "f(X_Stag)");


	//Test some derivatives
	//Double Derivative
	DDCentered(RES_PTR, FX_PTR, DELTA, SIZE);
	//PrintState(RES_PTR, SIZE-2, "Velocity, Second Derivative-Centered");
	Err= CheckDerDerErr(RES_PTR, XPTR, SIZE-2, OPT, 0);
	cout << setprecision(17) << "\nError of the second derivative: \t" << Err << "\n\n";
	myfile << Err << "\t\t";


	//Derivative Centered
	DCentered(RES_PTR, FX_PTR, DELTA, SIZE);
	//PrintState(RES_PTR, SIZE-2, "Velocity First Derivative-Centered");
	Err = CheckDerErr(RES_PTR, XPTR, SIZE-2, OPT, 0);
	cout << "\nError of DCentered:\t" << Err << "\n\n";
	myfile << Err << "\t\t";

	//Derivative, State to Velocity, Stagger up in size.
	//Take FXStaggered and set it up onto the non staggered grid
	DStaggerUp(RES_PTR, FXS_PTR, DELTA, SIZE-2);
	//PrintState(RES_PTR, SIZE-2,  "State to Velocity Derivative-Centered");
	Err = CheckDerErr(RES_PTR, XPTR, SIZE-2, OPT, 0);
	cout << "\nError of Stagger up to Velocity: \t " << Err << "\n\n";
	myfile << Err << "\t\t";

	//Derivative, Velocity to State, stagger down
	DStaggerDown(RESSTAG_PTR, FX_PTR, DELTA, SIZE-1);
	//PrintState(RESSTAG_PTR, SIZE-3, "Velocity to State Derivative-Centered");
	Err = CheckDerErr(RESSTAG_PTR, XSTAG_PTR, SIZE-3, OPT, 0);
	cout << "\nError of Stagger Down to State: \t " << Err << "\n\n";
	myfile << Err << "\t\t" << DELTA << endl;

	//Attempt to stagger the function F
	//StaggerV(FV_PTR, FVS_PTR, SIZE);
	//PrintState(FVS_PTR, SIZE-2,"Staggered averaged velocity function");

	//Clean up and free memory
	myfile.close();
	N_VDestroy_Serial(RES);
	N_VDestroy_Serial(XSTAG);
	N_VDestroy_Serial(X);
	N_VDestroy_Serial(FX);
	N_VDestroy_Serial(FXSTAG);
	N_VDestroy_Serial(RESSTAG);
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


double CheckDerErr(realtype *DATA, realtype *X , int SIZE, int OPT, int STAG)
{//Checks First Derivative Error
	N_Vector ErrVec =	N_VNew_Serial(SIZE); //Set the size
	realtype * ErrData =	NV_DATA_S(ErrVec);   //
	realtype ErrNorm = 0;
	int j = 0;

	for(int i = 0 ; i< SIZE; i++)
	{
		j=i+1-STAG;
		switch(OPT)
		{
		case 1:
			ErrData[i]=DATA[i] - 0;
			//OUT[i] = 1;
			break;
		case 2:
			ErrData[i] = DATA[i] - 1;
			//OUT[i] = X[i];
			break;
		case 3:
			ErrData[i] = DATA[i] - 2 * X[j];
			//OUT[i] = pow(X[i],2);
			break;
		case 4:
			ErrData[i] = DATA[i] - M_PI * cos (M_PI*X[j]);
			//OUT[i] = sin(M_PI*X[i]);
			break;
		case 5:
			ErrData[i] = DATA[i] - 3 * pow(X[j] , 2) ;
			//OUT[i] = pow(X[i], 3);
			break;
		}
	}
	//ErrNorm = sqrt( N_VDotProd(ErrVec, ErrVec) );
	ErrNorm = N_VMaxNorm(ErrVec);
	N_VDestroy_Serial (ErrVec);
	return ErrNorm;
}

double CheckDerDerErr(realtype *DATA, realtype *X , int SIZE, int OPT, int STAG)
{//Checks First Derivative Error
        N_Vector ErrVec =       N_VNew_Serial(SIZE); //Set the size
        realtype * ErrData =    NV_DATA_S(ErrVec);   //
        realtype ErrNorm = 0;
        int j = 0;

        for(int i = 0 ; i< SIZE; i++)
        {
                j=i+1-STAG;
                switch(OPT)
                {
                case 1:
                        ErrData[i]=DATA[i] - 0;
                        //OUT[i] = 1;
                        break;
                case 2:
                        ErrData[i] = DATA[i] - 0;
                        //OUT[i] = X[i];
                        break;
                case 3:
                        ErrData[i] = DATA[i] - 2 ;
                        //OUT[i] = pow(X[i],2);
                        break;
                case 4:
                        ErrData[i] = DATA[i] + M_PI * M_PI * sin (M_PI*X[j]);
                        //OUT[i] = sin(M_PI*X[i]);
                        break;
                case 5:
                        ErrData[i] = DATA[i] - 6 * pow(X[j] , 1) ;
                        //OUT[i] = pow(X[i], 3);
                        break;
                }
        }
        ErrNorm = sqrt( N_VDotProd(ErrVec, ErrVec) );
        N_VDestroy_Serial (ErrVec);
        return ErrNorm;
}


void PrintState(realtype * PTR, int SIZE, string NAME)
{
	cout << BAR << NAME << BAR <<endl;
	for (int i = 0 ; i < SIZE ; i++)
		cout << PTR[i] << "\t";
	cout << endl;

}
