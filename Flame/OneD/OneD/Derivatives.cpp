//
//Doing the second derivative in centered format as a matrix vector product.
//
#include "Derivatives.h"
#define BAR "===================="
using namespace std;


//================================================================
//Velocity staggered Derivative onto Scalars
//================================================================
void DStaggered(realtype * OUT, realtype * DATA, realtype DELTA, int SIZE)
{
	//Data is one size greater than the other grid
	// \ G | x | x | x | x | x | G /
	// Here we calculate the from x's to  |'s
	//Memory leaks in here fixed
	for (int i = 1 ; i < SIZE -1 ; i++)
		OUT[i-1] = ( DATA[i] - DATA[i-1] ) / (DELTA);
}

//=================================================================
//Scalars staggered derivative onto Velocity
//=================================================================
void DStaggeredDown(realtype * OUT, realtype * DATA, realtype DELTA, int SIZE)
{
	//Data is one size greater than the other grid
	// \ G | x | x | x | x | x | G /
	// Here we calculate the from x's to  |'s
	for ( int i = 1 ; i < SIZE -1 ; i++)
		OUT[i-1] = (DATA[i+1] - DATA[i]) / DELTA;
}

//=================================================================
//Second derivative, second order, scalar derivative. Use GHOST
//points to control different boundary conditions
//=================================================================
void DDCentered(realtype * OUT, realtype * DATA, realtype DELTA, int SIZE)
{
	for (int i =1; i < SIZE -1 ; i++)
		OUT[i-1] =(DATA[i+1]-2*DATA[i]+DATA[i-1]) / (DELTA*DELTA);
}


//==================================================================
//Second order first derivative calculations.
//Change the GHOST points to reflect the correct boundary condition.
//==================================================================
void DCentered(realtype * OUT, realtype * DATA, realtype DELTA, int SIZE)
{
	//cout << BAR << "Inside DCentered function!!!" << BAR << endl;
	for( int i =1 ; i < SIZE - 1 ; i++)
		OUT[i-1]=      (DATA[i+1]-DATA[i-1]) / (2*DELTA);

}



//=============================================================
//  __ __    ___   _____      __   __ ___   ___
// |  V  | ./ _ \.|_   _| ___ \ \ / /| __) / __)
// | .V. | | (_) |  | |  |___| \ v / | __)( (__
// |_| |_| |_/ \_|  |_|         \_/  |___) \___)
//============================================================

void MatVec(int number_of_equations, realtype * JacD, N_Vector x, N_Vector tmp, realtype * JV)
{
        realtype * TMP  = NV_DATA_S(tmp);
        for (int i=0; i< number_of_equations; i++)//for each state
        {
                for(int j=0; j<number_of_equations; j++)//marches across the column
                {
                        TMP[j]=JacD[j+i*number_of_equations];//Stored row major
                }
                JV[i]=N_VDotProd(tmp,x);
        }
}

//New material

void PrintVec(int SIZE, realtype * DATA, string OPT)
{
	cout<<"Printing " << OPT << endl;
	for (int i =0; i < SIZE; i++)
		cout << DATA[i] << "\t";
	cout << endl;

}
