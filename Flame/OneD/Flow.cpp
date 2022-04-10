#include "Derivatives.h"
//Recall all Data has the ghost points included
void StaggerV(realtype * V, realtype * Stag_V, int SIZE)
{	//Recall that V is SIZE+1
	//  / G | x | x | x | x | x | G \
	//  Moving the  |'s to the x's
	for(int i = 0; i < SIZE; i++)
		Stag_V[i]= ( V[i+1] + V[i] ) / 2;

}
//
//Calculates the advection term in one dimension
//Notes:
//The Stag_V will include staggered ghost points
//Param is either temp or species with ghost points
//The Ghost points are not changed or calculated, so we march from 1 <= i < SZ-1
void Advection(realtype * Stag_V, realtype * PARAMD, realtype * OUT, int SIZE)
{
	for( int i = 1; i < SIZE-1; i++)
		OUT[i-1]= Stag_V[i] * PARAMD [i];
}
