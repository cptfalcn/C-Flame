//========================================
//Header file for the Derivative functions
//========================================

#include <stdio.h>
#include <math.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.h>
#include <iostream>
#include <chrono>

using namespace std;

void DStaggered (realtype *, realtype * , realtype , int);
void DStaggeredDown (realtype *, realtype *, realtype, int);
void DDCentered	(realtype *, realtype * , realtype , int);
void DCentered	(realtype *, realtype * , realtype , int);
void MatVec	(int, realtype *, N_Vector, N_Vector, realtype *);
void PrintVec(int, realtype *, string);
