//==================================
//Header file for initial conditions
//==================================

#include <stdio.h>
#include <math.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.h>
#include <iostream>
#include <fstream>
#include <iomanip>


void 	SetIntCons(int, int, realtype*);
void	IntConKappa(realtype *, realtype);
void	IntConHydro(realtype *, int);
void	IntConGri30(realtype *, int);
void 	SetSuperInitCons(realtype *, realtype *, int, int);

//Testing
void	TestingInitCons(int, int, int, int, realtype, realtype *, realtype *);
void	SinWave(realtype *, int, int, realtype);
void	SinWave2(realtype*, int, int, realtype);
void	Constant(realtype, int);
void	ConditionsFromFile(realtype *, std :: string);
void	GaussWave(realtype *, int , int, realtype);
