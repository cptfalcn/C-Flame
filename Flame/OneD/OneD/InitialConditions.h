//==================================
//Header file for initial conditions
//==================================

#include <stdio.h>
#include <math.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.h>
#include <iostream>

void	IntConKappa(realtype *, realtype);
void	IntConHydro(realtype *, int);
void	IntConGri30(realtype *, int);
