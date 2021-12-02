#include <iostream>
#include <fstream>
#include <iomanip>
#include <sundials/sundials_types.h>
#include <sundials/sundials_nvector.h>
#include "Epic.h"

using namespace std;

void PrintExpData(realtype *, int, realtype);
void PrintFromPtr(realtype *, int);
void PrintDataToFile(ofstream &, realtype *, int, realtype, string, string, realtype);
void PrintExpParam(realtype, realtype, realtype, realtype, realtype, realtype, realtype, realtype, string);
void PrintSuperVector(realtype *, int ,int, string);
void PrintProfiling(IntegratorStats *, int, string, string);
void PrintMethod(string, string);
void PrintPreRun(realtype StepSize, realtype Delx, int Steps, realtype KryTol, realtype AbsTol,
			realtype RelTol, string Method, int num_pts, string Bar);

