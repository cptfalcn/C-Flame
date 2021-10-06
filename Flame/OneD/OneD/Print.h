#include <iostream>
#include <fstream>
#include <iomanip>
#include <sundials/sundials_types.h>

using namespace std;

void PrintExpData(realtype *, int, realtype);
void PrintFromPtr(realtype *, int);
void PrintDataToFile(ofstream &, realtype *, int, realtype);
void PrintExpParam(realtype, realtype, realtype, realtype, realtype, realtype, realtype, realtype);
