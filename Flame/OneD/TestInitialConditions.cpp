//Testing the initial conditions
#include "InitialConditions.h"
#include "Derivatives.h"
#include "Flow.h"
#define BAR "===================="
using namespace std;

int main(void)
{
	int SIZE = 		100;
	N_Vector TEST =		N_VNew_Serial(SIZE*SIZE);
	realtype *DATA =	NV_DATA_S(TEST);
	PrintVec(SIZE, DATA, "Initialization");

	IntConKappa(DATA, 1e-2);
	PrintVec(SIZE, DATA, "Kapila IC's");

	IntConHydro(DATA, 1);//Sample 1
	PrintVec(SIZE, DATA, "Hydro IC's");

	IntConGri30(DATA, 1);//Sample 1
	PrintVec(SIZE, DATA, "Gri30 IC's");

	N_VDestroy_Serial(TEST);
	return 0;

}
