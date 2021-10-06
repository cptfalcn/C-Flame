#include "Epi3_KIOPS.h"
#include "Epic.h"
#include <memory>
#define BAR "===================="

Epi2_KIOPS *    	CreateEPI2Integrator(CVRhsFn, CVSpilsJacTimesVecFn, void *, int, N_Vector,
		 		const int, int);
Epi3_KIOPS *    	CreateEPI3Integrator(CVRhsFn, CVSpilsJacTimesVecFn, void *, int, N_Vector,
				 const int, int);
EpiRK4SC_KIOPS* 	CreateEPIRK4SCIntegrator(CVRhsFn, CVSpilsJacTimesVecFn, void *, int, N_Vector,
				 const int, int);
void * 			CreateCVODE(CVRhsFn, CVSpilsJacTimesVecFn, void *, SUNMatrix, SUNLinearSolver,
				int, N_Vector, realtype, realtype, realtype, int);

void *			CreateCVODEKrylov(CVRhsFn, CVSpilsJacTimesVecFn, void *, SUNMatrix, SUNLinearSolver,
                                SUNNonlinearSolver, int, N_Vector, realtype, realtype, realtype, int);

IntegratorStats *	Integrate(int, realtype, realtype, realtype, int, N_Vector, realtype,
				int [], void *, void *, string, Epi2_KIOPS *, Epi3_KIOPS *,
				EpiRK4SC_KIOPS *, IntegratorStats *);

void * 			SelectIntegrator(CVRhsFn, CVSpilsJacTimesVecFn, void *,
				int, N_Vector, const int, int, string, SUNMatrix, SUNLinearSolver,
				realtype, realtype, realtype);
