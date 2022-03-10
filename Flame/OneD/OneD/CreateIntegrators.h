#include "Epi3_KIOPS.h"
#include "Epic.h"
#include "IntegratorStats.h"
#include "EpicSundials.h"
#include "EpicConst.h"
#include "EpicTypes.h"
#include <memory>
#define BAR "===================="

Epi2_KIOPS *    	CreateEPI2Integrator(CVRhsFn, CVSpilsJacTimesVecFn, void *, int, N_Vector,
		 		const int, int);
Epi3_KIOPS *    	CreateEPI3Integrator(CVRhsFn, CVSpilsJacTimesVecFn, void *, int, N_Vector,
				 const int, int);
EpiRK4SC_KIOPS* 	CreateEPIRK4SCIntegrator(CVRhsFn, CVSpilsJacTimesVecFn, void *, int, N_Vector,
				 const int, int);

void * 			CreateCVODE(CVRhsFn, CVSpilsJacTimesVecFn, CVLsJacFn, void *, SUNMatrix,
				SUNLinearSolver, SUNNonlinearSolver, int, N_Vector, realtype, realtype,
				realtype, int);

void *                  CreateCVODEAdp(CVRhsFn, CVSpilsJacTimesVecFn, CVLsJacFn, void *, SUNMatrix,
                                SUNLinearSolver, SUNNonlinearSolver, int, N_Vector, realtype, realtype,
                                realtype);

void *			CreateCVODEKrylov(CVRhsFn, CVSpilsJacTimesVecFn, void *, SUNMatrix, SUNLinearSolver,
                                SUNNonlinearSolver, int, N_Vector, realtype, realtype, realtype, int);

IntegratorStats *	Integrate(int, realtype, realtype, realtype, int, N_Vector, realtype,
				int [], void *, void *, string, Epi2_KIOPS *, Epi3_KIOPS *,
				EpiRK4SC_KIOPS *, IntegratorStats *);

//OneD stuff, outdated and being replaced.

void *  		SelectIntegrator(CVRhsFn, CVSpilsJacTimesVecFn, CVLsJacFn, void *, void *,
				int, N_Vector, const int, int, string, SUNMatrix, SUNLinearSolver,
				SUNNonlinearSolver, realtype, realtype, realtype);

IntegratorStats * 	IntegrateWrapper(int, realtype, realtype, realtype, int, N_Vector, realtype, int [],
						void * , void *, void *, string, IntegratorStats *, CVLsJacFn);


class Epi3VChem_KIOPS
{
public:
    	 Epi3VChem_KIOPS(CVRhsFn f, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
    	 Epi3VChem_KIOPS(CVRhsFn f, CVSpilsJacTimesVecFn jtv, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
    	~Epi3VChem_KIOPS();
    	//IntegratorStats *Integrate(realtype h, realtype t0, realtype tFinal, const int numBands, N_Vector y, realtype krylovTol, int startingBasisSizes[]);
    	//IntegratorStats *Integrate(realtype h, realtype t0, realtype tFinal, const int numBands, N_Vector, realtype krylovTol, int startingBasisSizes[]);
 	IntegratorStats* Integrate(const realtype hStart, const realtype hMax, const realtype absTol,
                        const realtype relTol, const realtype t0, const realtype tFinal,
			const int numBands, int basisSizes[], N_Vector y);
	void TestMatMult();
private:
    	CVRhsFn f;
    	CVSpilsJacTimesVecFn jtv;
    	EPICNumJacDelta delta;
    	void *userData;
    	Kiops *krylov;
    	static const int NumProjectionsPerStep = 1;
    	static const int MaxPhiOrder1 = 2;
    	const int NEQ;
    	IntegratorStats *integratorStats;

    	N_Vector fy;
        N_Vector Y1;
    	N_Vector fY1;
    	N_Vector hfy;
    	N_Vector hb1;
    	N_Vector hb2;
    	N_Vector r1;
        N_Vector r2;
        N_Vector r3;
        N_Vector Remainder;
    	N_Vector tmpVec;       // used as a temporary work vector by the jtv method
    	N_Vector zeroVec;      // zero-ed out vector, used by krylov
        realtype hMax;
        N_Vector Scratch1;

    	// Disallow copying by "poisoning" the copy constructor and assignment operator,
    	// i.e. declare private but provide no implementation.
    	Epi3VChem_KIOPS(const Epi3VChem_KIOPS &);  // no implementation
    	Epi3VChem_KIOPS & operator=(const Epi3VChem_KIOPS &);  // no implementation
        void Clean(N_Vector y, int Length);
	void CheckNanBlowUp(N_Vector, int);
	void MatMult(N_Vector R, N_Vector A, N_Vector B, int row, int col);
	void Phi2(N_Vector Jac, realtype h, N_Vector Result);
	int TrackProgress(realtype FinalTime, realtype TNext, realtype PercentDone, int ProgressDots);
};
