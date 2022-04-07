/* 
 * ===========================================================================================
 * 
 * This file is header file for EpiRK3.cpp. 
 * All constant declarations are contained in this file.
 * 
 * -------------------------------------------------------------------------------------------
 * 
 * Last revision on 02/08/2016 by Ilija Jegdic
 * 
 * ===========================================================================================
 */

#include <cvode/cvode.h>
#include <cvode/cvode_spils.h>
#include "Kiops.h"
#include "IntegratorStats.h"
#include "EpicSundials.h"
#include "EpicConst.h"
#include "EpicTypes.h"

class Epi3SC_KIOPS
{
public:
    Epi3SC_KIOPS(CVRhsFn f, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
    Epi3SC_KIOPS(CVRhsFn f, CVSpilsJacTimesVecFn jtv, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
    ~Epi3SC_KIOPS();
    //IntegratorStats *Integrate(realtype h, realtype t0, realtype tFinal, const int numBands, N_Vector y, realtype krylovTol, int startingBasisSizes[]);
    //IntegratorStats *Integrate(realtype h, realtype t0, realtype tFinal, const int numBands, N_Vector, realtype krylovTol, int startingBasisSizes[]);
    IntegratorStats* Integrate(const realtype hStart, const realtype hMax, const realtype absTol,
			const realtype relTol, const realtype t0, const realtype tFinal, const int numBands,
			int basisSizes[], N_Vector y);
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
    Epi3SC_KIOPS(const Epi3SC_KIOPS &);  // no implementation
    Epi3SC_KIOPS & operator=(const Epi3SC_KIOPS &);  // no implementation
	void Clean(N_Vector y, int Length);
};
