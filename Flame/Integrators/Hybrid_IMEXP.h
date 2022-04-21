/* 
 * ===========================================================================================
 * 
 * This file is header file for Hybrid_IMEXP.cpp. 
 * All constant declarations are contained in this file.
 * 
 * ===========================================================================================
 */

#include <cvode/cvode.h>
#include <cvode/cvode_spils.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include "Kiops.h"
#include "IntegratorStats.h"
#include "EpicSundials.h"
#include "EpicConst.h"
#include "EpicTypes.h"

class Hybrid_IMEXP
{
public:
    Hybrid_IMEXP(CVRhsFn f1, CVRhsFn f2, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
    Hybrid_IMEXP(CVRhsFn f1, CVSpilsJacTimesVecFn j1tv, CVRhsFn f2, CVSpilsJacTimesVecFn j2tv, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
    Hybrid_IMEXP(CVRhsFn f1, CVRhsFn f2, EPICNumJacDelta delta, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
    ~Hybrid_IMEXP();
    IntegratorStats *Integrate(realtype h, realtype t0, realtype tFinal, const int numBands, N_Vector y, realtype krylovTol, int startingBasisSizes[]);

private:
    CVRhsFn f1;
    CVRhsFn f2;
    CVSpilsJacTimesVecFn j1tv;
    CVSpilsJacTimesVecFn j2tv;
    EPICNumJacDelta delta;
    void *userData;
    Kiops *krylov;
    static const int MaxPhiOrder1 = 1;
    const int NEQ;
    IntegratorStats *integratorStats;

    SUNLinearSolver LS;

    N_Vector f1y;
    N_Vector half_hf1y;
    N_Vector f2y;
    N_Vector hf2y;
    N_Vector r1;
    N_Vector tmpVec;       // used as a temporary work vector by the jtv method
    N_Vector zeroVec;      // zero-ed out vector, used by krylov

    // Disallow copying by "poisoning" the copy constructor and assignment operator,
    // i.e. declare private but provide no implementation.
    Hybrid_IMEXP(const Hybrid_IMEXP &);  // no implementation
    Hybrid_IMEXP & operator=(const Hybrid_IMEXP &);  // no implementation
};
