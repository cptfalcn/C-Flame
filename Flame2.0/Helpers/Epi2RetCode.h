/* 
 * ===========================================================================================
 * 
 * This file is header file for EpiRK2.cpp. 
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
#include "KiopsRetCode.h"
#include "IntegratorStats.h"
#include "EpicSundials.h"
#include "EpicConst.h"
#include "EpicTypes.h"

class Epi2RetCode
{
public:
    Epi2RetCode(CVRhsFn f, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
    Epi2RetCode(CVRhsFn f, CVSpilsJacTimesVecFn jtv, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
    Epi2RetCode(CVRhsFn f, EPICNumJacDelta delta, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
    ~Epi2RetCode();
    IntegratorStats *Integrate(realtype h, realtype t0, realtype tFinal, const int numBands, N_Vector y, realtype krylovTol, int startingBasisSizes[]);

private:
    CVRhsFn f;
    CVSpilsJacTimesVecFn jtv;
    EPICNumJacDelta delta;
    void *userData;
    KiopsRetCode * krylov;
    static const int NumProjectionsPerStep = 1;
    static const int MaxPhiOrder1 = 2;
    const int NEQ;
    IntegratorStats *integratorStats;

    N_Vector fy;
    N_Vector hfy;
    N_Vector hb1;
    N_Vector hb2;
    N_Vector r1;
    N_Vector tmpVec;       // used as a temporary work vector by the jtv method
    N_Vector zeroVec;      // zero-ed out vector, used by krylov

    // Disallow copying by "poisoning" the copy constructor and assignment operator,
    // i.e. declare private but provide no implementation.
    Epi2RetCode(const Epi2RetCode &);  // no implementation
    Epi2RetCode & operator=(const Epi2RetCode &);  // no implementation
};
