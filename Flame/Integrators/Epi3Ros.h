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

class Epi3Ros_KIOPS
{
public:
    Epi3Ros_KIOPS(CVRhsFn f, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
    Epi3Ros_KIOPS(CVRhsFn f, CVSpilsJacTimesVecFn jtv, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
    ~Epi3Ros_KIOPS();
    //IntegratorStats *Integrate(realtype h, realtype t0, realtype tFinal, const int numBands, N_Vector y, realtype krylovTol, int startingBasisSizes[]);
    IntegratorStats *Integrate(realtype h, realtype t0, realtype tFinal, const int numBands, N_Vector, realtype krylovTol, int startingBasisSizes[]);

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

    // Disallow copying by "poisoning" the copy constructor and assignment operator,
    // i.e. declare private but provide no implementation.
    Epi3Ros_KIOPS(const Epi3Ros_KIOPS &);  // no implementation
    Epi3Ros_KIOPS & operator=(const Epi3Ros_KIOPS &);  // no implementation
};
