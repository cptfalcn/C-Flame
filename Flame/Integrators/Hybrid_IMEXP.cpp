/*
 * ===========================================================================================
 * 
 * This file contains implementation of the Hybrid IMEXP method. 
 * 
 * ===========================================================================================
 */

#include "Hybrid_IMEXP.h"
#include <math.h>
#include <vector>
#include <stdlib.h>
#include "EpicConst.h"

using std::vector; 

/*
 * ===========================================================================================
 * 
 * Declares all variables
 * 
 * Inputs:
 * f                RHS function from the test problem 
 * jtv              Jacobian times vector, i.e., Jtv function from the test problem. 
 * userData         structure UserData defined in the test problem. 
 * maxKrylovIters   maximum number of the Krylov iterations
 * tmpl             initial condition; it will be used to form other vectors necessary 
 *                  for the implementation
 * vecLength        
 * NEQ              number of processors multiplied by the number of nods in a subgrid 
 *                  (defined in header file of the test problem)
 * 
 * NOTE:  It's crucial that the tmpl vector have no NAN values in it.
 * 
 * ===========================================================================================
 */

Hybrid_IMEXP::Hybrid_IMEXP(CVRhsFn f1, CVSpilsJacTimesVecFn j1tv, CVRhsFn f2, CVSpilsJacTimesVecFn j2tv, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    this->f1 = f1;
    this->f2 = f2;
    this->j1tv = j1tv;
    this->j2tv = j2tv;
    this->userData = userData;

    krylov = new Kiops(2, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(1);

    f1y = N_VClone(tmpl);
    half_hf1y = N_VClone(tmpl);
    f2y = N_VClone(tmpl);
    hf2y = N_VClone(tmpl);
    r1 = N_VClone(tmpl);
    tmpVec = N_VClone(tmpl);
    zeroVec = N_VClone(tmpl);
    // Zero out zeroVec by using the trusted vector tmpl.
    // If zeroVec has random data in it, it may have NAN values.
    // zeroVec = 0.0 * zeroVec may accidentally be the calculation
    // 0.0 * NAN = NAN, and zeroVec will not end up zero.
    // Instead zeroVec = 0.0 * templ, since templ is trusted
    // to have no NAN values in it.
    N_VConst(0, f1y);
    N_VConst(0, half_hf1y);
    N_VConst(0, f2y);
    N_VConst(0, hf2y);
    N_VConst(0, r1);
    N_VConst(0, zeroVec);

    LS = SUNLinSol_SPGMR(tmpl, PREC_NONE, 20);
}

/*
 * ===========================================================================================
 * 
 * Declares all variables. This routine is invoked if Jacobian is not provided.
 * 
 * Inputs:
 * f                RHS function from the test problem 
 * userData         structure UserData defined in the test problem. 
 * maxKrylovIters   maximum number of the Krylov iterations
 * tmpl             initial condition; it will be used to form other vectors necessary 
 *                  for the implementation
 * vecLength        
 * NEQ              number of processors multiplied by the number of nods in a subgrid 
 *                  (defined in header file of the test problem)
 * 
 * NOTE:  It's crucial that the tmpl vector have no NAN values in it.
 * 
 * ===========================================================================================
 */

Hybrid_IMEXP::Hybrid_IMEXP(CVRhsFn f1, CVRhsFn f2, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    this->f1 = f1;
    this->f2 = f2;
    this->userData = userData;
    this->j1tv = nullptr;
    this->j2tv = nullptr;
    this->delta = nullptr;
    
    krylov = new Kiops(2, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(1);

    f1y = N_VClone(tmpl);
    half_hf1y = N_VClone(tmpl);
    f2y = N_VClone(tmpl);
    hf2y = N_VClone(tmpl);
    r1 = N_VClone(tmpl);
    tmpVec = N_VClone(tmpl);
    zeroVec = N_VClone(tmpl);
    // Zero out zeroVec by using the trusted vector tmpl.
    // If zeroVec has random data in it, it may have NAN values.
    // zeroVec = 0.0 * zeroVec may accidentally be the calculation
    // 0.0 * NAN = NAN, and zeroVec will not end up zero.
    // Instead zeroVec = 0.0 * templ, since templ is trusted
    // to have no NAN values in it.
    N_VConst(0, f1y);
    N_VConst(0, half_hf1y);
    N_VConst(0, f2y);
    N_VConst(0, hf2y);
    N_VConst(0, r1);
    N_VConst(0, zeroVec);

    LS = SUNLinSol_SPGMR(tmpl, PREC_NONE, 20);
}

/*
 * ===========================================================================================
 * 
 * Declares all variables. This routine is invoked if Jacobian is not provided.
 * 
 * Inputs:
 * f                RHS function from the test problem 
 * delta            Delta function used in numerical Jacobian
 * userData         structure UserData defined in the test problem. 
 * maxKrylovIters   maximum number of the Krylov iterations
 * tmpl             initial condition; it will be used to form other vectors necessary 
 *                  for the implementation
 * vecLength        
 * NEQ              number of processors multiplied by the number of nods in a subgrid 
 *                  (defined in header file of the test problem)
 * 
 * NOTE:  It's crucial that the tmpl vector have no NAN values in it.
 * 
 * ===========================================================================================
 */

Hybrid_IMEXP::Hybrid_IMEXP(CVRhsFn f1, CVRhsFn f2, EPICNumJacDelta delta, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    this->f1 = f1;
    this->f2 = f2;
    this->userData = userData;
    this->j1tv = nullptr;
    this->j2tv = nullptr;
    this->delta = delta;
    
    krylov = new Kiops(2, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(1);

    f1y = N_VClone(tmpl);
    half_hf1y = N_VClone(tmpl);
    f2y = N_VClone(tmpl);
    hf2y = N_VClone(tmpl);
    r1 = N_VClone(tmpl);
    tmpVec = N_VClone(tmpl);
    zeroVec = N_VClone(tmpl);
    // Zero out zeroVec by using the trusted vector tmpl.
    // If zeroVec has random data in it, it may have NAN values.
    // zeroVec = 0.0 * zeroVec may accidentally be the calculation
    // 0.0 * NAN = NAN, and zeroVec will not end up zero.
    // Instead zeroVec = 0.0 * templ, since templ is trusted
    // to have no NAN values in it.
    N_VConst(0, f1y);
    N_VConst(0, half_hf1y);
    N_VConst(0, f2y);
    N_VConst(0, hf2y);
    N_VConst(0, r1);
    N_VConst(0, zeroVec);

    LS = SUNLinSol_SPGMR(tmpl, PREC_NONE, 20);
}

Hybrid_IMEXP::~Hybrid_IMEXP()
{
    delete krylov;
    delete integratorStats;

    N_VDestroy(f1y);
    N_VDestroy(half_hf1y);
    N_VDestroy(f2y);
    N_VDestroy(hf2y);
    N_VDestroy(r1);
    N_VDestroy(tmpVec);
    N_VDestroy(zeroVec);
}

struct data_SUNLinSol
{
    JTimesV* jtimesv;
    realtype h;
};


int JTimesV_SUNLinSol(void* A_data, N_Vector v, N_Vector z) {
    data_SUNLinSol* data = static_cast<data_SUNLinSol*>(A_data);
    data->jtimesv->ComputeJv(v, z);
    N_VLinearSum(1.0, v, - data->h/2.0, z, z); 

    return 0;
}



double sumVector (N_Vector v) {
    N_Vector ones = N_VClone(v);
    N_VConst(1.,ones);

    return N_VDotProd(ones, v);
}

/*
 * ===========================================================================================
 * 
 * Function Integrate
 * 
 * Inputs:
 * h              time step size
 * t0             starting time
 * tFinal         ending time
 * numBands
 * y              initial condition
 * krylovTol      tolerance for Krylov Iteration
 * basisSizes[]   maximum basis size allowed for the Arnoldi iteration
 * 
 * Output:
 * Statistics
 * 
 * ===========================================================================================
 */
IntegratorStats *Hybrid_IMEXP::Integrate(realtype h, realtype t0, realtype tFinal, const int numBands, N_Vector y, realtype krylovTol, int basisSizes[])
{
    if( h < ZERO )
    {
        printf("Time step h is to small. \n");
        exit(EXIT_FAILURE);
    }
    if( tFinal < t0 )
    {
        printf("Starting time is larger the end time. \n");
        exit(EXIT_FAILURE);
    }
    
    realtype t = t0, hNew;
    vector<realtype> timeSteps;
    
    if( tFinal <= t0 + h )
    {
        hNew = tFinal - t;
        timeSteps.push_back(hNew);
        printf("Since tFinal - t is smaller then h, h is reduced to tFinal - t0; h = %e \n", hNew);
    }
    else
    {
        while( t + h < tFinal)
        {
            t += h;
            timeSteps.push_back(h);
        }
        
        if(tFinal - t < ZERO)
            timeSteps[timeSteps.size() - 1] = tFinal - t + h;
        else
            timeSteps.push_back(tFinal - t);
        t = t0;
        for(auto &k : timeSteps)
          k = h;
    }
    
    data_SUNLinSol A_data;
        
    // Main integration loop.
    for( auto &k : timeSteps )
    {
        h = k;
        // f is RHS function from the problem; y = u_n
        f1(t, y, f1y, userData); 
        N_VScale(h/2., f1y, half_hf1y);

        f2(t, y, f2y, userData); 
        N_VScale(h, f2y, hf2y);
        
        JTimesV j1timesv(j1tv, f1, delta, t, y, f1y, userData, tmpVec);
        JTimesV j2timesv(j2tv, f2, delta, t, y, f2y, userData, tmpVec);

        // Exponential part
        N_Vector stage1InputVecs[] = {half_hf1y, hf2y};
        N_Vector stage1OutputVecs[] = {r1}; 
        const double tau[] = {1.};
        krylov->Compute(2, stage1InputVecs, tau, 1, stage1OutputVecs, &j2timesv, h, krylovTol, basisSizes[0], &integratorStats->krylovStats[0]);
            
        N_VLinearSum(1.0, r1, 1.0, half_hf1y, r1); 

        // Linear solve
        A_data.h = h;
        A_data.jtimesv = &j1timesv;
        SUNLinSolSetATimes(LS, &A_data, &JTimesV_SUNLinSol);
        SUNLinSolInitialize(LS);
        SUNLinSolSolve(LS, NULL, tmpVec, r1, 1e-14);
          
        N_VLinearSum(1.0, y, 1.0, tmpVec, y); 
          
        t = t + h;
        integratorStats->Step();

    }

    return integratorStats;
}
