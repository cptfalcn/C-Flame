/*
 * ===========================================================================================
 *
 * This file contains implementation of the Epi2 method. 
 * All constant declarations are declared in header file.
 * 
 * u_{n+1} = u_n + phi_1(h J_n) h f(u_n) 
 * 
 * -------------------------------------------------------------------------------------------
 * 
 * Last revision on 07/05/2016 by Ilija Jegdic
 * 
 * ===========================================================================================
 */

#include "Epi2_KIOPS.h"
#include <math.h>
#include <vector>
#include <stdlib.h>
#include "EpicConst.h"
#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.h>

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

Epi2_KIOPS::Epi2_KIOPS(CVRhsFn f, CVSpilsJacTimesVecFn jtv, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    this->f = f;
    this->jtv = jtv;
    this->userData = userData;

    krylov = new Kiops(MaxPhiOrder1 + 1, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(NumProjectionsPerStep);

    fy = N_VClone(tmpl);
    hfy = N_VClone(tmpl);
    hb1 = N_VClone(tmpl);
    hb2 = N_VClone(tmpl);
    r1 = N_VClone(tmpl);
    tmpVec = N_VClone(tmpl);
    zeroVec = N_VClone(tmpl);
    // Zero out zeroVec by using the trusted vector tmpl.
    // If zeroVec has random data in it, it may have NAN values.
    // zeroVec = 0.0 * zeroVec may accidentally be the calculation
    // 0.0 * NAN = NAN, and zeroVec will not end up zero.
    // Instead zeroVec = 0.0 * templ, since templ is trusted
    // to have no NAN values in it.
    N_VConst(0, fy);
    N_VConst(0, hfy);
    N_VConst(0, hb1);
    N_VConst(0, hb2);
    N_VConst(0, r1);
    N_VConst(0, zeroVec);
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

Epi2_KIOPS::Epi2_KIOPS(CVRhsFn f, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    this->f = f;
    this->userData = userData;
    this->jtv = nullptr;
    this->delta = nullptr;
    
    krylov = new Kiops(MaxPhiOrder1 + 1, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(NumProjectionsPerStep);

    fy = N_VClone(tmpl);
    hfy = N_VClone(tmpl);
    hb1 = N_VClone(tmpl);
    hb2 = N_VClone(tmpl);
    r1 = N_VClone(tmpl);
    tmpVec = N_VClone(tmpl);
    zeroVec = N_VClone(tmpl);
    // Zero out zeroVec by using the trusted vector tmpl.
    // If zeroVec has random data in it, it may have NAN values.
    // zeroVec = 0.0 * zeroVec may accidentally be the calculation
    // 0.0 * NAN = NAN, and zeroVec will not end up zero.
    // Instead zeroVec = 0.0 * templ, since templ is trusted
    // to have no NAN values in it.
    N_VConst(0, fy);
    N_VConst(0, hfy);
    N_VConst(0, hb1);
    N_VConst(0, hb2);
    N_VConst(0, r1);
    N_VConst(0, zeroVec);

}


Epi2_KIOPS::~Epi2_KIOPS()
{
    delete krylov;
    delete integratorStats;
    
    N_VDestroy(fy);
    N_VDestroy(hfy);
    N_VDestroy(hb1);
    N_VDestroy(hb2);
    N_VDestroy(r1);
    N_VDestroy(tmpVec);
    N_VDestroy(zeroVec);

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

IntegratorStats *Epi2_KIOPS::Integrate(realtype h, realtype t0, realtype tFinal, const int numBands, N_Vector y, realtype krylovTol, int basisSizes[])
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
    
    if( tFinal < t0 + h )
    {
        hNew = tFinal - t;
        timeSteps.push_back(hNew);
        //printf("Since tFinal - t is smaller then h, h is reduced to tFinal - t0; h = %e \n", hNew);
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
    
    // Main integration loop.
    for( auto &k : timeSteps )
    {
        h = k;
        // f is RHS function from the problem; y = u_n
        f(t, y, fy, userData); // f(t, y) = fy
        N_VScale(h, fy, hfy);
        
        JTimesV jtimesv(jtv, f, delta, t, y, fy, userData, tmpVec);
        
        N_Vector stage1InputVecs[] = {zeroVec, hfy};
        N_Vector stage1OutputVecs[] = {r1}; 
        realtype timePts[] = {1.0};
        krylov->Compute(2, stage1InputVecs, timePts, 1, stage1OutputVecs, &jtimesv, h, krylovTol, basisSizes[0], &integratorStats->krylovStats[0]);

        N_VLinearSum(1.0, y, 1.0, r1, y); // r1 = y + r1 = y + phi( ) * hfy
	//This block added by me to make the problem stable
	//realtype *data_y=NV_DATA_S(y);
	//int length = N_VGetLength(y);
	//This was added by me, I should changed this back later.
	/*
	for (int i=0; i<length; i++){
		if(data_y[i]<0){
			data_y[i]=0;
		}
	}
	*/
	//End block
         t = t + h;
        integratorStats->Step();
    }

    return integratorStats;
}
