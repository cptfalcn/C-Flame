/*
 * ===========================================================================================
 * 
 * This file contains implementation of the EpiRK5SC method. 
 * EpiRK5SC method is method with constant time step and uses adaptive Krylov
 * method to compute Krylov iterations. 
 * Unlike it is suggested in a paper, this method is implemented such that U_2 and U_3 are 
 * computed vertically, while u_{n+1} is computed horizontaly. 
 * 5th order Runge-Kutta exponential propagation integrator.
 * EpiRK5SC method is part of Epic Package. 
 * All constant declarations are declared in header file.
 * 
 * U_2     = u_n + 1/2 phi_1(1/2 h J_n) h f(u_n)
 * U_3     = u_n + 9/10 phi_1(9/10 h J_n) h f(u_n)
 *               + (27/25 phi_3(1/2 h J_n) + 729/125 phi_3(9/10 h J_n)) h r(U_2)
 * u_{n+1} = u_n + phi_1(h J_n) h f(u_n) + (18 phi_3(h J_n) + 60 phi_4(h J_n)) h r(U_2)
 *               + (-250/81 phi_3(h J_n) + 500/27 phi_4 (h J_n)) h r(U_3)
 * 
 * -------------------------------------------------------------------------------------------
 * 
 * Last revision on 07/05/2016 by Ilija Jegdic
 * 
 * ===========================================================================================
 */

#include "EpiRK5SC.h"
#include <math.h>
#include <vector>
#include <stdlib.h>

using std::vector; 
using namespace EpiRK5SCNamespace;

/*
 * ===========================================================================================
 * 
 * Declares all variables
 * 
 * Inputs:
 * f                RHS function from the test problem 
 * jtv              Jacobian, i.e., Jtv function from the test problem
 * userData         structure State defined in the test problem. 
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

EpiRK5SC::EpiRK5SC(CVRhsFn f, CVSpilsJacTimesVecFn jtv, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    this->f = f;
    this->jtv = jtv;
    this->userData = userData;

    krylov1 = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder1 + 1, maxKrylovIters, tmpl, NEQ);
    krylov2 = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder2 + 1, maxKrylovIters, tmpl, NEQ);
    krylov3 = new AdaptiveKrylovPhi(MaxPhiOrder2 + 1, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(NumProjectionsPerStep);

    fy = N_VClone(tmpl);
    hfy = N_VClone(tmpl);
    hb1 = N_VClone(tmpl);
    hb2 = N_VClone(tmpl);
    r2 = N_VClone(tmpl);
    r3 = N_VClone(tmpl);
    r32 = N_VClone(tmpl);
    r33 = N_VClone(tmpl);
    scratchVec1 = N_VClone(tmpl);
    scratchVec2 = N_VClone(tmpl);
    scratchVec3 = N_VClone(tmpl);
    tmpVec = N_VClone(tmpl);
    zeroVec = N_VClone(tmpl);
    // Zero out zeroVec by using the trusted vector tmpl.
    // If zeroVec has random data in it, it may have NAN values.
    // zeroVec = 0.0 * zeroVec may accidentally be the calculation
    // 0.0 * NAN = NAN, and zeroVec will not end up zero.
    // Instead zeroVec = 0.0 * templ, since templ is trusted
    // to have no NAN values in it.
    N_VScale(0.0, tmpl, zeroVec);
    N_VConst(0, fy);
    N_VConst(0, hfy);
    N_VConst(0, hb1);
    N_VConst(0, hb2);
    N_VConst(0, r2);
    N_VConst(0, r3);
    N_VConst(0, r32);
    N_VConst(0, r33);
    N_VConst(0, scratchVec1);
    N_VConst(0, scratchVec2);
    N_VConst(0, scratchVec3);
}

/*
 * ===========================================================================================
 * 
 * Declares all variables
 * 
 * Inputs:
 * f                RHS function from the test problem 
 * userData         structure State defined in the test problem. 
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

EpiRK5SC::EpiRK5SC(CVRhsFn f, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    this->f = f;
    this->jtv = nullptr;
    this->delta = nullptr;
    this->userData = userData;

    krylov1 = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder1 + 1, maxKrylovIters, tmpl, NEQ);
    krylov2 = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder2 + 1, maxKrylovIters, tmpl, NEQ);
    krylov3 = new AdaptiveKrylovPhi(MaxPhiOrder2 + 1, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(NumProjectionsPerStep);

    fy = N_VClone(tmpl);
    hfy = N_VClone(tmpl);
    hb1 = N_VClone(tmpl);
    hb2 = N_VClone(tmpl);
    r2 = N_VClone(tmpl);
    r3 = N_VClone(tmpl);
    r32 = N_VClone(tmpl);
    r33 = N_VClone(tmpl);
    scratchVec1 = N_VClone(tmpl);
    scratchVec2 = N_VClone(tmpl);
    scratchVec3 = N_VClone(tmpl);
    tmpVec = N_VClone(tmpl);
    zeroVec = N_VClone(tmpl);
    // Zero out zeroVec by using the trusted vector tmpl.
    // If zeroVec has random data in it, it may have NAN values.
    // zeroVec = 0.0 * zeroVec may accidentally be the calculation
    // 0.0 * NAN = NAN, and zeroVec will not end up zero.
    // Instead zeroVec = 0.0 * templ, since templ is trusted
    // to have no NAN values in it.
    N_VScale(0.0, tmpl, zeroVec);
    N_VConst(0, fy);
    N_VConst(0, hfy);
    N_VConst(0, hb1);
    N_VConst(0, hb2);
    N_VConst(0, r2);
    N_VConst(0, r3);
    N_VConst(0, r32);
    N_VConst(0, r33);
    N_VConst(0, scratchVec1);
    N_VConst(0, scratchVec2);
    N_VConst(0, scratchVec3);
}

/*
 * ===========================================================================================
 * 
 * Declares all variables
 * 
 * Inputs:
 * f                RHS function from the test problem 
 * delta            Delta function used in numerical Jacobian
 * userData         structure State defined in the test problem. 
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

EpiRK5SC::EpiRK5SC(CVRhsFn f, EPICNumJacDelta delta, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    this->f = f;
    this->jtv = nullptr;
    this->delta = delta;
    this->userData = userData;

    krylov1 = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder1 + 1, maxKrylovIters, tmpl, NEQ);
    krylov2 = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder2 + 1, maxKrylovIters, tmpl, NEQ);
    krylov3 = new AdaptiveKrylovPhi(MaxPhiOrder2 + 1, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(NumProjectionsPerStep);

    fy = N_VClone(tmpl);
    hfy = N_VClone(tmpl);
    hb1 = N_VClone(tmpl);
    hb2 = N_VClone(tmpl);
    r2 = N_VClone(tmpl);
    r3 = N_VClone(tmpl);
    r32 = N_VClone(tmpl);
    r33 = N_VClone(tmpl);
    scratchVec1 = N_VClone(tmpl);
    scratchVec2 = N_VClone(tmpl);
    scratchVec3 = N_VClone(tmpl);
    tmpVec = N_VClone(tmpl);
    zeroVec = N_VClone(tmpl);
    // Zero out zeroVec by using the trusted vector tmpl.
    // If zeroVec has random data in it, it may have NAN values.
    // zeroVec = 0.0 * zeroVec may accidentally be the calculation
    // 0.0 * NAN = NAN, and zeroVec will not end up zero.
    // Instead zeroVec = 0.0 * templ, since templ is trusted
    // to have no NAN values in it.
    N_VScale(0.0, tmpl, zeroVec);
    N_VConst(0, fy);
    N_VConst(0, hfy);
    N_VConst(0, hb1);
    N_VConst(0, hb2);
    N_VConst(0, r2);
    N_VConst(0, r3);
    N_VConst(0, r32);
    N_VConst(0, r33);
    N_VConst(0, scratchVec1);
    N_VConst(0, scratchVec2);
    N_VConst(0, scratchVec3);
}

EpiRK5SC::~EpiRK5SC()
{
    delete krylov1;
    delete krylov2;
    delete krylov3;
    delete integratorStats;

    N_VDestroy(fy);
    N_VDestroy(hfy);
    N_VDestroy(hb1);
    N_VDestroy(hb2);
    N_VDestroy(r2);
    N_VDestroy(r3);
    N_VDestroy(r32);
    N_VDestroy(r33);
    N_VDestroy(scratchVec1);
    N_VDestroy(scratchVec2);
    N_VDestroy(scratchVec3);
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

IntegratorStats *EpiRK5SC::Integrate(realtype h, realtype t0, realtype tFinal, const int numBands, N_Vector y, realtype krylovTol, int basisSizes[])
{
    if( h < 10e-14 )
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
        //printf("Since tFinal - t is smaller then h, h is reduced to tFinal -  t0; h = %e \n", hNew);
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
        f(t, y, fy, userData);
        N_VScale(h, fy, hfy);
        JTimesV jtimesv(jtv, f, delta, t, y, fy, userData, tmpVec);

        // Stage 1.
        N_Vector stage1InputVecs[] = {zeroVec, hfy};
        N_Vector stage1OutputVecs[] = {r2, r3}; //r2 = phi(1/2 h J_n) hfy; r3 = phi(9/10 h J_n) hfy
        krylov1->Compute(numBands, Stage1NumInputVecs, stage1InputVecs, g1Times, g1NumTimes, stage1OutputVecs, &jtimesv, h, krylovTol, basisSizes[0], &integratorStats->krylovStats[0]);
	// computes phi_k(g_{i j} h J)*hfy and stores it in ri, ri is r(U_i)
	// in this case r2 = phi_1( 1/2 h J_n ) * hfy 
	// and          r3 = phi_1( 9/10 h J_n ) * hfy
	
	// States U_2 and U_3 are now being computed and results are stored in r2 and r3, respectively
	N_VLinearSum(1.0, y, a21, r2, r2); // r2 = y + a21 * r2 = y + a21 * phi( ) * hfy
        N_VLinearSum(1.0, y, a31, r3, r3); // r3 = y + a31 * r2 = y + a31 * phi( ) * hfy
	
	// Stage 3 - High-order part
	// In this part we compute residual of U_2, i.e., h*r(U_2). Result will be stored in hb1 variable
        N_VLinearSum(1.0, r2, -1.0, y, scratchVec1); // scratchVec1 = a21 * phi( ) * hfy
        jtimesv.ComputeJv(scratchVec1, scratchVec2); // scratchVec1 can be reused now, scratchVec2 is Jacobian of scratchVec1
        f(t, r2, hb1, userData); // f(t, r2) = hb1
        N_VLinearSum(h, hb1, -1.0, hfy, hb1); // hb1 = h * hb1 - hfy
        N_VLinearSum(1.0, hb1, -1.0*h, scratchVec2, hb1); // scratchVec2 can be reused now, hb1 contains residual of U_2
	// In order to compute residual we did the following
	// hb1 = hb1 - h*scratchVec2
	//     = h*hb1 - hfy - Jacobian(scratchVec1)
	//     = f(t, r1) - hfy - Jacobian(r1-y)
	//     = f(t, y + a21*phi( )*hfy) - hfy - Jacobian( a21*phi( )*hfy )*hfy
        //     = h*r(U_2) <- residual
	
	// scratchVec1 and scratchVec2 are free variables
	
	// In this part we finish conmputing state U_3
	N_Vector stage2InputVecs[] = {zeroVec, zeroVec, zeroVec, hb1};
        N_Vector stage2OutputVecs[] = {r32, r33}; // r32 = phi3(1/2 h J_n) hb1; r33 = phi(9/10 h J_n) hb1
        krylov2->Compute(numBands, Stage2NumInputVecs, stage2InputVecs, g2Times, g2NumTimes, stage2OutputVecs, &jtimesv, h, krylovTol, basisSizes[0], &integratorStats->krylovStats[1]); 
	// computes phi_k(g_{i j} h J)*hfy and stores it in ri, ri is r(U_i)
	// in this case r32 = phi_3( 1/2 h J_n ) * hb1 
	// and          r33 = phi_3( 9/10 h J_n ) * hb1
	// here verical method is used sice we don't have horizontal method for g \neq 1 at this time
	
	// state U_3
	N_VLinearSum(1.0, r3, a32, r32, r3); // r3 = r3 + a32 phi3(1/2 h J_n) hb1
	N_VLinearSum(1.0, r3, a33, r33, r3); // r3 = r3 + a33 phi(9/10 h J_n) hb1 = U_3
	
	// In this part we compute residual of U_3, i.e., h*r(U_3). Result will be stored in hb2 variable
	// We preform same caluculation we did for residual of U_2
        N_VLinearSum(1.0, r3, -1.0, y, scratchVec1); 
        jtimesv.ComputeJv(scratchVec1, scratchVec2); 
        f(t, r3, hb2, userData);
        N_VLinearSum(h, hb2, -1.0, hfy, hb2);
        N_VLinearSum(1.0, hb2, -1.0*h, scratchVec2, hb2); 
	// Similar like for hb1, hb2 is now h*r(U_3)
	
	// scratchVec1 and scratchVec2 are free variables
	
	// State u_{n+1} is now being computed
        N_VLinearSum(b2a, hb1, b2b, hb2, scratchVec1);
        N_VLinearSum(b3a, hb1, b3b, hb2, scratchVec2);
	
	// scratchVec1 and scratchVec2 contain linear combinations of residuals
	
	N_Vector stage3InputVecs[] = {zeroVec, hfy, zeroVec, scratchVec1, scratchVec2};
	krylov3->Compute(numBands, Stage3NumInputVecs, stage3InputVecs, scratchVec3, 1.0, &jtimesv, h, krylovTol, basisSizes[1], &integratorStats->krylovStats[2]);  
	// scratchVec3 now holds high-order portion, i.e., 
	// scratchVec3 = phi_1(h J_n) h f(u_n) + (18 phi_3(h J_n) + 60 phi_4(h J_n)) h r(U_2)
	//                                     + (-250/81 phi_3(h J_n) + 500/27 phi_4 (h J_n)) h r(U_3)
	
	N_VLinearSum(1.0, y, 1.0, scratchVec3, y); 
	// y is now new y, i.e., u_{n+1} = y
	
	t = t + h;
	//printf("Current time t = %g\n", t); // this option to be uncommented it you want to see present time
        integratorStats->Step();
    }
    
    return integratorStats;
}
