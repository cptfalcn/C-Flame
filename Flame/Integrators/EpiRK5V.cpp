/*
 * ===========================================================================================
 * 
 * This file contains implementation of the Epirk5P1 method using old Epic impelmentartion, 
 * i.e., the name EpiRK5V. 
 * EpiRK5V method is method with variable time step and uses adaptive Krylov
 * method to compute Krylov iterations. 
 * 5th order Runge-Kutta exponential propagation integrator.
 * EpiRK5V method is part of Epic Package. 
 * All constant declarations are declared in header file.
 * 
 * Higer order method (Epirk5P1)
 * U_1     = u_n + a_11 phi_1(g_11 h J_n) h f(u_n)
 * U_2     = u_n + a_21 phi_1(g_21 h J_n) h f(u_n) + a_22 phi_1(g_22 h J_n) h r(U_1)
 * u_{n+1} = u_n + b_1 phi_1(g_31 h J_n) h f(u_n) + b_2 phi_1(g_32 h J_n) h r(U_1) 
 *               + b_3 phi_3(g_33 h J_n) h (-2 r(U_1) + r(U_2))
 * 
 * Lower order method (4th order method)
 * U_1     = u_n + a_11 phi_1(g_11 h J_n) h f(u_n)
 * U_2     = u_n + a_21 phi_1(g_21 h J_n) h f(u_n) + a_22 phi_1(g_22 h J_n) h r(U_1)
 * u_{n+1} = u_n + b_1 phi_1(g_31 h J_n) h f(u_n) + b_2 phi_1(g_32 h J_n) h r(U_1) 
 *               + b_3 phi_3(g_33 h J_n) h (-2 r(U_1) + r(U_2))
 * 
 * Since g_22 for higer order method is arbitrary, in this code it is going to be 1/2. Once 
 * phi_1(g_22 h J_n) h r(U_1) for 5th order method is computed, it will be used for 4th order
 * method in place of b_2 phi_1(g_32 h J_n) h r(U_1).
 * 
 * Exact constants you can find in header file. 
 * 
 * -------------------------------------------------------------------------------------------
 * 
 * Last revision on 02/08/2016 by Ilija Jegdic
 * 
 * ===========================================================================================
 */

#include <sundials/sundials_math.h>
#include <math.h>
#include "EpiRK5V.h"
#include "EpicSundials.h"
#include "EpicConst.h"
#include <stdlib.h>

using namespace EpiRK5VNamespace;

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

EpiRK5V::EpiRK5V(CVRhsFn f, CVSpilsJacTimesVecFn jtv, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    this->f = f;
    this->jtv = jtv;
    this->delta = nullptr;
    this->userData = userData;

    krylov = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder1 + 1, maxKrylovIters, tmpl, NEQ);
    krylov2 = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder1 + 1, maxKrylovIters, tmpl, NEQ);
    krylov3 = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder2 + 1, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(NumProjectionsPerStep);
    
    fy = N_VClone(tmpl);
    hfy = N_VClone(tmpl);
    hb1 = N_VClone(tmpl);
    hb2 = N_VClone(tmpl);
    r1 = N_VClone(tmpl);
    r2 = N_VClone(tmpl);
    r22 = N_VClone(tmpl);
    r3 = N_VClone(tmpl);
    r32 = N_VClone(tmpl);
    r33 = N_VClone(tmpl);
    r3LowOrder = N_VClone(tmpl);
    r32LowOrder = N_VClone(tmpl);
    r33LowOrder = N_VClone(tmpl);
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

EpiRK5V::EpiRK5V(CVRhsFn f, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    this->f = f;
    this->jtv = nullptr;
    this->delta = nullptr;
    this->userData = userData;

    krylov = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder1 + 1, maxKrylovIters, tmpl, NEQ);
    krylov2 = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder1 + 1, maxKrylovIters, tmpl, NEQ);
    krylov3 = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder2 + 1, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(NumProjectionsPerStep);

    fy = N_VClone(tmpl);
    hfy = N_VClone(tmpl);
    hb1 = N_VClone(tmpl);
    hb2 = N_VClone(tmpl);
    r1 = N_VClone(tmpl);
    r2 = N_VClone(tmpl);
    r22 = N_VClone(tmpl);
    r3 = N_VClone(tmpl);
    r32 = N_VClone(tmpl);
    r33 = N_VClone(tmpl);
    r3LowOrder = N_VClone(tmpl);
    r32LowOrder = N_VClone(tmpl);
    r33LowOrder = N_VClone(tmpl);
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

EpiRK5V::EpiRK5V(CVRhsFn f, EPICNumJacDelta delta, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    this->f = f;
    this->jtv = nullptr;
    this->delta = delta;
    this->userData = userData;

    krylov = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder1 + 1, maxKrylovIters, tmpl, NEQ);
    krylov2 = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder1 + 1, maxKrylovIters, tmpl, NEQ);
    krylov3 = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder2 + 1, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(NumProjectionsPerStep);

    fy = N_VClone(tmpl);
    hfy = N_VClone(tmpl);
    hb1 = N_VClone(tmpl);
    hb2 = N_VClone(tmpl);
    r1 = N_VClone(tmpl);
    r2 = N_VClone(tmpl);
    r22 = N_VClone(tmpl);
    r3 = N_VClone(tmpl);
    r32 = N_VClone(tmpl);
    r33 = N_VClone(tmpl);
    r3LowOrder = N_VClone(tmpl);
    r32LowOrder = N_VClone(tmpl);
    r33LowOrder = N_VClone(tmpl);
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
}

EpiRK5V::~EpiRK5V()
{
    delete krylov;
    delete krylov2;
    delete krylov3;
    delete integratorStats;

    N_VDestroy(fy);
    N_VDestroy(hfy);
    N_VDestroy(hb1);
    N_VDestroy(hb2);
    N_VDestroy(r1);
    N_VDestroy(r2);
    N_VDestroy(r22);
    N_VDestroy(r3);
    N_VDestroy(r32);
    N_VDestroy(r33);
    N_VDestroy(r3LowOrder);
    N_VDestroy(r32LowOrder);
    N_VDestroy(r33LowOrder);
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
 * hStart         starting time step size
 * hMax           maximum time step size
 * absTol         absolute tolerance
 * relTol         relative tolerance
 * t0             starting time
 * tFinal         ending time
 * numBands
 * basisSizes[]   maximum basis size allowed for the Arnoldi iteration
 * y              initial condition
 * 
 * Output:
 * Statistics
 * 
 * ===========================================================================================
 */

IntegratorStats *EpiRK5V::Integrate(const realtype hStart, const realtype hMax, const realtype absTol, const realtype relTol, const realtype t0, const realtype tFinal, const int numBands, int basisSizes[], N_Vector y)
{
    if( hStart < ZERO )
    {
        printf("Initial time step h is to small. \n");
        exit(EXIT_FAILURE);
    }
    if( tFinal < t0 )
    {
        printf("Starting time is larger the end time. \n");
        exit(EXIT_FAILURE);
    }
    realtype t = t0;
    realtype hNew = hStart;
    realtype h = hStart;
    if (hMax < hStart)
    {
        hNew = hMax;
        h = hMax;
    }

    N_Vector yTemp = N_VClone(y);
    N_VConst(1.0, yTemp);
    realtype sqrtN = EPICRSqrt(N_VDotProd(yTemp, yTemp));
    realtype krylovTol = 0.1 * sqrtN * absTol;
//    realtype krylovTol = 1.0e-14;

    // Main integration loop.
    
    bool finalStep = false;
    while (t < tFinal)
    {
        realtype err = 5.0;
        while (err > 1)
        {
            h = hNew;
	    
	    // f is RHS function from the problem; y = u_n
	    f(t, y, fy, userData);
            N_VScale(h, fy, hfy);
	    JTimesV jtimesv(jtv, f, delta, t, y, fy, userData, tmpVec);
	    
	    // Stage 1
	    N_Vector stage1InputVecs[] = {zeroVec, hfy};
	    N_Vector stage1OutputVecs[] = {r1, r2, r3}; 	    
	    krylov->Compute(numBands, Stage1NumInputVecs, stage1InputVecs, g1Times, g1NumTimes, stage1OutputVecs, &jtimesv, h, krylovTol, basisSizes[0], &integratorStats->krylovStats[0]);
	    // computes phi_k(g_{i j} h J)*hfy and stores it in ri, ri is r(U_i)
	    // in this case r1 = phi_1( g_11 h J_n ) * hfy 
	    // and          r2 = phi_1( g_21 h J_n ) * hfy
	    // and          r3 = phi_1( g_31 h J_n ) * hfy
	    
	    // we will need r3 for lower order as well, so we store it in r3Loworder
	    N_VScale(1.0, r3, r3LowOrder);
	    
	    // State U_1 is being computed and result is stored in r1
	    N_VLinearSum(1.0, y, a11, r1, r1); // r1 = y + a11 * r1 = y + a11 * phi_1( ) * hfy 
	    
	    // Stage 3 - High-order part
	    // In this part we compute residual of U_1, i.e., h*r(U_1). Result will be stored in hb1 variable
	    N_VLinearSum(1.0, r1, -1.0, y, scratchVec1); // scratchVec1 = a11 * phi( ) * hfy
	    jtimesv.ComputeJv(scratchVec1, scratchVec2); // scratchVec1 can be reused now, scratchVec2 is Jacobian of scratchVec1
	    f(t, r1, hb1, userData); // f(t, r1) = hb1
	    N_VLinearSum(h, hb1, -1.0, hfy, hb1); // hb1 = h * hb1 - hfy
	    N_VLinearSum(1.0, hb1, -1.0*h, scratchVec2, hb1);  // scratchVec2 can be reused now, hb1 contains residual of U_1
	    // In order to compute residual we did the following
	    // hb1 = hb1 - h*scratchVec2
	    //     = h*hb1 - hfy - Jacobian(scratchVec1)
	    //     = f(t, r1) - hfy - Jacobian(r1-y)
	    //     = f(t, y + a11*phi( )*hfy) - hfy - Jacobian( a11*phi( )*hfy )*hfy
	    //     = h*r(U_1) <- residual
	    
	    // scratchVec1 and scratchVec2 are free variables

	    // Stage 2
	    N_Vector stage2InputVecs[] = {zeroVec, hb1};
	    N_Vector stage2OutputVecs[] = {r22, r32}; 
	    // r22 = phi_1(g22 h J_n) hb1;    r32 = phi_1(g32 h J_n) hb1; 
	    // r22 = phi_1(g22 h J_n) r(U_1); r32 = phi_1(g32 h J_n) r(U_1); 
	    
	    krylov2->Compute(numBands, Stage2NumInputVecs, stage2InputVecs, g2Times, g2NumTimes, stage2OutputVecs, &jtimesv, h, krylovTol, basisSizes[0], &integratorStats->krylovStats[1]);
	    // computes phi(g_{i j} h J)*hfy and stores it in r22
	    // in this case r22 = phi_1( g22 h J ) * r(U_1) = phi_1( g32LowOrder h J ) * r(U_1)
	    // in this case r32 = phi_1( g32 h J ) * r(U_1) 
	    
	    // We will need r22 in lower order so we store it in r32LowOrder
	    N_VScale(1.0, r22, r32LowOrder);
	    
	    // state U_2
	    // first part
	    N_VLinearSum(1.0, y, a21, r2, r2); // r2 = y + a21 * r2 = y + a21 * phi( ) * hfy
	    
	    // second part 
	    N_VLinearSum(1.0, r2, a22, r22, r2); // r2 is now whole U_2
	    // r2 = y + a21 * phi( ) * hfy + a22 * phi( ) * r(U_1)
	    
	    // In this part we compute residual of U_2, i.e., h*r(U_2). Result will be stored in hb1 variable
	    N_VLinearSum(1.0, r2, -1.0, y, scratchVec1); // scratchVec1 = r2 - y = a21 * phi( ) * hfy + a22 * phi( ) * r(U_1)
	    jtimesv.ComputeJv(scratchVec1, scratchVec2); // can now re-use scratchVec1
	    // scratchVec2 is Jacobian * scratchVec1
	    f(t, r2, hb2, userData); // f(t, r2) = hb2
	    N_VLinearSum(h, hb2, -1.0, hfy, hb2); // hb2 = h * hb2 - hfy
	    N_VLinearSum(1.0, hb2, -1.0*h, scratchVec2, hb2);
	    // Similar like for hb1, hb2 is now h*r(U_2);
	    
	    // scratchVec1 and scratchVec2 are free variables
	    
	    // State u_{n+1} is now being computed
	    N_VLinearSum(1.0, y, b1, r3, r3); // r3 = y + b1 * r3 = y + b1 * phi_1( ) * hfy
	    N_VLinearSum(1.0, r3, b2, r32, r3); // r3 = y + b1 * phi_1( ) hfy + b2 * phi_1( ) r(U_1)
	    
	    // Now to compute b_3 phi_3 (g_33 h J_n) h (-2 r(U_1) + r(U_2))
	    // Taken variables hb1 = h*r(U_1); hb2 = h*r(U_2); r3 and y
	    // Free variables scratchVec1 and scratchVec2
	    
	    N_VLinearSum(-2.0, hb1, 1.0, hb2, scratchVec1); // scratchVec1 = h (-2 r(Y_1) + r(Y_2))
	    
	    N_Vector stage3InputVecs[] = {zeroVec, zeroVec, zeroVec, scratchVec1}; // we can now reuse scratchVec1
	    N_Vector stage3OutputVecs[] = {r33, r33LowOrder}; 
	    
	    krylov3->Compute(numBands, Stage3NumInputVecs, stage3InputVecs, g3Times, g3NumTimes, stage3OutputVecs, &jtimesv, h, krylovTol, basisSizes[0], &integratorStats->krylovStats[2]);
	    // r33 = phi_3 (g_33 h J_n) h (-2 r(Y_1) + r(Y_2))
	    
	    N_VLinearSum(1.0, r3, b3, r33, scratchVec1); // scratchVec1 contains tentative y
	    // scratchVec1 = y + b1 * phi_1( ) hfy + b2 * phi_1( ) r(U_1) + r33
	    //                 + b3 phi_3 (g_33 h J_n) h (-2 r(Y_1) + r(Y_2))
	    
	    // Lower order part
	    // State u_{n+1} for low orderis being computed and results is stored in scratchVec2
	    N_VLinearSum( 1.0, y, b1, r3LowOrder, scratchVec2 );
	    N_VLinearSum( b2, r32LowOrder, 1.0,  scratchVec2,  scratchVec2 );
	    N_VLinearSum( b3, r33LowOrder, 1.0,  scratchVec2,  scratchVec2 );
	    // scratchVec2 = y + b1 r3LowOrder + b2 r32LowOrder + b3 r33LowOrder
	    
	    // At this moment scratchVec1 contains tentative y and scratchVec2 contains lower order y
	    
            // Estimate error.
            N_VLinearSum(1.0, scratchVec1, -1.0, scratchVec2, scratchVec2);
            N_VAbs(y, scratchVec3);
            N_VScale(relTol, scratchVec3, scratchVec3);
            N_VAddConst(scratchVec3, absTol, scratchVec3);
            N_VDiv(scratchVec2, scratchVec3, scratchVec2);  // can now re-use scratchVec3
            realtype norm = N_VDotProd(scratchVec2, scratchVec2);  // can now re-use scratchVec6
            norm = norm / NEQ;
            err = EPICRSqrt(norm);
            hNew = h * Fac * pow(err, Order);
            if (hNew > hMax)
            {
                hNew = hMax;
            }
            if (hNew < ZERO)
            {
                printf("There is possible singularity in the solution\n");
                exit(EXIT_FAILURE);
            }
            
            //printf("err = %f, hNew = %f\n", err, hNew);
        }

        // Copy tentative y (scratchVec1) to y.
        N_VLinearSum(1.0, scratchVec1, 0.0, y, y);
	// y is now new y, i.e., u_{n+1} = y
	
        integratorStats->Step();

        // If not final step, advance to next time step.
        if (finalStep)
        {
            break;
        }
        t = t + h;
//	printf("Current time t = %g\n", t); // this option to be uncommented it you want to see present time
        if (t + hNew >= tFinal)
        {
            hNew = tFinal - t;
            finalStep = true;
        }
    }

    return integratorStats;
}
