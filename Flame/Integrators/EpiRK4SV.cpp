/*
 * ===========================================================================================
 * 
 * This file contains implementation of the EpiRK4SV method. 
 * EpiRK4SV method is method with variable time step and uses adaptive Krylov
 * method to compute Krylov iterations. 
 * 4th order Runge-Kutta exponential propagation integrator.
 * EpiRK4SV method is part of Epic Package. 
 * All constant declarations are declared in header file.
 * 
 * Higer order method (EpiRK4SV)
 * U_2     = u_n + 1/2 phi_1(1/2 h J_n) h f(u_n)
 * U_3     = u_n + 2/3 phi_1(2/3 h J_n) h f(u_n)
 * u_{n+1} = u_n + phi_1(h J_n) h f(u_n) + (32 phi_3(h J_n) - 144 phi_4(h J_n)) h r(U_2)
 *               + (-27/2 phi_3(h J_n) + 81 phi_4 (h J_n)) h r(U_3)
 * 
 * Lower order method (3rd order method)
 * U_2     = u_n + 3/4 phi_1( 3/4 h J_n) h f(u_n)
 * u_{n+1} = u_n + phi_1( h J_n) h f(u_n) 
 *               + (p_{223} phi_3 (h J_n) + (128/9 - 4 p_{223}) phi_4(h J_n)) h r(U_2)
 * 
 * where p_{223} = 1
 * 
 * -------------------------------------------------------------------------------------------
 * 
 * Last revision on 02/13/2016 by Ilija Jegdic
 * 
 * ===========================================================================================
 */

#include <sundials/sundials_math.h>
#include <math.h>
#include "EpiRK4SV.h"
#include "EpicSundials.h"
#include "EpicConst.h"
#include <stdlib.h>

using namespace EpiRK4SVNamespace;

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

EpiRK4SV::EpiRK4SV(CVRhsFn f, CVSpilsJacTimesVecFn jtv, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    this->f = f;
    this->jtv = jtv;
    this->userData = userData;

    krylov = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder1 + 1, maxKrylovIters, tmpl, NEQ);
    krylov2 = new AdaptiveKrylovPhi(MaxPhiOrder2 + 1, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(NumProjectionsPerStep);

    fy = N_VClone(tmpl);
    hfy = N_VClone(tmpl);
    hb1 = N_VClone(tmpl);
    hb2 = N_VClone(tmpl);
    r1 = N_VClone(tmpl);
    r2 = N_VClone(tmpl);
    r3LowOrder = N_VClone(tmpl);
    scratchVec1 = N_VClone(tmpl);
    scratchVec2 = N_VClone(tmpl);
    scratchVec3 = N_VClone(tmpl);
    scratchVec4 = N_VClone(tmpl);
    scratchVec5 = N_VClone(tmpl);
    scratchVec6 = N_VClone(tmpl);
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

EpiRK4SV::EpiRK4SV(CVRhsFn f, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    this->f = f;
    this->jtv = nullptr;
    this->delta = nullptr;
    this->userData = userData;

    krylov = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder1 + 1, maxKrylovIters, tmpl, NEQ);
    krylov2 = new AdaptiveKrylovPhi(MaxPhiOrder2 + 1, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(NumProjectionsPerStep);

    fy = N_VClone(tmpl);
    hfy = N_VClone(tmpl);
    hb1 = N_VClone(tmpl);
    hb2 = N_VClone(tmpl);
    r1 = N_VClone(tmpl);
    r2 = N_VClone(tmpl);
    r3LowOrder = N_VClone(tmpl);
    scratchVec1 = N_VClone(tmpl);
    scratchVec2 = N_VClone(tmpl);
    scratchVec3 = N_VClone(tmpl);
    scratchVec4 = N_VClone(tmpl);
    scratchVec5 = N_VClone(tmpl);
    scratchVec6 = N_VClone(tmpl);
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

EpiRK4SV::EpiRK4SV(CVRhsFn f, EPICNumJacDelta delta, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    this->f = f;
    this->jtv = nullptr;
    this->delta = delta;
    this->userData = userData;

    krylov = new AdaptiveKrylovPhiMultipleTimes(MaxPhiOrder1 + 1, maxKrylovIters, tmpl, NEQ);
    krylov2 = new AdaptiveKrylovPhi(MaxPhiOrder2 + 1, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(NumProjectionsPerStep);

    fy = N_VClone(tmpl);
    hfy = N_VClone(tmpl);
    hb1 = N_VClone(tmpl);
    hb2 = N_VClone(tmpl);
    r1 = N_VClone(tmpl);
    r2 = N_VClone(tmpl);
    r3LowOrder = N_VClone(tmpl);
    scratchVec1 = N_VClone(tmpl);
    scratchVec2 = N_VClone(tmpl);
    scratchVec3 = N_VClone(tmpl);
    scratchVec4 = N_VClone(tmpl);
    scratchVec5 = N_VClone(tmpl);
    scratchVec6 = N_VClone(tmpl);
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

EpiRK4SV::~EpiRK4SV()
{
    delete krylov;
    delete integratorStats;

    N_VDestroy(fy);
    N_VDestroy(hfy);
    N_VDestroy(hb1);
    N_VDestroy(hb2);
    N_VDestroy(r1);
    N_VDestroy(r2);
    N_VDestroy(r3LowOrder);
    N_VDestroy(scratchVec1);
    N_VDestroy(scratchVec2);
    N_VDestroy(scratchVec3);
    N_VDestroy(scratchVec4);
    N_VDestroy(scratchVec5);
    N_VDestroy(scratchVec6);
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

IntegratorStats *EpiRK4SV::Integrate(const realtype hStart, const realtype hMax, const realtype absTol, const realtype relTol, const realtype t0, const realtype tFinal, const int numBands, int basisSizes[], N_Vector y)
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

            // Stage 1.
            N_Vector stage1InputVecs[] = {zeroVec, hfy};
            N_Vector stage1OutputVecs[] = {r1, r2, r3LowOrder};
            krylov->Compute(numBands, Stage1NumInputVecs, stage1InputVecs, g1Times, g1NumTimes, stage1OutputVecs, &jtimesv, h, krylovTol, basisSizes[0], &integratorStats->krylovStats[0]);
	    // computes phi_k(g_{i j} h J)*hfy and stores it in ri, ri is r(U_i)
	    // in this case r1         = phi_1( 1/2 h J_n ) * hfy 
	    // and          r2         = phi_1( 2/3 h J_n ) * hfy
	    // and          r3LowOrder = phi_1( 3/4 h J_n ) * hfy
	    
	    // States U_2 and U_3 are now being computed and results are stored in r1 and r2, respectively
	    N_VLinearSum(1.0, y, a21, r1, r1); // r1 = y + a21 * r1 = y + a21 * phi( ) * hfy
            N_VLinearSum(1.0, y, a31, r2, r2); // r2 = y + a31 * r1 = y + a31 * phi( ) * hfy

            // Stage 3 - High-order part
	    // In this part we compute residual of U_2, i.e., h*r(U_2). Result will be stored in hb1 variable
            N_VLinearSum(1.0, r1, -1.0, y, scratchVec1); // scratchVec1 = a21 * phi( ) * hfy
            jtimesv.ComputeJv(scratchVec1, scratchVec2); // scratchVec1 can be reused now, scratchVec2 is Jacobian of scratchVec1
            f(t, r1, hb1, userData); // f(t, r1) = hb1
            N_VLinearSum(h, hb1, -1.0, hfy, hb1); // hb1 = h * hb1 - hfy
            N_VLinearSum(1.0, hb1, -1.0*h, scratchVec2, hb1); // scratchVec2 can be reused now, hb1 contains residual of U_2
	    // In order to compute residual we did the following
	    // hb1 = hb1 - h*scratchVec2
	    //     = h*hb1 - hfy - Jacobian(scratchVec1)
	    //     = f(t, r1) - hfy - Jacobian(r1-y)
	    //     = f(t, y + a21*phi( )*hfy) - hfy - Jacobian( a21*phi( )*hfy )*hfy
	    //     = h*r(U_2) <- residual
	    
	    // scratchVec1 and scratchVec2 are free variables
	    
	    // In this part we compute residual of U_3, i.e., h*r(U_3). Result will be stored in hb2 variable
	    // We preform same caluculation we did for residual of U_2
            N_VLinearSum(1.0, r2, -1.0, y, scratchVec1);
            jtimesv.ComputeJv(scratchVec1, scratchVec2);
            f(t, r2, hb2, userData);
            N_VLinearSum(h, hb2, -1.0, hfy, hb2);
            N_VLinearSum(1.0, hb2, -1.0*h, scratchVec2, hb2);
	    // Similar like for hb1, hb2 is now h*r(U_3)
	
	    // scratchVec1 and scratchVec2 are free variables
	
	    // State u_{n+1} is now being computed
            N_VLinearSum(b2a, hb1, b2b, hb2, scratchVec1);
            N_VLinearSum(b3a, hb1, b3b, hb2, scratchVec2);
	
	    // scratchVec1 and scratchVec2 contain linear combinations of residuals

            N_Vector stage3InputVecs[] = {zeroVec, hfy, zeroVec, scratchVec1, scratchVec2};
            krylov2->Compute(numBands, Stage3NumInputVecs, stage3InputVecs, scratchVec3, 1.0, &jtimesv, h, krylovTol, basisSizes[1], &integratorStats->krylovStats[1]);
	    // scratchVec3 now holds high-order portion, i.e., 
	    // scratchVec3 = phi_1(h J_n) h f(u_n) + (32 phi_3(h J_n) - 144 phi_4(h J_n)) h r(U_2)
	    //                                     + (-27/2 phi_3(h J_n) + 81 phi_4 (h J_n)) h r(U_3)

            N_VLinearSum(1.0, y, 1.0, scratchVec3, scratchVec1);  
	    // scratchVec1 now holds tentative new y

            // Lower order
	    // State U_2 for low order is being computed and results is stored in r3LowOrder
	    N_VLinearSum(1.0, y, a21LowOrder, r3LowOrder, r3LowOrder); // r3LoweOrder = y + a21 * phi ( )*hfy
	    
	    // In this part we compute residual of U_2 for low order, i.e., h*r(U_2). Result will be stored in hb1 variable
	    N_VLinearSum(1.0, r3LowOrder, -1.0, y, scratchVec5); // scratchVec5 = a21 * phi( ) * hfy
            jtimesv.ComputeJv(scratchVec5, scratchVec6);  
            f(t, r3LowOrder, hb1, userData); // f(t, r3LowOrder) = hb1
            N_VLinearSum(h, hb1, -1.0, hfy, hb1); // hb1 = h * hb1 - hfy
            N_VLinearSum(1.0, hb1, -1.0*h, scratchVec6, hb1);// scratchVec6 can be reused now, hb1 contains residual of U_2
	    // In order to compute residual we did the following
	    // hb1 = hb1 - h*scratchVec6
	    //     = h*hb1 - hfy - Jacobian(scratchVec5)
	    //     = f(t, r3LowOrder) - hfy - Jacobian(r3LowOrder-y)
	    //     = f(t, y + a21LowOrder*phi( )*hfy) - hfy - Jacobian( a21LowOrder*phi( )*hfy )*hfy
	    //     = h*r(U_2) <- residual
	    
	    // State u_{n+1} low order is now being computed
            N_VScale(b2aLowOrder, hb1, scratchVec2);
	    N_VScale(b2bLowOrder, hb1, scratchVec5);
	
	    // scratchVec2 and scratchVec5 contain linear combinations of residuals

            N_Vector stage3LowOrderInputVecs[] = {zeroVec, hfy, zeroVec, scratchVec2, scratchVec5};
            krylov2->Compute(numBands, Stage3NumInputVecsLowOrder, stage3LowOrderInputVecs, scratchVec4, 1.0, &jtimesv, h, krylovTol, basisSizes[2], &integratorStats->krylovStats[2]);  
	    // scratchVec4 now holds low order portion, i.e., 
	    // scratchVec4 = phi_1( h J_n) h f(u_n) 
	    //                + (p_{223} phi_3 (h J_n) + (128/9 - 4 p_{223}) phi_4(h J_n)) h r(U_2)

            // Estimate error.
            N_VLinearSum(1.0, scratchVec3, -1.0, scratchVec4, scratchVec4);
            N_VAbs(y, scratchVec3);
            N_VScale(relTol, scratchVec3, scratchVec3);
            N_VAddConst(scratchVec3, absTol, scratchVec3);
            N_VDiv(scratchVec4, scratchVec3, scratchVec4);  // can now re-use scratchVec3
            realtype norm = N_VDotProd(scratchVec4, scratchVec4);  // can now re-use scratchVec6
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
