/*
 * ===========================================================================================
 *
 * This file contains implementation of the second-order Exponential Polynomial Block Method. 
 *
 * -------------------------------------------------------------------------------------------
 * 
 * Written by Tommaso Buvoli and Jared Stewart on 4/06/2022
 * 
 * ===========================================================================================
 */


#include "EpiP2_KIOPS.h"
#include <math.h>
#include <vector>
#include <stdlib.h>
#include "EpicConst.h"
#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.h>
#include <iostream>

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

EpiP2_KIOPS::EpiP2_KIOPS(CVRhsFn f, CVSpilsJacTimesVecFn jtv, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength,  int q, int kappa, void (*userPostProcess)(N_Vector)) : NEQ(vecLength)
{
    this->constuctor_helper(q, kappa, f, userData, maxKrylovIters,  tmpl, vecLength);
    this->jtv = jtv;
    this->userPostProcess = userPostProcess;
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

EpiP2_KIOPS::EpiP2_KIOPS(CVRhsFn f, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength,  int q, int kappa, void (*userPostProcess)(N_Vector)) : NEQ(vecLength)
{
    this->constuctor_helper(q, kappa, f, userData, maxKrylovIters,  tmpl, vecLength);
    this->jtv = nullptr;
    this->userPostProcess = userPostProcess;
}

void EpiP2_KIOPS::constuctor_helper(int q, int kappa, CVRhsFn f, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength)
{
    this->q = q; // TODO: add q,kappa to parameters
    this->kappa = kappa;
    
    this->f = f;
    this->userData = userData;

    krylov = new Kiops(this->q, maxKrylovIters, tmpl, NEQ);
    integratorStats = new IntegratorStats(this->q + this->kappa);
    fy = N_VClone(tmpl);
    tmpVec = N_VClone(tmpl);
    tmpVecJTV = N_VClone(tmpl);
    zeroVec = N_VClone(tmpl);
    // zero-out to avoid potential NANs
    N_VConst(0, fy);
    N_VConst(0, zeroVec);

	this->initWeights(this->q);
    this->initZ(this->q);
    for(int i = 0; i < this->q; i ++)
	{
		Y.push_back(N_VClone(zeroVec));
        Y_tmp.push_back(N_VClone(zeroVec));
        R.push_back(N_VClone(zeroVec));
        B.push_back(N_VClone(zeroVec));
	}
}

EpiP2_KIOPS::~EpiP2_KIOPS()
{
    delete krylov;
    delete integratorStats;
    
    N_VDestroy(fy);
    N_VDestroy(tmpVec);
    N_VDestroy(zeroVec);

    for(int i = 0; i < this->q; i ++)
	{
		N_VDestroy(Y[i]);
        N_VDestroy(Y_tmp[i]);
        N_VDestroy(R[i]);
        N_VDestroy(B[i]);
	}

    Y.clear();
    Y_tmp.clear();
    R.clear();
    B.clear();
}

void EpiP2_KIOPS::initWeights(int q)
{
    // initialize empty weight matrix
    for(int i = 0; i < q; i ++)
    {
        weights.push_back(vector<realtype> (this->q));
    }

    //weights[i,j] - weights for the ith phi function
    if( q == 3 )
    {
        weights[0][0] = 0;
        weights[0][1] = 0;
        weights[0][2] = 0;

        weights[1][0] = 0;
        weights[1][1] = (1.0 + sqrt(3.0))/2.0;
        weights[1][2] = (1.0 - sqrt(3.0))/2.0;

        weights[2][0] = 0;
        weights[2][1] = -1.0 * sqrt(3.0)/2.0;
        weights[2][2] = sqrt(3.0)/2.0;
    }
    else
    {
        // TODO: add code for arbitrary q
        std::cout << "Q > 3 is not yet supported" << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void EpiP2_KIOPS::initZ(int q)
{
    for(int i = 0; i < this->q; i ++)
    {
        z.push_back(0);
    }

    if( q == 3 )
    {
        z[0] = 0;
        z[1] = 1 - 1/sqrt(3.0);
        z[2] = 1 + 1/sqrt(3.0);
    }
    else
    {
        // TODO: add code for arbitrary q
        std::cout << "q > 3 is not yet supported" << std::endl;
        std::exit(EXIT_FAILURE);
    }

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
 * y              initial condition
 * krylovTol      tolerance for Krylov Iteration
 * basisSizes[]   maximum basis size allowed for the Arnoldi iteration
 * 
 * Output:
 * Statistics
 * 
 * ===========================================================================================
 */

IntegratorStats *EpiP2_KIOPS::Integrate(realtype h, realtype t0, realtype tFinal, N_Vector y, realtype krylovTol, int basisSizes[])
{
    realtype alpha_base = 2.0;
    realtype r = h / alpha_base;

	// initial conditions
    this->fillYs(y);	
    for(int i = 0; i < this->q + this->kappa; i ++)
    {
        this->step(r, 0, 0.0, krylovTol, basisSizes, integratorStats);
        this->postProcess();       
    }


    //timestepping
    int Nt = round((tFinal - t0) / h); // TODO: User should be allowed to specify h that does not evenly divide
    int t  = t0;

    for(int i = 0; i < Nt; i ++)
    {
        this->step(r, t, alpha_base, krylovTol, basisSizes, integratorStats);
        this->postProcess();

        for(int j = 0; j < this->kappa; j ++)
        {
            this->step(r, t, 0.0, krylovTol, basisSizes, integratorStats);
        }

        t = t + r * alpha_base;
        integratorStats->Step();
    }

    N_VScale(1.0, this->Y[0], y);

    return integratorStats;
}

void EpiP2_KIOPS::step(realtype r, realtype t0, realtype alpha, 
        realtype krylovTol, int basisSizes[], IntegratorStats* integratorStats)
{
	JTimesV JTV = this->evalRs(r, t0);
    this->buildB(r);

    // create timepoints this.z + alpha
    vector<realtype> output_times(this->z);
    for(int i = 0; i < this->q; i ++)
    {
        output_times[i] += alpha;
    }
    
    krylov->Compute(this->q, &B[0], &output_times[0], this->q, &this->Y_tmp[0], &JTV, r, krylovTol,
			basisSizes[0], &integratorStats->krylovStats[0]);

    for(int i = this->q - 1; i > -1; i --) // run loop backwards to avoid overwriting Y[0] (which is needed by other outputs)
    {
        N_VLinearSum(1.0, this->Y[0], 1.0, this->Y_tmp[i], this->Y[i]); // y_i = y_0 + y_tmp[i]
    }

}

void EpiP2_KIOPS::buildB(realtype r)
{
    // Form KIOPS input paramete Br:
    // B[i] = \sum_{j=0}^{q-1} w[i,j] R[j]
    for(int i = 0; i < this->q; i ++)
    {
        N_VScale(this->weights[i][0], R[0], B[i]);  // B[i] = self->weights[i,j] * R[0]
    
        for(int j = 1; j < this->q; j++)
        {
            N_VLinearSum(1.0, B[i], this->weights[i][j], R[j], B[i] );  // B[i] = B[i]+ weights[ij]*R[j]
        }
    }
    
    N_VLinearSum(1.0, B[1] , r, this->fy, B[1]);                        // B[1] = B[1] + r * this->fy
}

void EpiP2_KIOPS::fillYs(N_Vector Y)
{
	for(int i = 0; i < this->q; i ++)
    {
       N_VScale(1.0, Y, this->Y[i]);
    }
}

// Evaluates F(y_n), forms function for Jacobian vector product and evaluate the remainder
// at y_j=y(t_0+r(Z_j)), j = 1,2,3. Linearization happens at y_1
JTimesV EpiP2_KIOPS::evalRs(realtype r, realtype t0)
{
	this-> f(t0, this->Y[0], this->fy, this->userData);
    JTimesV JTV(this->jtv, this->f, this->delta, t0, this->Y[0], this->fy, this->userData, this->tmpVecJTV);

	for(int i = 0; i < this->q; i ++)
    {
        evalR(r, t0+r*this->z[i], this->Y[i], this->R[i], &JTV, this->Y[0]);
    }

    return JTV;
}

// R(y) = F(y) - [F_n + J_n(Y - Y_n)]
void EpiP2_KIOPS::evalR(realtype r, realtype t, N_Vector Y, N_Vector R, JTimesV* JTV, N_Vector Yn)
{
    // REMARK: this function assumes that this->fy already contains F(y_n)
    N_VLinearSum(1.0, Y, -1.0, Yn, this->tmpVec); // tmpVec = Y - Yn
	JTV->ComputeJv(this->tmpVec, R);	          // R = J * tmpVec = J (Y - Y_n)
    N_VLinearSum(-1.0, R, -1.0, this->fy, R);     // R = -R - f_n =  - [F_n + J_n(Y - Y_n)]
    this->f(t,Y,this->tmpVec,this->userData);     // tmpVec = F(y)
    N_VLinearSum(1.0, R, 1.0, this->tmpVec, R);	  // R = F(y) + R = F(y) - [F_n + J_n(Y - Y_n)]
    N_VScale(r, R, R);					          // R = r*R
}

void EpiP2_KIOPS::postProcess()
{
    for(int i = 0; i < this->q; i ++)
    {
        this->userPostProcess(this->Y[i]);
    }
}

// void EpiP2_KIOPS::print(N_Vector v, string name)
// {
//     std :: cout << "--- NVector: " << name << std :: endl;
//     N_VPrint_Serial(v);
// }

// void EpiP2_KIOPS::printv(vector<N_Vector> v, string name)
// {
//     std :: cout << "--- Vector<N_Vector>: " << name << std :: endl;
//     for(int i = 0; i < v.size(); i ++)
//     {
//         std :: cout << "N_Vector " << i << std::endl;
//         N_VPrint_Serial(v[i]);
//     }
    
// }