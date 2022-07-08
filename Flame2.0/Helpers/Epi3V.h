#include <stdio.h>
#include <math.h>
#include "Epic.h"
#include <nvector/nvector_serial.h>
#include <cvode/cvode.h>
#include <sundials/sundials_nvector.h>
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix  */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver  */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */
#include <sunlinsol/sunlinsol_spgmr.h>  /* SPGMR SUNLinearSolver                */
#include <sunlinsol/sunlinsol_spbcgs.h>       /* access to SPBCGS SUNLinearSolver            */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver         */
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TChem_Impl_IgnitionZeroD_Problem.hpp" // here is where Ignition Zero D problem is implemented

#include <memory>
#define BAR "===================="
class NewKiops{
    public:
        NewKiops(int maxNumVectors, int m_max, N_Vector templateVector, int vecLength);
        ~NewKiops();
        int ComputeKry(const int numVectors, N_Vector* inputVectors, const realtype timePoints[], const int numTimePoints, N_Vector* outputVectors, JTimesV* jtimesv, realtype h, realtype tol, int& m, KrylovStats* krylovStats);
    private:
        const int MaxNumInputVectors;
        const int MaxPhiOrder;
        const int VecLength;  

        const int M_max;
        const int M_min;

        const int Orth_len;

        const int MatrixSize;
        const int PhiMatrixSize;

        N_Vector* V;
        realtype* V_aug;

        realtype* H;
        realtype* phiMatrix;
        realtype* phiMatrixSkipped;

        N_Vector w;
        realtype* w_aug; 

        N_Vector scratchVector;
        N_Vector* B;

       ExpPade* expm; 
};



class Epi3VChem_KIOPS
{
	public:
		//Create without jtv or jacf
		Epi3VChem_KIOPS(CVRhsFn f, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
		//Create with f, jtv, ToDo: finite difference jac and jtv
		Epi3VChem_KIOPS(CVRhsFn f, CVSpilsJacTimesVecFn jtv, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
		//Create with f, jtv, jacf.
		Epi3VChem_KIOPS(CVRhsFn f, CVSpilsJacTimesVecFn jtv, CVLsJacFn jacf, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);
		~Epi3VChem_KIOPS();
		//IntegratorStats *Integrate(realtype h, realtype t0, realtype tFinal, const int numBands, N_Vector y, realtype krylovTol, int startingBasisSizes[]);
 		IntegratorStats* Integrate(const realtype hStart, const realtype hMax, const realtype absTol,
			const realtype relTol, const realtype t0, const realtype tFinal,
			const int numBands, int basisSizes[], N_Vector y);

		IntegratorStats* NewIntegrate(const realtype hStart, const realtype hMax, const realtype absTol,
			const realtype relTol, const realtype t0, const realtype tFinal,
			int basisSizes[], N_Vector y);
		
	private:
    	CVRhsFn f;
    	CVSpilsJacTimesVecFn jtv;
		CVLsJacFn jacf;
    	EPICNumJacDelta delta;
    	void *userData;
    	Kiops *krylov;
		NewKiops *NewKrylov;
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
    	Epi3VChem_KIOPS(const Epi3VChem_KIOPS &);  // no implementation
    	Epi3VChem_KIOPS & operator=(const Epi3VChem_KIOPS &);  // no implementation
        void Clean(N_Vector y, int Length);
		void CheckNanBlowUp(N_Vector, int);
		// void MatMult(N_Vector R, N_Vector A, N_Vector B, int row, int col);
		// void Phi2(N_Vector Jac, realtype h, N_Vector Result);
		int TrackProgress(realtype FinalTime, realtype TNext, realtype PercentDone, int ProgressDots);
		
};


