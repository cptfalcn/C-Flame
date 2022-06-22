//The return code version of KIOPS header
#include <stdio.h>
#include <math.h>
#include "Epic.h"
#include <nvector/nvector_serial.h>
//#include <cvode/cvode.h>
#include <sundials/sundials_nvector.h>
// #include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix  */
// #include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver  */
// #include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
// #include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */
// #include <sunlinsol/sunlinsol_spgmr.h>  /* SPGMR SUNLinearSolver                */
// #include <sunlinsol/sunlinsol_spbcgs.h>       /* access to SPBCGS SUNLinearSolver            */
// #include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver         */
#include <iostream>
#include <fstream>
#include <iomanip>
//#include "TChem_Impl_IgnitionZeroD_Problem.hpp" // here is where Ignition Zero D problem is implemented
#include <memory>
#define BAR "===================="
class KiopsRetCode{
    public:
        KiopsRetCode(int maxNumVectors, int m_max, N_Vector templateVector, int vecLength);
        ~KiopsRetCode();
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
