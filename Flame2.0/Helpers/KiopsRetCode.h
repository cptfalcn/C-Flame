//The return code version of KIOPS header
#include <stdio.h>
#include <math.h>
#include "Epic.h"
#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>

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
