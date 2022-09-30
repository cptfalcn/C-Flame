//Starts the return code version of KIOPS
#include "KiopsRetCode.h"
int fac1(int n) { return (n == 1 || n == 0) ? 1 : fac1(n - 1) * n; }


KiopsRetCode::KiopsRetCode(int maxNumVectors, int m_max, N_Vector templateVector, int vecLength) 
    : MaxNumInputVectors(maxNumVectors),
      MaxPhiOrder(maxNumVectors-1),
      VecLength(vecLength),
      M_max(m_max),
      M_min(min(vecLength, 10)),
      Orth_len(vecLength),
      MatrixSize(m_max + 1),
      PhiMatrixSize(MatrixSize + 1)
{
    V = N_VCloneVectorArray(MatrixSize, templateVector);
    if ( V == NULL ) {
        throw std::bad_alloc();
    }

    V_aug = new realtype[MatrixSize*MaxPhiOrder];
    if ( V_aug == NULL ) {
        throw std::bad_alloc();
    }

    scratchVector = N_VClone(templateVector);
    if ( scratchVector == NULL ) {
        throw std::bad_alloc();
    }

    B = N_VCloneVectorArray(MaxPhiOrder, templateVector);
    if ( B == NULL ) {
        throw std::bad_alloc();
    }

    H = new realtype[MatrixSize*MatrixSize];
    if ( H == NULL ) {
        throw std::bad_alloc();
    }
    for (int i = 0; i < MatrixSize*MatrixSize; i++) {
        H[i] = 0.0;
    }

    phiMatrix = new realtype[PhiMatrixSize*PhiMatrixSize];
    if ( phiMatrix == NULL ) {
        throw std::bad_alloc();
    }
    
    phiMatrixSkipped = new realtype[MatrixSize*MatrixSize];
    if ( phiMatrixSkipped == NULL ) {
        throw std::bad_alloc();
    }
    
    w_aug = new realtype[MaxPhiOrder];
    if ( w_aug == NULL ) {
        throw std::bad_alloc();
    }

    expm = new ExpPade(PhiMatrixSize);
}

KiopsRetCode::~KiopsRetCode() {
    N_VDestroyVectorArray(V, M_max);
    N_VDestroyVectorArray(B, MaxPhiOrder);
    N_VDestroy(scratchVector);
    delete[] H;
    delete[] phiMatrix;
    delete[] phiMatrixSkipped;
    delete[] w_aug;
    delete expm;
}


//==================
//Kiops compute phi
//==================
int KiopsRetCode::ComputeKry(const int numVectors, N_Vector* inputVectors, const realtype timePoints[], 
	const int numTimePoints, N_Vector* outputVectors, JTimesV* jtimesv, realtype h, 
	realtype tol, int& m, KrylovStats* krylovStats) 
{
	int p = numVectors - 1;
	//Set optional dump filestreams
	// ofstream Hfile;
	// ofstream Phifile;
	// //Open said filestreams
	// Hfile.open("FailedIsoHMat.txt", std::ios_base::app);
	// Phifile.open("FailedIsoPhiMat.txt", std::ios_base::app);
    // We only allow m to vary between mmin and mmax
    m = max(M_min, min(M_max, m));

    // Initial condition
    N_VScale(1.0, inputVectors[0], outputVectors[0]);
    
    realtype normU = 0.0;
    for (int i = 1; i < numVectors; i++) {
        realtype sum = N_VL1Norm(inputVectors[i]);
        normU = max(normU, sum);
    }

    realtype nu, mu;
    if (normU != 0.0) {
        double ex = ceil(log2(normU));
        nu = pow(2, -ex);
        mu = pow(2, ex);
    } else {
        nu = 1;
        mu = 1;
    }

    for (int i = 0; i < p; i++) {
        N_VScale(nu, inputVectors[p-i], B[i]);
    }

    realtype sign = timePoints[numTimePoints-1] > 0 ? 1 : -1;
    realtype t_now = 0.0;
    realtype t_out = abs(timePoints[numTimePoints-1]);
    realtype tau = t_out;
    bool happy = false;
    int totalMatrixExponentials = 0;

    // Setting the safety factors and tolerance requirements
    realtype gamma, gamma_mmax;
    if (t_out > 1.0) {
        gamma = 0.2;
        gamma_mmax = 0.1;
    } else {
        gamma = 0.9;
        gamma_mmax = 0.6;
    }

    realtype delta = 1.4;

    // Used in the adaptive selection
    int oldm = -1;
    double oldtau = -1;
    double omega = -1;
    bool orderold = true;
    bool kestold = true;

    int l = 0;
    int j = 0;
    double beta;
    int ireject = 0;
    int reject = 0;

    while (t_now < t_out) {
        if (j == 0) {
            for (int k = 1; k < p; k++) {
                int i = p - k;
                w_aug[k - 1] = pow(t_now, i) / fac1(i) * mu;
            }
            w_aug[p - 1] = mu;

            double sum = N_VDotProd(outputVectors[l],outputVectors[l]);
            for (int i = 0; i < p; i++) {
                sum += w_aug[i] * w_aug[i];
            }
            beta = sqrt(sum);

            N_VScale(1.0/beta, outputVectors[l], V[0]);
            for (int i = 0; i < p; i++) {
                V_aug[i] = w_aug[i] / beta;
            }
        }

        // Incomplete orthogonalization process
        while (j < m) 
		{
            j++;

            // Augmented matrix - vector product
            jtimesv->ComputeJv(V[j-1], V[j]);
            N_VLinearCombination(p, V_aug + (j-1)*p, B, scratchVector);
            N_VLinearSum(h, V[j], 1.0, scratchVector, V[j]);

            for (int k = 0; k < p - 1; k++) {
                V_aug[j*p + k] = V_aug[(j-1)*p + k + 1];
            }
            V_aug[j*p + p - 1] = 0.0;

            // Modified Gram-Schmidt
            for (int i = max(0, j - Orth_len); i < j; i++) {
                int idx = i + (j - 1) * MatrixSize;
                H[idx] = N_VDotProd(V[i], V[j]);
                for (int k = 0; k < p; k++) {
                    H[idx] += V_aug[i*p + k] * V_aug[j*p + k];
                }

                N_VLinearSum(1.0, V[j], -H[idx], V[i], V[j]);
                for (int k = 0; k < p; k++) {
                    V_aug[j*p + k] -= H[idx] * V_aug[i*p + k];
                }
            }

            double nrm = N_VDotProd(V[j],V[j]);
            for (int k = 0; k < p; k++) {
                nrm += V_aug[j*p + k] * V_aug[j*p + k];
            }
            nrm = sqrt(nrm);
			if(isnan(nrm))
            {
				//cout << "norm is NaN\n";
                return 1;
            }
            // Happy breakdown
            if (nrm < tol) {
                happy = true;
                break;
            }

            H[j + (j - 1) * MatrixSize] = nrm;
            N_VScale(1.0 / nrm, V[j], V[j]);
            for (int k = 0; k < p; k++) {
                V_aug[j*p + k] /= nrm;
            }
        }


        // setup PhiMatrix
        for (int i = 0; i < j; i++) {
           for (int k = 0; k < j; k++) {
                phiMatrix[k + i * (j + 1)] = tau * sign * H[k + i * MatrixSize];
            }
            phiMatrix[j + i * (j + 1)] = 0.0;
        }
        
        phiMatrix[j * (j + 1)] = tau * sign;
        for (int k = 1; k < j + 1; k++) {
            phiMatrix[k + j * (j + 1)] = 0.0;
        }
        
        
        // Compute the exponential of the augmented matrix
        expm->Compute(j + 1, phiMatrix);
        totalMatrixExponentials++;
        
        double m_new, tau_new;
        if (happy) {
            // Happy breakdown; wrap up
            omega = 0;
            happy = false;
            m_new = m;
            tau_new = min(t_out - (t_now + tau), tau);
        } else 
		{
            // Local truncation error estimation
            double err = abs(beta * H[j + (j - 1) * MatrixSize] * phiMatrix[(j - 1) + j * (j + 1)]);
			//======================
			//Additional error state check
			//======================
			if(isnan(err))//New, put a better condition in here
			{
				// cout << "Dumping the phi matrix to PhiMat.txt\n";
				// cout << "Dumping the H matrix to HMat.txt\n"; 
				// std :: cout << "Error State has occurred: Err NaN\n";
				// std :: cout << "Beta: " << beta << "\n";
				// std :: cout << "Phi Mat value: " <<phiMatrix[j-1 + j * (j+1)] << "\n";
				// std :: cout << "H Mat value: " << H[j + (j-1) * MatrixSize] << "\n";
				// cout << expm->Compute(j + 1, phiMatrix) << endl;
				// cout << j << endl;
				return 1;
			}
			//===============
			//End error state check
			//===============
            // Error for this step
			double oldomega = omega;
			omega = t_out * err / (tau * tol);
            // Estimate order
			double order;
			if (m == oldm && tau != oldtau && ireject >= 1) {
				order = max(1.0, log(omega / oldomega) / log(tau / oldtau));
				orderold = false;
			} else if (orderold || ireject == 0) {
                orderold = true;
                order = j / 4;
            } else {
                orderold = true;
            }

            // Estimate k
            double kest;
            if (m != oldm && tau == oldtau && ireject >= 1) {
                kest = max(1.1, pow(omega / oldomega, 1 / (oldm - m)));
                kestold = false;
            } else if (kestold || ireject == 0) {
                kestold = true;
                kest = 2;
            } else {
                kestold = true;
            }

            //double remaining_time = omega > delta ? t_out - t_now : t_out - (t_now + tau);
            
            double remaining_time = t_out - t_now;
            if( omega - delta < 0)
            {
                remaining_time -= tau;
            }

            // Krylov adaptivity
            double same_tau = min(remaining_time, tau);
            double tau_opt = tau * pow(gamma / omega, 1.0 / order);

			// if(isnan(tau_opt)||tau_opt<1e-10)		//Added by me to try to stop endless cycling conditions.
			// {
			// 	tau_opt= tau/100;	//Same as above.
			// }
            tau_opt = min(remaining_time, max(tau / 5.0, min(5.0 * tau, tau_opt)));
	    	int m_opt = ceil(j + log(omega / gamma) / log(kest));
            m_opt = max(M_min, min(M_max, max((int)floor(3.0 / 4.0 * m),
                            min(m_opt, (int)ceil(4.0 / 3.0 * m)))));

            if (j == M_max) {
                if (omega > delta) {
                    m_new = j;
                    tau_new = tau * pow(gamma_mmax / omega, 1 / order);
                    tau_new = min(t_out - t_now, max(tau / 5, tau_new));
                } else {
                    tau_new = tau_opt;
                    m_new = m;
                }
            } else {
                m_new = m_opt;
                tau_new = same_tau;
            }


        }

        // Check error against target
        if (omega <= delta) {
            // Yep, got the required tolerance; update
            krylovStats->NewProjection(j, ireject);
            reject += ireject;

            // Udate for t in the interval (t_now, t_now + tau)
            int blownTs = 0;
            double nextT = t_now + tau;
            for (int k = l; k < numTimePoints; k++) {
                if (abs(timePoints[k]) < abs(nextT)) {
                    blownTs++;
                }
            }

            if (blownTs != 0) {
                // Copy current w to w we continue with.
                N_VScale(1.0, outputVectors[l], outputVectors[l+ blownTs]);

                for (int k = 0; k < blownTs; k++) {
                    double tauPhantom = timePoints[l + k] - t_now;

                    // setup PhiMatrixSkipped
                    for (int i = 0; i < j; i++) {
                        for (int k = 0; k < j; k++) {
                            phiMatrixSkipped[k + i*j] = tauPhantom * sign * H[k + i * MatrixSize];
                        }
                    }
                    
                    // Compute the exponential of the augmented matrix
                    expm->Compute(j, phiMatrixSkipped);
                    totalMatrixExponentials++;

                    N_VLinearCombination(j, phiMatrixSkipped, V, outputVectors[l+k]);
                    N_VScale(beta, outputVectors[l+k], outputVectors[l+k]);
                }

                // Advance l.
                l = l + blownTs;
            }

            // Using the standard scheme
            N_VLinearCombination(j, phiMatrix, V, outputVectors[l]);
            N_VScale(beta, outputVectors[l], outputVectors[l]);
            
            // Update t
            t_now += tau;

            j = 0;
            ireject = 0;

        } else {
            // Nope, try again
            ireject = ireject + 1;
        }

        oldtau = tau;
        tau = tau_new;
	oldm = m;
	m = m_new;
		if (ireject > 2*VecLength)
		{
			// std::cout<<"==========KIOPS has stalled, rejects========="<<std::endl;
			return 1;
		}
		if (tau_new==0&&t_now!=t_out)
		{
			// std::cout << "\n============KIOPS stalled, tau==========" << std::endl;
			return 1;
		}
		/**/
    }
    krylovStats->numMatrixExponentials = totalMatrixExponentials;

	return 0;
}
