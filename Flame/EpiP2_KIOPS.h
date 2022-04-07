/* 
 * ===========================================================================================
 * 
 * This file is header file for EpiRK2.cpp. 
 * All constant declarations are contained in this file.
 * 
 * -------------------------------------------------------------------------------------------
 * 
 * Last revision on 02/08/2016 by Ilija Jegdic
 * 
 * ===========================================================================================
 */

#include <cvode/cvode.h>
#include <cvode/cvode_spils.h>
#include "Kiops.h"
#include "IntegratorStats.h"
#include "EpicSundials.h"
#include "EpicConst.h"
#include "EpicTypes.h"
#include <vector>

class EpiP2_KIOPS
{
public:
    EpiP2_KIOPS(CVRhsFn f, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength, int q = 3, int kappa = 0, void (*func)(N_Vector) = [](N_Vector n) {});
    EpiP2_KIOPS(CVRhsFn f, CVSpilsJacTimesVecFn jtv, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength, int q = 3, int kappa = 0, void (*postProcess)(N_Vector) = [](N_Vector n) {});
    EpiP2_KIOPS(CVRhsFn f, EPICNumJacDelta delta, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength, int q = 3, int kappa = 0, void (*postProcess)(N_Vector) = [](N_Vector n) {});
    ~EpiP2_KIOPS();
    IntegratorStats *Integrate(realtype h, realtype t0, realtype tFinal, N_Vector y, realtype krylovTol, int startingBasisSizes[]);
    
private:
    CVRhsFn f;
    CVSpilsJacTimesVecFn jtv;
    void (*userPostProcess)(N_Vector);
    EPICNumJacDelta delta;
    void *userData;
    Kiops *krylov;
    const int NEQ;
    IntegratorStats *integratorStats;

    N_Vector fy;           // stores f(y_n)
    N_Vector zeroVec;      // zero-ed out vector, used by krylov
    N_Vector tmpVec;       // used as a temporary work vector by remainder function
    N_Vector tmpVecJTV;    // extra temp vector for jacobian times transpose function
	//New stuff
	int q;
    int kappa;
    vector<vector<realtype>> weights;

	std::vector<N_Vector> Y;
    std::vector<N_Vector> Y_tmp;
	std::vector<N_Vector> R;
    std::vector<N_Vector> B;
    std::vector<realtype> z;

    // Disallow copying by "poisoning" the copy constructor and assignment operator,
    // i.e. declare private but provide no implementation.
    EpiP2_KIOPS(const EpiP2_KIOPS &);  // no implementation
    EpiP2_KIOPS & operator=(const EpiP2_KIOPS &);  // no implementation
	void constuctor_helper(int q, int kappa, CVRhsFn f, void *userData, int maxKrylovIters, N_Vector tmpl, const int vecLength);

    void step(realtype r, realtype t, realtype alpha, 
        realtype krylovTol, int basisSizes[], IntegratorStats* integratorStats);
	void fillYs(N_Vector Y);
	void buildB(realtype r);
    JTimesV evalRs(realtype r, realtype t0);
	void evalR (realtype r, realtype t, N_Vector Y, N_Vector R, JTimesV * jtv, N_Vector Yn);
    void initWeights(int q);
    void initZ(int q);
    void postProcess();
    // void print(N_Vector v, string name);
    // void printv(vector<N_Vector> v, string name);
};
