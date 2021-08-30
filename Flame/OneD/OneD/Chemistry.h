#include <stdio.h>
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include "TChem_Util.hpp"
#include "TChem_KineticModelData.hpp"
#include "TChem_Impl_IgnitionZeroD_Problem.hpp" // here is where Ignition Zero D problem is implemented

//#define TCHEMPB TChem::Impl::IgnitionZeroD_Problem<TChem::KineticModelConstData<Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> >>
#define TCHEMPB TChem::Impl::IgnitionZeroD_Problem      <TChem::KineticModelConstData   <Kokkos::Device <Kokkos::Serial, Kokkos::HostSpace>  >  >
#define BAR "===================="

using real_type = double;
using ordinal_type = int;

//=====================
//Prototypes & Classes
//=====================
//Add an additional class elements to pass to epic.
class myPb : public TCHEMPB{
        public:
        //members
        ordinal_type  num_equations;
        N_Vector Jac;
};

//Definitions

int CHEM_RHS_TCHEM(realtype , N_Vector , N_Vector, void *);
int CHEM_COMP_JAC(N_Vector u, void* pb);
int CHEM_COMP_JAC_CVODE(realtype, N_Vector, N_Vector, SUNMatrix, void *, N_Vector, N_Vector, N_Vector);
int CHEM_JTV(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
int RHS_KAPPA(realtype , N_Vector , N_Vector, void *);
int Jtv_KAPPA(N_Vector , N_Vector , realtype , N_Vector , N_Vector , void * , N_Vector);
void MatrixVectorProduct(int, realtype *, N_Vector, N_Vector, realtype *);
//void PrintExpData(realtype *, int, realtype);
