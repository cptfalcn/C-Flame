#ifndef CHEM_H
#define CHEM_H
// your declarations (and certain types of definitions) here
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
#define WORK TChem::KineticModelConstData<Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> >

using real_type = double;
using ordinal_type = int;
using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;

//=====================
//Prototypes & Classes
//=====================
//Add an additional class elements to pass to epic.

class myPb : public TCHEMPB{
        public:
        //members
        ordinal_type  num_equations;
        N_Vector Jac;
	//realtype *JacArr;
};

class myPb2{
	public:
	TCHEMPB pb;
	ordinal_type  num_equations;
	int NumGridPts;
	N_Vector Jac;
	SUNMatrix Mat;
	myPb2(ordinal_type, real_type_1d_view_type, WORK,int);
	~myPb2();//destructor
	void PrintGuts(void);
	void PrintJac();
	N_Vector Velocity;
};

//End class stuff
//Definitions
int CHEM_RHS_TCHEM(realtype , N_Vector , N_Vector, void *);
int CHEM_COMP_JAC(N_Vector u, void* pb);
int CHEM_COMP_JAC_CVODE(realtype, N_Vector, N_Vector, SUNMatrix, void *, N_Vector, N_Vector, N_Vector);
int CHEM_JAC_VOID(realtype, N_Vector, N_Vector, SUNMatrix, void *, N_Vector, N_Vector, N_Vector);
int CHEM_JTV(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
int RHS_KAPPA(realtype , N_Vector , N_Vector, void *);
int Jtv_KAPPA(N_Vector , N_Vector , realtype , N_Vector , N_Vector , void * , N_Vector);
void MatrixVectorProduct(int, realtype *, N_Vector, N_Vector, realtype *);

//One-D versions
	//RHS Functions
	int SUPER_CHEM_RHS_TCHEM	(realtype, N_Vector, N_Vector, void *);
	int SUPER_RHS_DIFF		(realtype, N_Vector, N_Vector, void *);
	int SUPER_RHS_ADV		(realtype, N_Vector, N_Vector, void *);

	//Jacobian function
	int SUPER_CHEM_JAC_TCHEM	(realtype, N_Vector, N_Vector, SUNMatrix, void *, N_Vector,
					N_Vector, N_Vector);

	//JtV functions
	int SUPER_CHEM_JTV		(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);

//One-D helpers
int CleverMatVec(int, int, int, realtype*, realtype *, realtype *);
int SUPER_2_VEC(int, realtype *, realtype *, int, int);
int VEC_2_SUPER(int, realtype *, realtype *, int , int);
int Jac_2_SuperJac(int, realtype*, realtype*, int, int);
int SuperJac_2_Jac(int, realtype *, realtype*, int, int);
int Clean(int, realtype *);
//Debugging
realtype CompareJacobians(void *, void *, N_Vector, N_Vector, N_Vector, N_Vector, SUNMatrix);

#endif
