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
#include "InitialConditions.h"
#include <chrono>
#include <string>

#include "TChem_Impl_ReactionRates.hpp"

#include "TChem_NetProductionRatePerMass.hpp"
#include "TChem_SpecificHeatCapacityPerMass.hpp"

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
//ZeroD class
// class myPb : public TCHEMPB{
//         public:
//         //members
//         ordinal_type  num_equations;
//         N_Vector Jac;
// };
class myPb : public TCHEMPB{
	public:
	//members
	ordinal_type  num_equations;
	N_Vector Jac;
	realtype t;
	realtype MaxStepTaken;
	realtype MinStepTaken;
	realtype ignTime;
	SUNMatrix	Mat;
	int 			Movie;
	std :: string	dumpJacFile;
	//functions
	ordinal_type	get_num_equations(void)
	{
		return this->num_equations;
	}
};
//OneD class
class myPb2{
	public:
	TCHEMPB pb;
	ordinal_type	num_equations;
	realtype		TubeLength;			//Marked for removal
	int 		NumGridPts;
	int 		vecLength;
	realtype 	delx;
	N_Vector 	Ghost;
	std::vector<N_Vector> Jacs;
	N_Vector 	Jac;
	N_Vector 	Tmp;
	//New grids for transport
	N_Vector	CpGrid;
	N_Vector	RhoGrid;
	N_Vector	DiffGrid;

	SUNMatrix	Mat;
	myPb2(ordinal_type, real_type_1d_view_type, WORK, int, N_Vector, realtype);
	~myPb2();//destructor
	void PrintGuts(void);				//Marked for removal
	void PrintJac();
	void SetGhost(N_Vector);
	void ScaleP(int, int);
	void SetVels(int, int);
	void SetVelAve();
	void UpdateOneDVel(N_Vector);
	void VelIntegrate(realtype *, N_Vector, realtype, realtype);
	void RunTests(N_Vector State);
	void CheckHeating(N_Vector State, realtype t);
	void SetLeftDiff(N_Vector State);
	void GhostChem(realtype t, N_Vector y, N_Vector ydot, void * pb);
	void SetGhostPVel(N_Vector y, int Experiment, int SampleNum, realtype VelVal);
	void SetAdvDiffReacPow(realtype, realtype, realtype, realtype, bool);
	void TempGradient(N_Vector State, N_Vector Gradient);
	void SetTransportGrid(N_Vector State);
	void SetGradient(N_Vector Input, N_Vector GradientVec, realtype Ghost);

	int			dumpJac;				//Do we want to dump the Jacobian
	std :: string	dumpJacFile;				//Where to dump it
	int			Movie;					//Do you want a movie?Y/n
	//Timers
	realtype	dataMoveTime;
	realtype	jacTime;
	realtype	jtv_Chem;
	realtype	jtv_Adv;
	realtype	jtv_Diff;
	realtype	rhsTime;
	realtype	rhs_Chem;
	realtype	rhs_Adv;
	realtype	rhs_Diff;
	realtype	jacMakeTime;
	//Lookup table timers
	realtype	rhsDiffLookTime;
	realtype	rhsChemLookTime;
	realtype	jtvDiffLookTime;
	realtype	jtvChemLookTime;
	//total integration time
	realtype	integrateTime;
	//Timestep data
	realtype	MaxStepTaken;
	realtype	MinStepTaken;
	//Storage vectors
	N_Vector	OMEGA;
	N_Vector	CP;
	N_Vector	CPPoly;
	N_Vector	TempTable;
	N_Vector	RhoTable;
	N_Vector	DiffTable;
	N_Vector	RHO;
	N_Vector	Vel;
	N_Vector	Scrap;
	N_Vector	VelScrap;
	N_Vector	SmallScrap;
	N_Vector	VelAve;
	realtype	PMultiplier;
	realtype	LeftDiff;
	N_Vector	SmallChem;
	N_Vector	MolarWeights;
	//Gradients
	N_Vector	CpGrad;
	N_Vector	RhoGrad;
	N_Vector	DiffGrad;
	N_Vector	TempGrad;

	//Parameter selectors
	realtype	Adv;
	realtype	Diff;
	realtype	React;
	realtype	Power;
	realtype	t;
	realtype	ignTime;
	bool		VelUp;
	int			FlameFrontLocation;
	realtype	HeatingRightGhost;
	int 		HeatingOn;
	//realtype	Lambda;
	KineticModelData	kmd;

	void VerifyHeatingExp(N_Vector State, N_Vector State0, realtype tElapsed);
	void CheckNaN(N_Vector, int);				//Checks the vector for NaN
	realtype ComputeCpTemp(N_Vector y);
	//void FetchTherm(N_Vector State);			//Sets cp and cpmm(the wanted quantity).
	void VerifyTempTable(N_Vector State);			//Verify Read in and run a temperature test.
	int  TempTableLookUp(realtype Temp, N_Vector TempTable);//Lookup a table
};

//End class stuff
//Definitions
//Old Prototype functions, where they are called is uncertain
int CHEM_RHS_TCHEM(realtype , N_Vector , N_Vector, void *);
int CHEM_COMP_JAC(N_Vector u, void* pb);
int CHEM_COMP_JAC_CVODE(realtype, N_Vector, N_Vector, SUNMatrix, void *, N_Vector, N_Vector, N_Vector);
int CHEM_JTV(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
//Zero-D Functions
int CHEM_RHS_TCHEM_V2(realtype , N_Vector , N_Vector, void *); 	//Rename required
int CHEM_COMP_JAC_V2(N_Vector u, void * pb);
int CHEM_COMP_JAC_CVODE_V2(realtype, N_Vector, N_Vector, SUNMatrix, void *, N_Vector, N_Vector, N_Vector); //Rename required
int CHEM_JAC_VOID(realtype, N_Vector, N_Vector, SUNMatrix, void *, N_Vector, N_Vector, N_Vector);
int CHEM_JTV_V2(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector); //Rename required


//Kapila problem functions
int RHS_KAPPA(realtype , N_Vector , N_Vector, void *);
int Jtv_KAPPA(N_Vector , N_Vector , realtype , N_Vector , N_Vector , void * , N_Vector);

//Helpers
void MatrixVectorProduct(int, realtype *, N_Vector, N_Vector, realtype *);

//One-D versions
	//RHS Functions
	int SUPER_RHS			(realtype, N_Vector, N_Vector, void *);
	int SUPER_CHEM_RHS_TCHEM	(realtype, N_Vector, N_Vector, void *);
	int SUPER_RHS_DIFF_CP		(realtype, N_Vector, N_Vector, void *);
	int SUPER_RHS_ADV_VEL		(realtype, N_Vector, N_Vector, void *); //Depreciated
	int SUPER_RHS_ADV_UPW		(realtype, N_Vector, N_Vector, void *);
	int SUPER_RHS_HEATING		(realtype, N_Vector, N_Vector, void *);

	//Jacobian function
	int SUPER_CHEM_JAC_TCHEM	(realtype, N_Vector, N_Vector, SUNMatrix, void *, N_Vector,
						N_Vector, N_Vector);

	//JtV functions
	int SUPER_JTV				(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
	int SUPER_CHEM_JTV			(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
	int SUPER_ADV_VEL_JTV		(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
	int SUPER_ADV_UPW_JTV		(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);
	int SUPER_DIFF_CP_JTV		(N_Vector, N_Vector, realtype, N_Vector, N_Vector, void*, N_Vector);


//One-D helpers
//int CleverMatVec(int, int, int, realtype*, realtype *, realtype *);
int SUPER_2_VEC(int, realtype *, realtype *, int, int);
int VEC_2_SUPER(int, realtype *, realtype *, int , int);

int SuperJac_2_Jac(int, realtype *, realtype*, int, int);
int Clean(int, realtype *);
int RescaleTemp(realtype, realtype *, int);
//Debugging
realtype CompareJacobians(void *, void *, N_Vector, N_Vector, N_Vector, N_Vector, SUNMatrix);

#endif
