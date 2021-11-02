#include "Chemistry.h"
//====================================================
//Note the following:
//All Sundials Matrices are stored column major.
//All Kokkos 2-d are stored row major.
//====================================================
//Class stuff
using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;

//==================================
//New problem class
//==================================
myPb2::myPb2(ordinal_type num_eqs, real_type_1d_view_type work, WORK kmcd, int GridPts)
{
	this->num_equations= num_eqs;
	this->pb;//create the problem
	this->pb._p = 101325;
	this->pb._work = work;
	real_type_1d_view_type fac("fac", 2*num_eqs);
	this->pb._fac = fac;
	this->pb._kmcd = kmcd;
	this->NumGridPts=GridPts;
	this->Jac=N_VNew_Serial(num_eqs*num_eqs*GridPts);
	this->Mat=SUNDenseMatrix(num_eqs*GridPts,num_eqs*GridPts);
	this->Velocity = N_VNew_Serial( (GridPts+1)); //This is staggered bigger
}

myPb2::~myPb2()
{
	N_VDestroy_Serial(this->Jac);
	SUNMatDestroy(this->Mat);
	N_VDestroy_Serial(this->Velocity);
}

void myPb2::PrintGuts(void)
{
	std :: cout << BAR << "General Parameters" << BAR << std :: endl;
	std :: cout << "Number of grid points & ghosts:" << this->NumGridPts << "\t\t";
	std :: cout << "Jacobians stored:" << N_VGetLength(this->Jac)/(this->num_equations*
					this->num_equations) << std :: endl;
}

//=====================================
//Zero-D Functions
//=====================================
int CHEM_RHS_TCHEM(realtype t, N_Vector u, N_Vector udot, void * pb)
{
        //TCHEMPB *pbPtr{ static_cast<TCHEMPB*>(pb)} ;//recast the type here.
        myPb * pbPtr{static_cast<myPb *> (pb)};//Recast
	ordinal_type number_of_equations = pbPtr->num_equations;
        realtype *y= NV_DATA_S(u);
        realtype *dy=NV_DATA_S(udot);
        /// we use std vector mimicing users' interface
        using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
        using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
        // output rhs vector and x via a kokkos wrapper.
        real_type_1d_view_type x(y,     number_of_equations);
        real_type_1d_view_type f(dy,    number_of_equations);
        //=================================
        // Compute right hand side vector
        //=================================
        auto member =  Tines::HostSerialTeamMember();
        pbPtr->computeFunction(member, x ,f);
        return 0;
}

        //===============================================================
        //   _____   ___     ____   ___    _____   _____    ___   __   _
        //  |__ __| /   \   / ___) / _ \  |     \ |__ __|  /   \ |  \ | |
        //  _ | |  | (x) | | /    / / \ \ |  x  /   | |   | (x) ||   \| |
        // / (| |  |  n  | | \___ \ \_/ / |  x  \  _| |   |  n  || |\   |
        // \____/  |_| |_|  \____) \___/  |_____/ |_____| |_| |_||_| \__|
        //===============================================================
int CHEM_COMP_JAC(N_Vector u, void* pb)
{
        //problem_type problem;
        myPb * pbPtr{static_cast<myPb *> (pb)};//Recast
        ordinal_type number_of_equations=pbPtr->num_equations;//Get number of equations
        //======================================
        //Set the necessary vectors and pointers
        //======================================
        realtype *y= NV_DATA_S(u);
        realtype *JacD= NV_DATA_S(pbPtr->Jac);
        //=============================================
        //We use std vector mimicing users' interface
        //=============================================
        using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
        using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
        using real_type_2d_view_type = Tines::value_type_2d_view<real_type,host_device_type>;

        //============================================
        // output rhs vector and J matrix wrappers.
        //============================================

        real_type_1d_view_type x(y   ,   number_of_equations);
        real_type_2d_view_type J(JacD,   number_of_equations, number_of_equations);

        //===============================
        /// Compute Jacobian
        //===============================
        auto member =  Tines::HostSerialTeamMember();
        pbPtr->computeJacobian(member, x ,J);
        return 0;
}

int CHEM_JTV(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void* pb, N_Vector tmp)
{
        //problem_type problem;
        myPb * pbPtr{static_cast<myPb *> (pb)};//Recast
        ordinal_type number_of_equations=pbPtr->num_equations;//Get number of equations
        //======================================
        //Set the necessary vectors and pointers
        //======================================
        realtype *y= NV_DATA_S(u);
        realtype *JacD= NV_DATA_S(pbPtr->Jac);
        realtype *JV=NV_DATA_S(Jv);
        //===================
        //Compute JV
        //===================
        MatrixVectorProduct(number_of_equations, JacD, v, tmp, JV);
        return 0;
}

int CHEM_JAC_VOID(realtype t, N_Vector u, N_Vector fy, SUNMatrix Jac,
               void * pb, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	return 0;
}

int CHEM_COMP_JAC_CVODE(realtype t, N_Vector u, N_Vector fy, SUNMatrix Jac,
               void * pb, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	//problem_type problem;
	myPb * pbPtr{static_cast<myPb *> (pb)};//Recast
        ordinal_type number_of_equations=pbPtr->num_equations;//Get number of equations
        //======================================
        //Set the necessary vectors and pointers
        //======================================
        realtype *y= NV_DATA_S(u);
        realtype *JacD= NV_DATA_S(pbPtr->Jac);
	realtype *SUNJACPTR = SM_DATA_D(Jac);
        //=============================================
        //We use std vector mimicing users' interface
        //=============================================
        using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
        using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
        using real_type_2d_view_type = Tines::value_type_2d_view<real_type,host_device_type>;

        //============================================
        // output rhs vector and J matrix wrappers.
        //============================================
        real_type_1d_view_type x(y   ,   number_of_equations);
        real_type_2d_view_type J(JacD,   number_of_equations, number_of_equations);

        //===============================
        /// Compute Jacobian
        //===============================
        auto member =  Tines::HostSerialTeamMember();
        pbPtr->computeJacobian(member, x ,J);//Do not transpose, let Jtv handle this.
	/**/
	for (int i = 0 ; i < number_of_equations; i ++)//Need to swap orientation
	{//row
		for (int j = 0 ; j < number_of_equations; j++)//column
			SUNJACPTR[i+j*number_of_equations] = JacD[j + i *number_of_equations];//Swap
	}
	/**/
	return 0;
}

int RHS_KAPPA(realtype t, N_Vector u, N_Vector udot, void *userData)
{
        realtype EPS=1e-2;
        realtype *y = NV_DATA_S(u);
        realtype *dy = NV_DATA_S(udot);
        double Omega= .5 * y[0] * y[1] * exp ( (y[2]-1.0 ) /  (EPS * y[2]  ) );
        dy[0] = -Omega;
        dy[1] = Omega-y[1];
        dy[2] = y[1];
//      cout << dy[0] <<"\t\t" << dy[1] << "\t\t" << dy[2] <<endl;
        return 0;
}

int Jtv_KAPPA(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void *userData, N_Vector tmp)
{
        realtype EPS=1e-2;
        realtype *JV = NV_DATA_S(Jv);
        realtype *Y = NV_DATA_S(u);
        realtype *V = NV_DATA_S(v);
        double EXP= exp ( (Y[2] - 1.0 ) / (EPS * Y[2]));
        double DEXP=1/(EPS * Y[2]*Y[2]);
        JV[0] = -.5*EXP* (V[0]*Y[1]+V[1]*Y[0]+V[2]*Y[0]*Y[1]*DEXP);
        JV[1] =-JV[0]-V[1];
        JV[2] = V[1];
        return 0;
}


/*
//Matrix Vector Product
//Needs to be rewritten
//=============================================================
//  __ __    ___   _____      __   __ ___   ___
// |  V  | ./ _ \.|_   _| ___ \ \ / /| __) / __)
// | .V. | | (_) |  | |  |___| \ v / | __)( (__
// |_| |_| |_/ \_|  |_|         \_/  |___) \___)
//============================================================
//We cover our bases by doing MatVec directly for EPI methods.
//By doing this, we don't have to transpose the data.
*/
void MatrixVectorProduct(int number_of_equations, realtype * JacD, N_Vector x, N_Vector tmp, realtype * JV)
{
        realtype * TMP  = NV_DATA_S(tmp);
	realtype * X 	= NV_DATA_S(x);
        for (int i=0; i< number_of_equations; i++)//for each state
        {
                for(int j=0; j<number_of_equations; j++)//marches across the column
                {
                        TMP[j]=JacD[j+i*number_of_equations];//Stored row major
			//JV[i] += X[j] * JacD[i * number_of_equations + j];
                }
                JV[i]=N_VDotProd(tmp,x);
        }
}


//=======================================================
//  ___  _  _  ___  ___  ___
// / __/| || || _ \|  _|| _ \
// \__ \| || ||  _/|  _||   <
// /___/ \__/ |_|  |___||_|_|
//
//========================================================
//Bugged, but it runs without throwing an error.
//These use the second version of the problem class, myPb2.
//========================================================
int SUPER_CHEM_RHS_TCHEM(realtype t, N_Vector State, N_Vector StateDot, void * pb)
{
	myPb2 * pbPtr{static_cast<myPb2 *> (pb)};//Recast
	int Length = N_VGetLength(State);
	int num_eqs = pbPtr->num_equations;
	int grid_sz = pbPtr->NumGridPts;
	N_Vector yTemp = N_VNew_Serial(num_eqs);
	N_Vector fTemp = N_VNew_Serial(num_eqs);
	realtype * DATA = NV_DATA_S(yTemp);
	realtype * STATEDATA = NV_DATA_S(State);
	realtype * FDATA = NV_DATA_S(fTemp);
	realtype * STATEDOTDATA = NV_DATA_S(StateDot);

	for(int i = 0 ; i < grid_sz ; i ++ )
	{//March over copies/grid points and grab a data set.
		SUPER_2_VEC(i, DATA, STATEDATA, num_eqs, grid_sz);//Appears to be working
        	using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
        	using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
        	// output rhs vector and x via a kokkos wrapper.
        	real_type_1d_view_type x(DATA,     num_eqs);
        	real_type_1d_view_type f(FDATA,	    num_eqs);
        	//===============================
        	// Compute right hand side vector
        	//===============================
        	auto member =  Tines::HostSerialTeamMember();
        	pbPtr->pb.computeFunction(member, x ,f);
		VEC_2_SUPER(i, FDATA, STATEDOTDATA, num_eqs, grid_sz);//Set FDATA into STATEDOTDATA
	}
	N_VDestroy_Serial(yTemp);
	N_VDestroy_Serial(fTemp);
	return 0;
}

//===============================
//Super Jac Via TChem
//===============================
int SUPER_CHEM_JAC_TCHEM(realtype t, N_Vector State, N_Vector StateDot, SUNMatrix Jac, void * pb, N_Vector tmp1,
				N_Vector tmp2, N_Vector tmp3)
{
	myPb2 * pbPtr{static_cast<myPb2 *> (pb)};//Recast
        int Length = N_VGetLength(State);
        int num_eqs = pbPtr->num_equations;
        int grid_sz = Length/num_eqs; //This give a grid size/ # copies
        N_Vector yTemp = N_VNew_Serial(num_eqs);
	N_Vector JacTemp = N_VNew_Serial(num_eqs*num_eqs);
        realtype * DATA = NV_DATA_S(yTemp);
        realtype * STATEDATA = NV_DATA_S(State);
        realtype * JACDATA = NV_DATA_S(pbPtr->Jac);
        realtype * STATEDOTDATA = NV_DATA_S(StateDot);
	realtype * TEMPJACDATA = NV_DATA_S(JacTemp);
	using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
        using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
        using real_type_2d_view_type = Tines::value_type_2d_view<real_type,host_device_type>;
	int Block = num_eqs * num_eqs ;

        for(int i = 0 ; i < grid_sz ; i ++ )
        {//March over copies/grid points and grab the data needed
		SUPER_2_VEC(i, DATA, STATEDATA, num_eqs, grid_sz);
        	//output Jac vector and x via a kokkos wrapper.
		real_type_1d_view_type x(DATA,     num_eqs);
        	real_type_2d_view_type J(TEMPJACDATA,     num_eqs, num_eqs);
		//=================================
		// Compute Jacobian
        	//=================================
        	auto member =  Tines::HostSerialTeamMember();
        	pbPtr->pb.computeJacobian(member, x ,J);
		//Jac_2_SuperJac(i, JACDATA, TEMPJACDATA, num_eqs, grid_sz);//Bugged
		/**/
		for( int j = 0 ; j < num_eqs; j ++)
                	for(int k = 0; k < num_eqs; k++)
                        	JACDATA[i*Block + j * num_eqs + k] = TEMPJACDATA[j*num_eqs + k];
		/**/
        }

        N_VDestroy_Serial(yTemp);
	N_VDestroy_Serial(JacTemp);
	return 0;
}
//=========================================
//Computes the TChem JtV in Sundials format
//=========================================
int SUPER_CHEM_JTV(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void* pb, N_Vector tmp)
{
	myPb2 * pbPtr{static_cast<myPb2 *> (pb)};//Recast
	//Set necessary temp Data
	int num_eqs		= pbPtr->num_equations;
	int num_grid		= pbPtr->NumGridPts;
	realtype * JACDATA	= NV_DATA_S(pbPtr->Jac);
	realtype * VDATA 	= NV_DATA_S(v);
	realtype * JVDATA 	= NV_DATA_S(Jv);
	int Block 		= num_eqs*num_eqs;
	N_Vector SmallJac 	= N_VNew_Serial(Block);//It might be possible to use "tmp".
	realtype * SmallData 	= NV_DATA_S(SmallJac);
	N_Vector SmallJV 	= N_VNew_Serial(num_eqs*num_grid);
	realtype * SmallJVData	= NV_DATA_S(SmallJV);
	N_Vector SmallV		= N_VNew_Serial(num_eqs*num_grid);
	realtype * SmallVData	= NV_DATA_S(SmallV);
	N_Vector SmallTmp 	= N_VNew_Serial(num_eqs*num_grid);

	for( int i = 0; i < num_grid; i ++ )
	{//Over each grid
		SuperJac_2_Jac(i, SmallData, JACDATA, num_eqs, Block);
		SUPER_2_VEC(i, SmallJVData, JVDATA, num_eqs, num_grid);
		SUPER_2_VEC(i, SmallVData, VDATA, num_eqs, num_grid);
		MatrixVectorProduct(num_eqs, SmallData, SmallV, SmallTmp, SmallJVData);
		VEC_2_SUPER(i, SmallJVData, JVDATA, num_eqs, num_grid);
		//CleverMatVec(i, num_eqs, Block, JVDATA, VDATA, JACDATA);
	}
	N_VDestroy_Serial(SmallJac);
	N_VDestroy_Serial(SmallTmp);
	N_VDestroy_Serial(SmallJV);
	return 0;
}

//============================
//Vector Shenanagins
//==================================================================================================
//Do the clever MatVec product
//Tested on two problems, Jac is block diagonal ones, and Jac is block diagonal 0,1,2,...,n in rows.
//==================================================================================================
int CleverMatVec(int i, int num_eqs, int Block, realtype * JV, realtype * V, realtype * Jac)
{
	for(int j = 0 ; j < num_eqs; j ++)//row
		for( int k = 0; k < num_eqs; k++)//column
			JV[i * num_eqs + j ] += V[i * num_eqs + k] * Jac[i*Block + j * num_eqs + k];
	return 0;
}

//=================
//SUPERRHS_2_RHS
//=================
int  SUPER_2_VEC(int i, realtype * VECDATA, realtype* SUPERDATA, int num_eqs, int grid_sz)
{
	for(int j = 0; j < num_eqs; j++)//March over each equation of the state.
		VECDATA[j] = SUPERDATA[i + j * grid_sz];
	return 0;
}

//Vec_2_Super
int VEC_2_SUPER(int i, realtype * VECDATA, realtype *SUPERDATA, int num_eqs, int grid_sz)
{
	for(int j = 0; j < num_eqs; j++)
		SUPERDATA[i + j * grid_sz] = VECDATA[j];//Working?
		//SUPERDATA[j + i * grid_sz] = VECDATA[j];
	return 0;
}

//Jac into SuperJac
int Jac_2_SuperJac(int i, realtype * Jac, realtype * SuperJac, int num_eqs, int grid_sz)
{
	int Block = num_eqs * grid_sz;
	for( int j = 0 ; j < num_eqs; j ++)
		for(int k = 0; k < num_eqs; k++)
			SuperJac[i*Block + j * num_eqs + k] = Jac[j*num_eqs + k];
	return 0;
}

int SuperJac_2_Jac(int i, realtype * Jac, realtype * SuperJac, int num_eqs, int Block)
{
	for ( int j = 0 ; j < num_eqs; j ++)
		for (int k = 0 ; k < num_eqs; k ++)
			Jac[j*num_eqs + k] = SuperJac[i * Block + j * num_eqs + k];
	return 0;
}
//=======================================
//Clean the vector of non-physical values
//=======================================
int Clean(int length, realtype * Data)
{
	for( int i = 0 ; i < length ; i ++)
	{
		if(Data[i]<0)
			Data[i]=0;
	}
	return 0;
}


//============================================
//General Debugging functions
//=============================================
//Error Checking the Jacobians
//Send in the two Jacobian functions and compare their two norm.
// if it is off by a large amount, something is up... Something is up.
//======================================================
realtype CompareJacobians(void * ZeroDVersion, void * OneDVersion, N_Vector y, N_Vector State,
			N_Vector yDot, N_Vector StateDot, SUNMatrix SuperJac)
{	//Get this figured out, have a seg fault do to the Error calculation loop.
	realtype Error=0;
	realtype t=0;
	CHEM_COMP_JAC(y, ZeroDVersion);
	SUPER_CHEM_JAC_TCHEM(t, State, State, SuperJac, OneDVersion, StateDot, StateDot, StateDot);
	CHEM_RHS_TCHEM(t, y, yDot, ZeroDVersion);
	SUPER_CHEM_RHS_TCHEM(t, State, StateDot, OneDVersion);
	myPb2 * OneDProblem{static_cast<myPb2 *> (OneDVersion)};//Recast
	myPb *  ZeroDProblem{static_cast<myPb *> (ZeroDVersion)};//Recast
	std :: cout << "ZeroD Max Jac Entry: " << N_VMaxNorm(ZeroDProblem->Jac);
	std :: cout << "\t\tOneD Max Jac Entry: "  << N_VMaxNorm(OneDProblem->Jac)<< std :: endl;
	std :: cout << "ZeroD Max Rhs Entry: " << N_VMaxNorm(yDot);
	std :: cout << "\t\tOneD Max Rhs Entry: " << N_VMaxNorm(StateDot) << std :: endl;
	//Define some more stuff now that the problems have been type cast
	int num_eqs=  OneDProblem->num_equations;
	int num_pts = OneDProblem->NumGridPts;
	int vecLength = num_eqs*num_pts;
	int Block = vecLength * vecLength;
	N_Vector StretchedRhs = N_VNew_Serial(vecLength);
	N_Vector StretchedJac = N_VNew_Serial(Block);
	realtype * StretchedData = NV_DATA_S(StretchedRhs);
	realtype * StretchedJacData= NV_DATA_S(StretchedJac);
	realtype * RhsData = NV_DATA_S(yDot);
	realtype * OneDJacD= NV_DATA_S(OneDProblem->Jac);
	realtype * ZeroDJacD=NV_DATA_S(ZeroDProblem->Jac);
	if (num_pts != 1)
	{
		for( int i = 0 ; i < num_pts ; i ++ )
		{
			VEC_2_SUPER(i, RhsData, StretchedData, num_eqs, num_pts);
			Jac_2_SuperJac(i, ZeroDJacD, OneDJacD, num_eqs, num_pts);
		}
	}
	N_Vector ErrJac = N_VNew_Serial(Block);
	N_Vector ErrRhs = N_VNew_Serial(vecLength);
	//std :: cout << "Number of grid points: " << num_pts << std :: endl;
	if (num_pts != 1 )
	{
		N_VLinearSum(1, StretchedRhs, -1, StateDot, ErrRhs);
	}
	else
	{
		N_VLinearSum(1, yDot, -1, StateDot, ErrRhs);
		N_VLinearSum(1, OneDProblem->Jac, -1, ZeroDProblem->Jac, ErrJac);
	}
	realtype * DummyJacD= NV_DATA_S(ErrJac);
	if(N_VMaxNorm(ErrRhs) != 0 || N_VMaxNorm(ErrJac) !=0)
	{
		std :: cout << "Max Rhs Error: " << N_VMaxNorm(ErrRhs) << "\t\t";
		std :: cout << "Max Jacobian Error: " << N_VMaxNorm(ErrJac) << std :: endl;
	}
	else
	{
		std :: cout << "These quantities agree completely!" << std :: endl;
	}
	//Take out the trash
	N_VDestroy_Serial(ErrJac);
	N_VDestroy_Serial(ErrRhs);
	return 0;
}

