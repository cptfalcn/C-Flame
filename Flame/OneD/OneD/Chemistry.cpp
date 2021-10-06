#include "Chemistry.h"
//====================================================
//Note the following:
//All Sundials Matrices are stored column major.
//All Kokkos 2-d are stored row major.
//====================================================
//Class stuff
using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;

//New problem class
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
}

myPb2::~myPb2()
{
	N_VDestroy_Serial(this->Jac);
}

void myPb2::PrintGuts(void)
{
	std :: cout << BAR <<"Zero-D Ignition Parameters" << BAR << std :: endl;
	std :: cout << "pressure: " << this->pb._p <<  "\t\t";
        std :: cout << "workspace: " << sizeof(this->pb._work) << std :: endl;
        std :: cout << "Fac ??: " << sizeof(this->pb._fac) << "\t\t" ;
        std :: cout << "Kinetic model: " << std :: endl;
	std :: cout << BAR << "General Parameters" << BAR << std :: endl;
	std :: cout << "Number of grid points & ghosts:" << this->NumGridPts << "\t\t";
	std :: cout << "Jacobians stored:" << N_VGetLength(this->Jac)/(this->num_equations*
					this->num_equations) << std :: endl;
}


//Functions
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
        realtype *TMP=NV_DATA_S(tmp);
        //===================
        //Compute JV
        //===================
        MatrixVectorProduct(number_of_equations, JacD, v, tmp, JV);
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
        pbPtr->computeJacobian(member, x ,J);
	/*
	for (int i = 0 ; i < number_of_equations; i ++)//Need to swap orientation
	{//row
		for (int j = 0 ; j < number_of_equations; j++)//column
			SUNJACPTR[i+j*number_of_equations] = JacD[j + i *number_of_equations];//Swap
	}
	*/
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
        for (int i=0; i< number_of_equations; i++)//for each state
        {
                for(int j=0; j<number_of_equations; j++)//marches across the column
                {
                        TMP[j]=JacD[j+i*number_of_equations];//Stored row major
                }
                JV[i]=N_VDotProd(tmp,x);
        }
}

//========================================================
//Not yet tested, but it runs without throwing an error.
//These use the second version of the problem class, myPb2.
//========================================================
int SUPER_CHEM_RHS_TCHEM(realtype t, N_Vector State, N_Vector StateDot, void * pb)
{
	myPb2 * pbPtr{static_cast<myPb2 *> (pb)};//Recast
	int Length = N_VGetLength(State);
	int num_eqs = pbPtr->num_equations;
	int grid_sz = Length/num_eqs; //This give a grid size/ # copies
	N_Vector yTemp = N_VNew_Serial(num_eqs);
	N_Vector fTemp = N_VNew_Serial(num_eqs);
	realtype * DATA = NV_DATA_S(yTemp);
	realtype * STATEDATA = NV_DATA_S(State);
	realtype * FDATA = NV_DATA_S(fTemp);
	realtype * STATEDOTDATA = NV_DATA_S(StateDot);

	for(int i = 0 ; i < grid_sz ; i ++ )
	{//March over copies/grid points and grab the data needed
		for(int j = 0; j < num_eqs; j++)
		{//March over each equation, State is organized by supervector
			DATA[j]=STATEDATA[i + j * grid_sz];//std :: cout << DATA[j] << "\t\t";
		}
	//run TCHEM to get the point's RHS
        using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
        using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
        // output rhs vector and x via a kokkos wrapper.
        real_type_1d_view_type x(DATA,     num_eqs);
        real_type_1d_view_type f(FDATA,	    num_eqs);
        //=================================
        // Compute right hand side vector
        //=================================
        auto member =  Tines::HostSerialTeamMember();
        pbPtr->pb.computeFunction(member, x ,f);
	//This is sorted by copies of state->STATEDOTDATA[  FData[j]<-This is entire state aap
		for(int j = 0; j < num_eqs; j++)
			STATEDOTDATA[i + j * grid_sz] = FDATA[j];

	}
	//std :: cout << std :: endl;
	//std :: cout << BAR << "FDATA" << BAR << std :: endl;
	/*
	for(int i = 0 ; i <num_eqs; i++)
		std:: cout << FDATA[i] << "\t\t";
	std :: cout << std :: endl << BAR << "Testing the SuperVector" << BAR << std :: endl;
	for( int i = 0 ; i < num_eqs * grid_sz; i++)
		std :: cout << STATEDATA[i] << "\t\t";
	std :: cout << std :: endl << BAR << "State Derivative" << BAR << std :: endl;
	for( int i = 0 ; i < num_eqs *grid_sz; i++)
		std :: cout << STATEDOTDATA[i] << "\t\t";
	std :: cout << std :: endl;
	*/
	N_VDestroy_Serial(yTemp);
	N_VDestroy_Serial(fTemp);
	return 0;
}

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

        for(int i = 0 ; i < grid_sz ; i ++ )
        {//March over copies/grid points and grab the data needed
                for(int j = 0; j < num_eqs; j++)//March over each equation of the state.
                        DATA[j]=STATEDATA[i + j * grid_sz];//std :: cout << DATA[j] << "\t\t";

        	//output Jac vector and x via a kokkos wrapper.
        	real_type_1d_view_type x(DATA,     num_eqs);
        	real_type_2d_view_type J(TEMPJACDATA,     num_eqs, num_eqs);
		//=================================
		// Compute Jacobian
        	//=================================
        	auto member =  Tines::HostSerialTeamMember();
        	pbPtr->pb.computeJacobian(member, x ,J);

		for (int j = 0 ; j < num_eqs; j ++)
        	{//row
                	for (int k = 0 ; k < num_eqs; k++)//column
                        	JACDATA[i * num_eqs * num_eqs + k + j * num_eqs]
						= TEMPJACDATA[ k + j * num_eqs];
        	}
        }

        N_VDestroy_Serial(yTemp);
	return 0;
}

int SUPER_CHEM_JTV(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void* pb, N_Vector tmp)
{
	myPb2 * pbPtr{static_cast<myPb2 *> (pb)};//Recast
	//Set necessary temp Data
	//SmallJac=NV_New_Serial(num_eqs);

	//for( int i = 0; i < num_grid; i ++ )
	{
		//Do small Jac JtV
	}
	return 0;
}
