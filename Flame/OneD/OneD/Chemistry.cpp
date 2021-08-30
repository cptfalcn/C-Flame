#include "Chemistry.h"


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
        //g_timer.reset();
        pbPtr->computeFunction(member, x ,f);
        //RHSTime+=g_timer.seconds();
        //RHSCnt++;
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
        //g_timer.reset();
        //auto Start=std::chrono::high_resolution_clock::now();
        pbPtr->computeJacobian(member, x ,J);
        //auto Stop=std::chrono::high_resolution_clock::now();
        //auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
        //JacTime+=Pass.count()/1e9;
        //JacTime+=g_timer.seconds();
        //JacCnt++;
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
        //g_timer.reset();
        MatrixVectorProduct(number_of_equations, JacD, v, tmp, JV);
        //MatVecTime+=g_timer.seconds();
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
	realtype *ADATA = SM_DATA_D(Jac);
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
        //g_timer.reset();
        //auto Start=std::chrono::high_resolution_clock::now();
        pbPtr->computeJacobian(member, x ,J);
	for (int i = 0 ; i < number_of_equations; i ++)
	{
		for (int j = 0 ; j < number_of_equations; j++)
			//ADATA[i+(number_of_equations-1)*j] = JacD[j+i*(number_of_equations-1)];//Original
			ADATA[i+(number_of_equations-1)*j] = JacD[j+i*(number_of_equations-1)];//Original
	}

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
