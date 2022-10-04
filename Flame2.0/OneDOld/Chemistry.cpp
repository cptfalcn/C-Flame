#include "Chemistry.h"
//====================================================
//Note the following:
//All Sundials Matrices are stored column major.
//All Kokkos 2-d are stored row major.
//====================================================
//Class stuff
using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
using real_type_2d_view_type = Tines::value_type_2d_view<real_type,host_device_type>;

//==================================
//New problem class
//============================================================================================
//  .d8888b.                             888                              888                    
// d88P  Y88b                            888                              888                    
// 888    888                            888                              888                    
// 888         .d88b.  88888b.  .d8888b  888888 888d888 888  888  .d8888b 888888 .d88b.  888d888 
// 888        d88""88b 888 "88b 88K      888    888P"   888  888 d88P"    888   d88""88b 888P"   
// 888    888 888  888 888  888 "Y8888b. 888    888     888  888 888      888   888  888 888     
// Y88b  d88P Y88..88P 888  888      X88 Y88b.  888     Y88b 888 Y88b.    Y88b. Y88..88P 888     
//  "Y8888P"   "Y88P"  888  888  88888P'  "Y888 888      "Y88888  "Y8888P  "Y888 "Y88P"  888     
//===========================================================================================                                                                                                                                                                                                            
myPb2::myPb2(ordinal_type num_eqs, real_type_1d_view_type work, WORK kmcd, int GridPts,
			N_Vector y,realtype Delx)
{
	this->num_equations= num_eqs;
	//Create TChem Pb
	this->pb;								//create the problem
	this->pb._p 	= 101325;				//Standard (1 atm) pressure
	this->pb._work	= work;
	this->pb._kmcd 	= kmcd;
	this->GasConst  = 8.31446261815324;		//wikipedia gas constant
	//End TChem problem
	//Critical parameters 
	this->NumGridPts= GridPts;
	this->vecLength = num_eqs * GridPts;
	sunindextype Block = this->vecLength * this->vecLength;
	real_type_1d_view_type fac("fac", 2*vecLength);
	this->pb._fac   = fac;
	this->delx 		= Delx;
	this->LeftDiff	= 0;

	//Large N_Vectors
	this->Jac		= N_VNew_Serial(Block);				//Remove this
	//this->Mat		= SUNDenseMatrix(num_eqs*GridPts,num_eqs*GridPts);//Change to some tiny number
	this->Mat		= SUNDenseMatrix(1,1);				//Hard code to a small number
	this->Ghost 	= N_VNew_Serial(num_eqs);			//One hanging ghost NEQ
	this->PointScrap= N_VClone(this->Ghost);			//Scrap work for a single pnt of num_eqs
	//Set the ghost
	N_VScale(1.0 , y, this->Ghost);						//Use y to get the ghost points.

	//Velocity data
	this->VelAve	= N_VNew_Serial(GridPts);			//Averaged velocity
	this->Vel 		= N_VNew_Serial(GridPts+1);		//This is staggered bigger
	//Data Tables:	All fixed sizes
	this->CPPoly	= N_VNew_Serial(500);				//CP table
	this->TempTable	= N_VNew_Serial(500);				//Temperature reference table
	this->RhoTable	= N_VNew_Serial(500);				//Rho table
	this->DiffTable	= N_VNew_Serial(500*num_eqs);		//Each scalar has its own table
	this->RhoDiffTab= N_VNew_Serial(500*num_eqs);		//Each scalar has its own table
	//Scratch work vectors
	this->SmallScrap= N_VNew_Serial(GridPts);			//Used for any small grid
	this->VelScrap	= N_VNew_Serial(GridPts+1);			//Used for any velocity grid
	this->SmallChem	= N_VNew_Serial(num_eqs);			//Used for a single point of chem
	this->Tmp 		= N_VNew_Serial(vecLength);			//TempScratch vector of full size
	this->Scrap 	= N_VNew_Serial(vecLength);			//Temp storage, clean before accessing.
	this->Scrap2 	= N_VClone(this->Scrap);			//Another temp storage vector.
	this->ScalarGradient= N_VClone(this->Scrap);
	//Options
	this->HeatingOn = 1;
	this->dumpJac	= 0;								//Do we dump the Jac for visualization
	//Transport grids
	this->CpGrid	= N_VNew_Serial(GridPts);		//The value of Cp based on grid position.
	this->RhoGrid	= N_VNew_Serial(GridPts);		//The value of Rho based on grid position.
	std :: cout << "Vec length: " << vecLength << std :: endl;
	this->DiffGrid		= N_VClone(this->Tmp);
	//this->DiffGrid	= N_VNew_Serial(vecLength);			//Diffusion coefficient based on grid position.
	this->MolarWeights = N_VNew_Serial(GridPts);		//Only grab the species molar wieghts;
	//Gradient Grids
	this->CpGrad	= N_VClone(this->CpGrid);
	this->RhoGrad	= N_VClone(this->RhoGrid);
	this->DiffGrad	= N_VClone(this->DiffGrid);
	this->TempGrad	= N_VClone(this->CpGrid);

	//Jacobian Array and prototype
	N_Vector LittleJac	= N_VNew_Serial(num_eqs*num_eqs);//JacArray Prototype

	for(int i = 0; i < GridPts; i ++)
	{
		this->Jacs.push_back(N_VClone(LittleJac));
	}

	this->GasWeight		= N_VClone(this->CpGrid);
	this->MolarWeights 	= N_VClone(this->Ghost);
	realtype * MWPtr 	= NV_DATA_S(this->MolarWeights);
	MWPtr[0] 			= 0;
	for(int i = 0; i < num_eqs - 1; i++)
	{
		//std :: cout << this->pb._kmcd.sMass[i] << std :: endl;
		MWPtr[i+1] = kmcd.sMass[i];
		//std :: cout << MWPtr[i+1] << std :: endl;
	}
}

//====================================================
//  #     #                                           
//  #     # ###### #      #####  ###### #####   ####  
//  #     # #      #      #    # #      #    # #      
//  ####### #####  #      #    # #####  #    #  ####  
//  #     # #      #      #####  #      #####       # 
//  #     # #      #      #      #      #   #  #    # 
//  #     # ###### ###### #      ###### #    #  ####  
//=========================================================
//Run after construction to set problem specific parameters
//=========================================================
//============================================================================//
//The RHS and JtV requires constants to scale their term as follows:		  //
// AdvPart * adv + DiffPart * diff + ReacPart * reac + Heating*pow			  //
//Additionally we set if the Velocity is update or not.						  //
//depreciation:																  //
//The VelUp will be removed in the future so that it is always being updated. //
//============================================================================//
void myPb2::SetAdvDiffReacPow(realtype adv, realtype diff, realtype reac, realtype pow, bool up)
{
	this->Adv=adv;
	this->Diff=diff;
	this->React=reac;
	this->Power=pow;
	this->VelUp=up;
}

//============================================//
//Sets the ghost point based off the vector y.
//Main uses the initial condition y
//============================================//
void myPb2::SetGhost(N_Vector y)
{
	N_VScale(1.0, y, this->Ghost);
}



//=================
//Destructor
//Called by main
//=================
myPb2::~myPb2()
{
	N_VDestroy_Serial(this->Jac);
	SUNMatDestroy(this->Mat);
	N_VDestroy_Serial(this->Vel);
	N_VDestroy_Serial(this->Ghost);
	N_VDestroy_Serial(this->Tmp);
	N_VDestroy_Serial(this->SmallChem);
	N_VDestroy_Serial(this->SmallScrap);//Problematic line?
	N_VDestroy_Serial(this->VelScrap);
	N_VDestroy_Serial(this->Scrap);
	N_VDestroy_Serial(this->TempTable);
	N_VDestroy_Serial(this->RhoTable);
	N_VDestroy_Serial(this->CPPoly);
	N_VDestroy_Serial(this->DiffTable);
	N_VDestroy_Serial(this->CpGrid);
	N_VDestroy_Serial(this->RhoGrid);
	N_VDestroy_Serial(this->DiffGrid);
	for( int i = 0 ; i < this->NumGridPts; i ++)
		N_VDestroy(Jacs[i]);
}

//==============================
//Scale the pressure
//==============================
//Different experiment/samples can be run at different pressures.
//The higher the pressure the more energy.
//	Called by main after contructor during setup.
//================================================================
void myPb2::ScaleP(int Experiment, int Sample)
{
	realtype scale = 1.0;
	if(Experiment == 2 && Sample == 4)
		scale = 4.0;
	this->pb._p=scale*this->pb._p;
}

//===============================
//Create the Velocity grids
//===============================
//This is run during setup.  The user can provide a constant velocity field value.
//We set the velocity field to that, then call Set VelAve.
//		SetVelAve is called after an integration step
//		This gives access to the velocity at the cell center.
//		Called by main during setup
//==================================================================================
void myPb2::SetVels(int velLength, int velOpt)
{
	realtype * VelData 		= NV_DATA_S(this->Vel);
	realtype * VelAveData	= NV_DATA_S(this->VelAve);
	N_VConst(velOpt, this->Vel);
	this->SetVelAve();
}
//Attach the RHS function handles
void myPb2::Set_RHS(CVRhsFn OneD_RHS_Chem, CVRhsFn OneD_RHS_Adv, CVRhsFn OneD_RHS_Diff, CVRhsFn OneD_RHS_Heat)
{
	this->RHS_Adv		=OneD_RHS_Adv;
	this->RHS_Diff		=OneD_RHS_Diff;
	this->RHS_Chem		=OneD_RHS_Chem;
	this->RHS_Heat 		=OneD_RHS_Heat;
}

//Attach the JtV function handles.
void myPb2::Set_Jtv(CVLsJacTimesVecFn OneD_JtV_Adv, CVLsJacTimesVecFn OneD_Jtv_Chem, CVLsJacTimesVecFn OneD_JtV_CrossDiff,CVLsJacTimesVecFn OneD_JtV_Diff)
{
	this->JtV_Adv		=  	OneD_JtV_Adv;
	this->JtV_Diff		=  	OneD_JtV_Diff;
	this->JtV_CrossDiff	= 	OneD_JtV_CrossDiff;
	this->JtV_Chem		= 	OneD_Jtv_Chem;
	

}

//===============================
//Set Vel Average
//===============================
//Velocity is staggered above the scalar grid at its boundaries.
//So we have to average the velocity onto the scalar grid.
//The spacing of the velocity is superior: ith and ith+1 velocity pts associate to ith scalar 
//		Called by : MAIN after each velocity update step.
void myPb2::SetVelAve()
{
	realtype * VelData 		= N_VGetArrayPointer(this->Vel);
	realtype * VelAveData	= N_VGetArrayPointer(this->VelAve);
	for(int i = 0 ; i < this->NumGridPts; i++)		//Before used GridPts-1
		VelAveData[i]   	= (VelData[i] + VelData[i+1])/2;
}


//===========================
//Get TempGradient
//===========================
//Note here we step up to the Velocity grid.
//The Scalar grid is inferior to the velocity grid.
//So the scalar ith and ith-1 associate to the ith velocity.
//Because of the inferiority, this prelcudes the boundary pts.
//Called by: myPb2::UpdateOneDVel (in a comment)
//Status: Dangling
void myPb2::TempGradient(N_Vector State, N_Vector Gradient)
{
	//Declare
	realtype * SD		= N_VGetArrayPointer(State);
	realtype * GD		= N_VGetArrayPointer(this->Ghost);
	realtype * GradD	= N_VGetArrayPointer(Gradient);
	int	End				= this->NumGridPts;
	N_VScale(0.0, Gradient, Gradient);				//Zero out
	//Loop
	GradD[0]=(SD[0]-GD[0])/this->delx;				// Boundary
	for( int i = 1 ; i < End; i ++)
		GradD[i] = ( SD[i] - SD[i - 1])/this->delx;
}

//Called by: ???
//Status:  Most likely dangling
void myPb2::SetGradient(N_Vector Input, N_Vector GradientVec, realtype Ghost)
{
	//Declare
	realtype * SD		= N_VGetArrayPointer(Input);
	realtype * GD		= N_VGetArrayPointer(this->Ghost);
	realtype * GradD	= N_VGetArrayPointer(GradientVec);
	int	End				= this->NumGridPts;
	N_VScale(0.0, GradientVec, GradientVec);				//Zero out
	//Loop
	GradD[0]=(SD[0]-Ghost)/this->delx;				// Boundary
	for( int i = 1 ; i < End-1; i ++)
		GradD[i] = ( SD[i] - SD[i - 1])/this->delx;
}

//        d88P   d88P  .d8888b.                             888                              888                    
//       d88P   d88P  d88P  Y88b                            888                              888                    
//      d88P   d88P   888    888                            888                              888                    
//     d88P   d88P    888         .d88b.  88888b.  .d8888b  888888 888d888 888  888  .d8888b 888888 .d88b.  888d888 
//    d88P   d88P     888        d88""88b 888 "88b 88K      888    888P"   888  888 d88P"    888   d88""88b 888P"   
//   d88P   d88P      888    888 888  888 888  888 "Y8888b. 888    888     888  888 888      888   888  888 888     
//  d88P   d88P       Y88b  d88P Y88..88P 888  888      X88 Y88b.  888     Y88b 888 Y88b.    Y88b. Y88..88P 888     
// d88P   d88P         "Y8888P"   "Y88P"  888  888  88888P'  "Y888 888      "Y88888  "Y8888P  "Y888 "Y88P"  888   




//=============================================
//Verification tests for the Simpson Integrator
//=============================================
//Called by:  Main in testing mode.
void myPb2::RunTests(N_Vector State)
{
	//Solve with V(0)=1, totally arbitrary.
	N_VConst(1,this->Vel);
	N_Vector fAve 	= N_VNew_Serial(this->NumGridPts+1);
	realtype * FAD	= N_VGetArrayPointer(fAve);
	N_VConst(1,fAve);

	this->VelIntegrate( FAD, State, 1,1);
	std:: cout << "Updated Const Velocity" << std :: endl;
	N_VPrint_Serial(this->Vel);
	std :: cout << "End Const Experiment" << std :: endl;

	N_VConst(1,this->Vel);
	SinWave(FAD, this->NumGridPts, this->NumGridPts, this->delx);
	this->VelIntegrate(FAD, State, FAD[0], FAD[this->NumGridPts-1]);
	std :: cout << "Updated Cos Velocity" << std :: endl;
	N_VPrint_Serial(this->Vel);
	std :: cout << "End Cos Experiment" << std :: endl;


	SinWave2(FAD, this->NumGridPts, this->NumGridPts, this->delx);
	this->VelIntegrate(FAD, State, FAD[0], FAD[this->NumGridPts-1]);
	std :: cout << "Updated Sin Velocity" << std :: endl;
	N_VPrint_Serial(this->Vel);
	std :: cout << "End Sine Experiment" << std :: endl;
	N_VDestroy(fAve);
}

//===========================
//Additional Set up  Wrapper
//===========================
//Calls:   	this->SetGhost()
//			this->ScaleP()
//			this->SetVels()
//Called by:  Main during setup.
void myPb2::SetGhostPVel(N_Vector y, int Experiment, int SampleNum, realtype VelVal)
{
	this->SetGhost(y);
	this->ScaleP(Experiment, SampleNum);
	if(SampleNum ==10|| SampleNum == 12 )//AdvDiffCoS problem
		this->SetVels(this->NumGridPts+1, 1);
	else
		this->SetVels(this->NumGridPts+1,VelVal);
}

//Called by: ???
void myPb2::CheckNaN(N_Vector State, int vecLength)
{
	realtype * DATA = N_VGetArrayPointer(State);
	for( int i = 0 ; i < vecLength; i ++)
	{
		if(isnan(DATA[i]) || abs(DATA[i]) >1e5)
		{
			if( isnan(DATA[i]) )
				std :: cout << "Data has NaN: " <<  i << DATA[i] << "\n";
			if(abs(DATA[i]) > 1e5)
				std :: cout << "Data is out of control: "<< i << DATA[i] <<"\n";
			exit(EXIT_FAILURE);
		}
	}
}



//End Class functions
//===========================================================
//||  <:::::>  <::::>  <:::::>     <:::>         <::::>    ||
//||     .*/   |::     |::  ::>   <>   <>        |::  `>   ||
//||    .*/    |::::>  |:::::<   <>     <>  <::> |::   *>  ||
//||   .*/     |::     |:: ':>    <>   <>        |::  .>   ||
//||  <:::::>  <::::>  |::   :>    <:::>         <::::>    ||
//===========================================================
//Called by myPb to solve the single boundary point chemistry problem.
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
//Called by:  ????
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

/*
//=============================================================
//  __ __    ___   _____      __   __ ___   ___
// |  V  | ./ _ \.|_   _| ___ \ \ / /| __) / __)
// | .V. | | (_) |  | |  |___| \ v / | __)( (__
// |_| |_| |_/ \_|  |_|         \_/  |___) \___)
//============================================================
//We cover our bases by doing MatVec directly for EPI methods.
//By doing this, we don't have to transpose the data.
*/


void MatVecProdFast(int len, N_Vector Jac, N_Vector x, N_Vector tmp, realtype * JV)
{
	realtype * X 			= NV_DATA_S(x);
	realtype * JacD			= NV_DATA_S(Jac);
	for (int i=0; i< len; i++)//for each state
	{
		for(int j=0; j<len; j++)//marches across the column
		{
			JV[i] += X[j] * JacD[j+ i * len];
		}
	}
}

void MatVecProdParallel(int len, N_Vector Jac, N_Vector x, N_Vector tmp, realtype * JV)
{
#pragma omp parallel num_threads(4)
{	
	//realtype * TMP  		= NV_DATA_S(tmp);
	realtype * X 			= NV_DATA_S(x);
	realtype * JacD			= NV_DATA_S(Jac);
	double* PrivateJV = (double*)calloc(len, sizeof(double));
	for (int i=0; i< len; i++)//for each state
	{
		for(int j=0; j<len; j++)//marches across the column
		{
			PrivateJV[i] += X[j] * JacD[j+ i * len];
		}
	}
	#pragma omp critical
    {
        for(int i =0; i<len; i++) 
			JV[i] += PrivateJV[i];
    }
	free(PrivateJV);
}
// #pragma omp parallel num_threads(4)
// {
//     int y,i;
//     double* results_private = (double*)calloc(matrix_size, sizeof(double));
//     for(y = 0; y < matrix_size ; y++) {
//         #pragma omp for
//         for(i = 0; i < matrix_size; i++) {
//             results_private[y] += vector[i]*matrix[i][y];   
//         }
//     }
//     #pragma omp critical
//     {
//         for(y=0; y<matrix_size; y++) results[y] += results_private[y];
//     }
//     free(results_private);
// }
}



//=================================================
//||      <:::>   :>   ::  <::::>      <::::>    ||
//||     <>   <>  |:>  ::  |::         |::  `>   ||
//||    <>     <> |: > ::  |::::> <::> |::   *>  ||
//||     <>   <>  |:  >::  |::         |::  .>   ||
//||      <:::>   |:   >:  <::::>      <::::>    ||
//=================================================
//===============================
// ______  _____  ______
//|_    _||  _  ||   ___|
// _|  |  |     ||  |___
//|____|  |__|__||______|
//================================
//Called directly by EPIP2 Marked to be removed
int SUPER_CHEM_JAC_TCHEM(realtype t, N_Vector State, N_Vector StateDot, SUNMatrix Jac, void * pb, N_Vector tmp1,
				N_Vector tmp2, N_Vector tmp3)
{
	myPb2 * pbPtr{static_cast<myPb2 *> (pb)};//Recast
	auto Start=std::chrono::high_resolution_clock::now();
	int Length = N_VGetLength(State);
	int num_eqs = pbPtr->num_equations;
	int grid_sz = pbPtr->NumGridPts;	//Length/num_eqs; //This give a grid size/ # copies
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
		for( int j = 0 ; j < num_eqs; j ++)
			for(int k = 0; k < num_eqs; k++)
				JACDATA[i*Block + j * num_eqs + k] = TEMPJACDATA[j*num_eqs + k];
	}
	auto Stop=std::chrono::high_resolution_clock::now();
	auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
	pbPtr->jacMakeTime+=Pass.count()/1e9;
	N_VDestroy_Serial(yTemp);
	N_VDestroy_Serial(JacTemp);
	return 0;
}



//=================
//MOVERS
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
		SUPERDATA[i + j * grid_sz] = VECDATA[j];
		//SUPERDATA[j + i * grid_sz] = VECDATA[j];
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

//==============================
//Needed to get the inlet ghost.
//==============================
void myPb2::GhostChem(realtype t, N_Vector y, N_Vector ydot, void * pb)
{
	myPb2 * problem{static_cast<myPb2 *> (pb)};
	N_VScale(0.0, problem->SmallChem, problem->SmallChem);
	CHEM_RHS_TCHEM(t, y, ydot, pb);		//uses the zerod problem
}


//Will be removed in the final version
//================================
//Check heating
//================================
void myPb2::CheckHeating(N_Vector State, realtype t)
{
	if(this->HeatingOn ==1)
		if(t>1e-4)
		{
		//if(N_VMaxNorm(State)> 2600 && t > 1e-6)//2600 originally//usually 2100
			this->HeatingOn = 0;
			std :: cout <<"x";
		}

}

//Will be removed in final version
void myPb2::VerifyHeatingExp(N_Vector State, N_Vector State0, realtype tElapsed)
{
	if( this->Adv == 0.0 && this->Diff == 0.0 && this->React == 0.0 && this->Power!=0)
	{
		//N_VPrint_Serial(State0);
		CheckHeating(State0, 0);
		if(this->HeatingOn==1)
			std:: cout << "Verifying Heating\n";
		else
		{
			std :: cout << "Error with the Heating Check\n";
		 	exit(EXIT_FAILURE);
		}
	N_Vector Powered 	= N_VClone(State0);
	N_VScale(1.0, State0, Powered);
	this->RHS_Heat(tElapsed, State0, Powered, this);
	//SUPER_RHS_HEATING(tElapsed, State0, Powered, this);
	N_VScale(tElapsed*this->Power, Powered, Powered);		//RHS*Time
	N_VLinearSum(1.0, State0, 1.0, Powered, Powered);		//S(0)+RHS*Time
	N_VLinearSum(1.0, State, -1.0, Powered, Powered);
	std :: cout << "Powered Error: " << sqrt(N_VDotProd(Powered,Powered)) << std :: endl ;
	}
	else
	{
		std :: cout << "Heating check skipped\n";
	}


}
//==============================================
//  _       ___     ___   _    _  _   _   ____ 
// | |     /   \   /   \ | |_/ / | | | | |    \
// | |    |     | |     ||   /   | | | | |  __/
// | |___ |     | |     ||  _ \  | |_| | | |
// |_____| \___/   \___/ |_| \_\  \___/  |_|
//===============================================

//To be removed from release version
void myPb2::VerifyTempTable(N_Vector State)
{
	realtype * Data 		= N_VGetArrayPointer(State);
	realtype * LookupTemp	= N_VGetArrayPointer(this->TempTable);
	realtype * LookupCp		= N_VGetArrayPointer(this->CPPoly);
	realtype * LookupDiff	= N_VGetArrayPointer(this->DiffTable);
	realtype * LookupRho	= N_VGetArrayPointer(this->RhoTable);
	realtype Temp			= Data[0];			//Just check against temp
	std :: cout << "Temperature: " << Temp << "\n";
	int TempInd				= 0;
	realtype indErr			= 5;
	int lookupLen			= N_VGetLength(this->TempTable);
	for(int i = 0; i < lookupLen; i++)
	{
		if(abs(Temp-LookupTemp[i])< indErr)
		{
			indErr=abs(Temp-LookupTemp[i]);
			TempInd=i;
		}
	}
	if( TempInd == this->TempTableLookUp(Temp, this->TempTable) )
	std :: cout << "Table lookup functions agree\n";
	std :: cout << "Temp Lookup index: " << TempInd << "\n" << "Extrapolated cp: " << LookupCp[TempInd];
	std :: cout << "Extrapolated density: " << LookupRho[TempInd] << "\n";
	for( int i = 0 ; i < this->num_equations; i ++)
	{
		std :: cout << "Extrapolated D_" <<i << "\t" << LookupDiff[ TempInd + i*500] << "\n";
	}
}


//================================================================================
//Sets the Transport properties based off the table provided to us by Dr. Bisetti
//Pushed to production
//================================================================================
int myPb2::TempTableLookUp(realtype Temp, N_Vector TempTable)
{
	realtype * LookupTemp	= N_VGetArrayPointer(TempTable);
	int len					= N_VGetLength(TempTable);
	int TempInd				= 0;
	realtype indErr			= 5;
	for(int i = 0 ; i < len ; i ++)
	{
		if(abs(Temp-LookupTemp[i])< indErr)
		{
			indErr=abs(Temp-LookupTemp[i]);
			TempInd=i;
		}
	}
	return TempInd;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                                                                                          
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                                                                                          
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                                                                                          
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                                                                                          
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                                                                                          
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                                                                                          
//  ad88888ba                           88                                                   
// d8"     "8b                          88                                                   
// Y8,                                  88                                                   
// `Y8aaaaa,     ,adPPYba,  ,adPPYYba,  88  ,adPPYYba,  8b,dPPYba,                           
//   `"""""8b,  a8"     ""  ""     `Y8  88  ""     `Y8  88P'   "Y8                           
//         `8b  8b          ,adPPPPP88  88  ,adPPPPP88  88                                   
// Y8a     a8P  "8a,   ,aa  88,    ,88  88  88,    ,88  88                                   
//  "Y88888P"    `"Ybbd8"'  `"8bbdP"Y8  88  `"8bbdP"Y8  88                                                                                                                            
//   ,ad8888ba,                                    88  88                                    
//  d8"'    `"8b                                   88  ""                             ,d     
// d8'                                             88                                 88     
// 88             8b,dPPYba,  ,adPPYYba,   ,adPPYb,88  88   ,adPPYba,  8b,dPPYba,   MM88MMM  
// 88      88888  88P'   "Y8  ""     `Y8  a8"    `Y88  88  a8P_____88  88P'   `"8a    88     
// Y8,        88  88          ,adPPPPP88  8b       88  88  8PP"""""""  88       88    88     
//  Y8a.    .a88  88          88,    ,88  "8a,   ,d88  88  "8b,   ,aa  88       88    88,    
//   `"Y88888P"   88          `"8bbdP"Y8   `"8bbdP"Y8  88   `"Ybbd8"'  88       88    "Y888  
// 88b           d88                        88               88                              
// 888b         d888                        88               88                              
// 88`8b       d8'88                        88               88                              
// 88 `8b     d8' 88   ,adPPYba,    ,adPPYb,88  88       88  88   ,adPPYba,                  
// 88  `8b   d8'  88  a8"     "8a  a8"    `Y88  88       88  88  a8P_____88                  
// 88   `8b d8'   88  8b       d8  8b       88  88       88  88  8PP"""""""                  
// 88    `888'    88  "8a,   ,a8"  "8a,   ,d88  "8a,   ,a88  88  "8b,   ,aa                  
// 88     `8'     88   `"YbbdP"'    `"8bbdP"Y8   `"YbbdP'Y8  88   `"Ybbd8"'                                                                                
//==================================================================================================
//==================================================================================================
// 
//  _____             _                                         _        _              _    ___    ___      __  ___    ___      __ ____   ____  
// /__   \  ___  ___ | |_  ___    ___   ___   _ __ ___   _ __  | |  ___ | |_   ___   __| |  / _ \  ( _ )    / / / _ \  / _ \    / /|___ \ |___ \ 
//   / /\/ / _ \/ __|| __|/ __|  / __| / _ \ | '_ ` _ \ | '_ \ | | / _ \| __| / _ \ / _` | | | | | / _ \   / / | | | || (_) |  / /   __) |  __) |
//  / /   |  __/\__ \| |_ \__ \ | (__ | (_) || | | | | || |_) || ||  __/| |_ |  __/| (_| | | |_| || (_) | / /  | |_| | \__, | / /   / __/  / __/ 
//  \/     \___||___/ \__||___/  \___| \___/ |_| |_| |_|| .__/ |_| \___| \__| \___| \__,_|  \___/  \___/ /_/    \___/    /_/ /_/   |_____||_____|
//                                                      |_|                                                                                                                  
//The test call is:
//./OneDIgnCons.x --chemfile=gri3.0/chem.inp --thermfile=gri3.0/therm.dat --StepSize=2e-7 --FinalTime=4e-7 --MyFile="Scrap.txt"  --KrylovTol=1e-7 \
//--UseJac=1 --SampleNum=4 --Method="EPI2" --Experiment=2 --NumPts=100 --Pow=1.0 --VelUp=1 --Profiling=1 --Movie=0 --absTol=1e-10 --relTol=1e-8 --Delx=1e-2

void myPb2::Set_ScalarGradient(N_Vector State)
{
	realtype * uData        = NV_DATA_S(State);//Stuff comes from here.
	realtype * resultData	= NV_DATA_S(this->ScalarGradient);
    int numPts 				= this->NumGridPts;
    int grid 				= 0;
	int TI					= 0;
    realtype * Ghost 		= NV_DATA_S(this->Ghost);
	realtype x				= 0;
	realtype denom			= 1.0/(2 * this->delx);
    int vecLength 			= this->num_equations * this->NumGridPts;
	//std :: cout << this->delx << std :: endl;
	//Start main loop
	for (int i = 0; i < vecLength; i++)
	{
		TI		= i % numPts;
		grid 	= floor(i/this->NumGridPts);
		x 		= TI*this->delx + this->delx/2.0;
		if(TI==0)//Left
		{
			//std :: cout << "Ghost: " <<Ghost[grid];
			resultData[i]		= (uData[i+1] - Ghost[grid])* denom;
			//std :: cout << " x: " << x <<" Boundary derivative: " << resultData[i] << "\t\t Data @ x" << uData[i] << "\n\n";  //<< "\t\t Data @ i+1: " << uData[i+1] <<std :: endl;
		}
		else if(TI==this->NumGridPts-1)//Right
		{
			resultData[i]		= (uData[i] - uData[i-1]) * denom;
			//std :: cout << "x: " << x <<" Derivative: " << resultData[i] << "\t\t Data @ x " << uData[i] << "\n\n"; //"\t\t Data @ i-1: " << uData[i-1] << "\n\n";
		}
		else//middle
		{
			resultData[i]		= (uData[i+1] - uData[i-1]) * denom;
			//std :: cout << "x: " << x <<" Derivative: " << resultData[i] << "\t\t Data @ x " << uData[i] << "\n\n"; //"\t\t Data @ i-1: " << uData[i-1] << "\t\t Data @ i+1: " << uData[i+1] << std :: endl;
		}
    }
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//This derivative is on the scalar grid
//Should be called only once per step.
//Set Gradients of the transport properties.  Called by cross diffusion.
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void myPb2::Set_TransportGradient(N_Vector State)
{
	//realtype * uData        = NV_DATA_S(State);//Stuff comes from here.
	realtype * resultData	= NV_DATA_S(this->ScalarGradient); //Target for result
	int index				= 0;
	realtype 	denom		= 1.0/(2*this->delx);
	//Used for ghost data
	realtype * LookupTemp	= NV_DATA_S(this->TempTable);
	realtype * LookupCp		= N_VGetArrayPointer(this->CPPoly);
	realtype * LookupDiff	= N_VGetArrayPointer(this->DiffTable);
	realtype * LookupRho	= N_VGetArrayPointer(this->RhoTable);
	//Grids resolved for that spots current temperature (user need make sure it is updated).
	realtype * CpData		= NV_DATA_S(this->CpGrid);
	realtype * RhoData		= NV_DATA_S(this->RhoGrid);
	realtype * DiffData 	= NV_DATA_S(this->DiffGrid);
	realtype * data 		= NV_DATA_S(State);
	realtype * GhostP		= NV_DATA_S(this->Ghost);
	//Get the transport scalars derivative pointers.
	index				= this->TempTableLookUp(GhostP[0], this->TempTable);//Ghost temp
	//std :: cout << "o\n";
	realtype * GradDiff	= NV_DATA_S(this->DiffGrad);  //This line randomly cannot access DiffGrad when gridpts <50
	//std :: cout << "x\n";

	realtype * GradT	= NV_DATA_S(this->TempGrad);
	realtype * GradRho	= NV_DATA_S(this->RhoGrad);
	realtype * GradCp	= NV_DATA_S(this->CpGrad);


	//Run the boundary data
	//Left boundary
	GradRho	[0]	=	(RhoData[1] - LookupRho[index])	*denom; // Boundary GradRho
	GradCp	[0] = 	(CpData[1]- LookupCp[index])	*denom;	// Boundary Cp
	//Right boundary
	int End= this->NumGridPts-1;
	GradRho	[End]	=	(RhoData[End] - RhoData[End-1])	*denom; // Boundary GradRho
	GradCp	[End] = 	(CpData[End]- CpData[End-1])	*denom;	// Boundary Cp

	//std:: cout << CpData[1] << " " << LookupCp[index] << std :: endl;
	
	//This needs to be fixed. Does the full set of boundary 
	for(int i = 0 ; i < this->num_equations; i ++)
	{	//============================================v jth diffusion @ i+1              v Ghost of jth diffusion
		GradDiff[this->NumGridPts*i] 		=	( DiffData[i*this->NumGridPts + 1] - LookupDiff[i*500 + index])*denom; //Left
		//============================================v jth diffusion @ end v jth diffusion @ end-1	
		GradDiff[this->NumGridPts*(i+1)-1] 	=	( DiffData[(i+1)*this->NumGridPts-1]- DiffData[(i+1)*this->NumGridPts-2] )*denom; //Right
		
		// GradDiff[this->num_equations*(i+1)-1] 	=	( DiffData[(i+1)*500-1]- DiffData[(i+1)*500-2] )*denom; //Right
	}

	//Do the large nasty grid.
	int gridjump 		= this->NumGridPts;
	//Loop over interior
	for( int i = 1 ; i < this->NumGridPts-1; i ++)
	{
		GradRho[i]	= ( RhoData[i+1] - RhoData[i-1])/(2*this->delx);
		GradCp[i]	= ( CpData[i+1] - CpData[i-1])/(2*this->delx);
		for( int j = 0; j < this->num_equations ; j ++)
		{
			GradDiff[j*gridjump + i] = (DiffData[j* gridjump+ i +1] - DiffData[j*gridjump+ i - 1 ])*denom;
		}
	}

}

//---------------------------------------------------------------------------------------------------------                                                                                          
//---------------------------------------------------------------------------------------------------------                                                                                          
//---------------------------------------------------------------------------------------------------------                 
//Gas weight is 1/(1/W) when returned
void myPb2::Set_GasWeightBisetti(N_Vector State)
{
	realtype * GWPtr 	= NV_DATA_S(this->GasWeight);
	realtype * SDPtr 	= NV_DATA_S(State);
	realtype * SMPtr 	= NV_DATA_S(this->MolarWeights);
	N_VConst(0.0, this->GasWeight);
	//March through each grid point and grab each mass fraction and multiply by the weight of each molecule
	for(int i = 0; i < this->NumGridPts; i ++)
	{
		for(int j = 1; j < this->num_equations; j++)
		{	//GW = 1/ Sigma Y_i/W_i

			GWPtr[i]	+= SDPtr[i + j * this->NumGridPts]/ SMPtr[j];
			// std :: cout << SDPtr[i + j * this->NumGridPts] << std :: endl;
			// std :: cout << SMPtr[j] << std :: endl;
			
		}
		//Then invert the entry.
		GWPtr[i] = 1/GWPtr[i];
		//std :: cout << GWPtr[i] << std :: endl;
		//Outputs uniform gas weight.
	}


	//N_VPrint_Serial(this->GasWeight);  //I don't understand why this line prints 0's

}

//88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
//```````````````````````````````````````````````````````````````````````````````````````````````````````````
// oooooo     oooo           oooo                       o8o      .                                         
//  `888.     .8'            `888                       `"'    .o8                                         
//   `888.   .8'    .ooooo.   888   .ooooo.   .ooooo.  oooo  .o888oo oooo    ooo                           
//    `888. .8'    d88' `88b  888  d88' `88b d88' `"Y8 `888    888    `88.  .8'                            
//     `888.8'     888ooo888  888  888   888 888        888    888     `88..8'                             
//      `888'      888    .o  888  888   888 888   .o8  888    888 .    `888'                              
//       `8'       `Y8bod8P' o888o `Y8bod8P' `Y8bod8P' o888o   "888"     .8'                               
//                                                                   .o..P'                                
//                                                                   `Y8P'                                 
                                                                                                        
// oooooooooo.    o8o                                                                                      
// `888'   `Y8b   `"'                                                                                      
//  888      888 oooo  oooo    ooo  .ooooo.  oooo d8b  .oooooooo  .ooooo.  ooo. .oo.    .ooooo.   .ooooo.  
//  888      888 `888   `88.  .8'  d88' `88b `888""8P 888' `88b  d88' `88b `888P"Y88b  d88' `"Y8 d88' `88b 
//  888      888  888    `88..8'   888ooo888  888     888   888  888ooo888  888   888  888       888ooo888 
//  888     d88'  888     `888'    888    .o  888     `88bod8P'  888    .o  888   888  888   .o8 888    .o 
// o888bood8P'   o888o     `8'     `Y8bod8P' d888b    `8oooooo.  `Y8bod8P' o888o o888o `Y8bod8P' `Y8bod8P' 
//                                                    d"     YD                                            
//                                                    "Y88888P'                                            
                                                                                                        
// oooooooooooo                                       .    o8o                                             
// `888'     `8                                     .o8    `"'                                             
//  888         oooo  oooo  ooo. .oo.    .ooooo.  .o888oo oooo   .ooooo.  ooo. .oo.                        
//  888oooo8    `888  `888  `888P"Y88b  d88' `"Y8   888   `888  d88' `88b `888P"Y88b                       
//  888    "     888   888   888   888  888         888    888  888   888  888   888                       
//  888          888   888   888   888  888   .o8   888 .  888  888   888  888   888                       
// o888o         `V88V"V8P' o888o o888o `Y8bod8P'   "888" o888o `Y8bod8P' o888o o888o                      
//__________________________________________________________________________________________________________
//8888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888                                                                                                        
//Checking part by part

//=====================================================================
//Verified on both Sine and Cosine waves to show second order in space.
//=====================================================================
//Velocity integrate function
//Solves:	V(x) = integral( f(x), [0,x] ) + V(0);
//Pass the off centered values as an N_Vector and it will calculate a second order integral.
//Called by:  Main, during the integration loop after Set_VelocityDivergence
//==============================================
void myPb2::VelIntegrate(realtype * fMid, N_Vector State, realtype LeftGhost, realtype RightGhost)
{
	realtype * VelData	= N_VGetArrayPointer(this->Vel);
	realtype * ScrapData	= N_VGetArrayPointer(this->SmallScrap);
	realtype * VelIntData	= N_VGetArrayPointer(this->VelScrap);
	int End 		= this->NumGridPts;
	//Clean up Data
	N_VScale(0.0, this->VelScrap, this->VelScrap);
	N_VScale(0.0, this->SmallScrap, this->SmallScrap);

	//Integrate each block's approximate area given we know only x's data 'exactly'
	//|-x-|-x-|-x-|-x-|-x-| x's are staggered off main grid, and |'s are on the grid.
	for(int i = 1 ; i < End-1; i++)
		ScrapData[i] 	= this->delx/ 6 * ( (fMid[i+1] + fMid[i-1])/2 + 5 * fMid[i] );
	//Set Boundaries
	ScrapData[0]	=	this->delx/6 * ( (fMid[1] + LeftGhost)/2 + 5*fMid[0]);
	ScrapData[End-1]=	this->delx/6 * ( (RightGhost + fMid[End-2])/2 + 5*fMid[End-1] );

	//Perform Box summation for: integral from 0 to x of (fAve(x))dx
	for(int i = 0; i < End; i ++)
		for(int j = 0; j <= i ; j ++)
			VelIntData[i+1]+=	(ScrapData[j]);

	//Finalize
	N_VAddConst(this->VelScrap, VelData[0], this->Vel);//New
}

void myPb2::Set_VelocityDivergence(N_Vector State)
{
	//ToDo: Set Global Boundary Left and right

	//Top Declarations==================================
	int Length 			= N_VGetLength(State);
	realtype * TmpData 	= NV_DATA_S(this->Tmp);
	realtype * GhostD	= NV_DATA_S(this->Ghost);
	//This houses the velocity RHS boundary data
	realtype * VelBndPtr= NV_DATA_S(this->SmallScrap);
	realtype * CpPtr	= NV_DATA_S(this->CpGrid);
	N_VScale(0.0, this->Tmp, this->Tmp);
	N_VScale(0.0, this->SmallScrap, this->SmallScrap);
	realtype GasConst  	= this->GasConst;  //8.31446261815324;
	this->Set_GasWeightBisetti(State);

	realtype * SD			= N_VGetArrayPointer(State);
	realtype * TmpD			= N_VGetArrayPointer(Tmp);			//Small Tmp data
	realtype * SCP			= N_VGetArrayPointer(this->SmallChem);
	realtype * GD			= N_VGetArrayPointer(this->Ghost);
	realtype * TG			= N_VGetArrayPointer(TempGrad);
	//Lookup data tables pointers

	realtype LeftTemp		= 0;
	realtype RightTemp		= 0;
	int End 				= this->NumGridPts-1;			//The final vector entry
	//See Dr. Bisetti's book, chapter 18 for full details.

	//===================================================
	//Given the temperature, set all the grids
	this->SetTransportGrid(State);
	this->Set_ScalarGradient(State);
	this->Set_TransportGradient(State);
	//Species Terms==========================================
	//WARNING:  This must be run first, writes directly to the VelAve
	//=======================================================
	this->Set_VelocityDivergence_SpecDiff(State);
	//N_VPrint_Serial(this->VelAve);
	
	//Temp Terms============================================
	this->Set_VelocityDivergence_TempReac(State);				//Data Written to small scrap
	N_VLinearSum(1.0, this->SmallScrap, 1.0, this->VelAve, this->VelAve);
	//Set the boundary
	SUPER_2_VEC(End, NV_DATA_S(this->PointScrap), SD, this->num_equations, this->NumGridPts);//copy right ghost state from far right

	this->GhostChem(0, this->PointScrap, this->SmallChem, this);			//Right Floating Chem using Scratch
	
	this->VelAveLeftBnd= NV_DATA_S(SmallChem)[0]/GD[0];


	//N_VPrint_Serial(this->VelAve);
	this->Set_VelocityDivergence_TempDiff(State);				//Data written to small scrap
	N_VLinearSum(1.0, this->SmallScrap, 1.0, this->VelAve, this->VelAve);
	//N_VPrint_Serial(this->VelAve);
	//Set the boundary
	this->VelAveRightBnd = VelBndPtr[End];

	//Heating needs to be treated differently
	this->Set_VelocityDivergence_TempHeat(State);
	N_VLinearSum(1.0, this->SmallScrap, 1.0, this->VelAve, this->VelAve);
	//Set boundary
	if(this->HeatingOn==1)
	{
		this->RHS_Heat(0, State, this->Tmp, this);
		this->HeatingRightGhost	= this->Power*this->HeatingRightGhost;	//Scale on/off the ghost point
		this->VelAveRightBnd+= this->HeatingRightGhost/SD[End];			//Call the end boundary.
	}
	//N_VPrint_Serial(this->VelAve);

	//Finalize==============================================
	N_VScale(0.0, this->Tmp, this->Tmp);							//Clean Tmp		
}

//====================================================
//  #     #                                           
//  #     # ###### #      #####  ###### #####   ####  
//  #     # #      #      #    # #      #    # #      
//  ####### #####  #      #    # #####  #    #  ####  
//  #     # #      #      #####  #      #####       # 
//  #     # #      #      #      #      #   #  #    # 
//  #     # ###### ###### #      ###### #    #  ####  
//=========================================================
//                                                                                 
//Temperature Diffusion and Cross Diffusion
// D_T/T nabla^2 T + 1/(rho Cp T) nabla(lambda) nabla(T)
void myPb2::Set_VelocityDivergence_TempDiff(N_Vector State)
{
	int Length 			= N_VGetLength(State);
	realtype * SData 	= NV_DATA_S(State);
	realtype * TmpData 	= NV_DATA_S(this->Tmp);
	realtype * GhostD	= NV_DATA_S(this->Ghost);
	//This houses the velocity RHS boundary data
	realtype * VelBndPtr= NV_DATA_S(this->SmallScrap);
	realtype * CpPtr	= NV_DATA_S(this->CpGrid);
	realtype * RhoPtr 	= NV_DATA_S(this->RhoGrid);
	realtype * GasWPtr  = NV_DATA_S(this->GasWeight);
	N_VScale(0.0, this->Tmp, this->Tmp);
	N_VScale(0.0, this->SmallScrap, this->SmallScrap);

	//Get D_T/T * nabla^2 T
	this->RHS_Diff(t, State, this->Tmp, this);
	for( int i = 0; i < this->NumGridPts; i ++)
	{
		VelBndPtr[i]	= TmpData[i]/SData[i];
	}
	//Start cross diffusion
	this->RHS_CrossDiff(t, State, this->Tmp, this);
	//Take the data and modify the Temperature diffusion by R/cp (dimensionless quantity).
	//p = rho R T
	//So R/cp = p/(rho * cp * T)
	//[Grad lambda dot Grad T] = M/(T^3 THETA) * THETA/ L = M/(T^3 L)
	// note [p] = M/(L T^2), so dividing by pressure gives 1/T
	// So the prefactor looks like : 1/ (rho R T ) * R/ cp = 1/(rho cp T )
	//But some terms are already included in RHS_CrossDiff, namely 1/(rho cp)
	for( int i = 0; i < this->NumGridPts; i ++)
	{
		VelBndPtr[i]	+= TmpData[i]/SData[i];	
	}
}

                                                                                       
//Reactive term for Temperature:  omegadot/( cp rho T)
void myPb2::Set_VelocityDivergence_TempReac(N_Vector State)
{
	realtype * CpPtr	= NV_DATA_S(this->CpGrid);
	realtype * RhoPtr	= NV_DATA_S(this->RhoGrid);
	realtype * SSPtr 	= NV_DATA_S(this->SmallScrap);
	realtype * TmpPtr 	= NV_DATA_S(this->Tmp);
	realtype * Data 	= NV_DATA_S(State);

	//Generate the chemistry state data
	this->RHS_Chem(t, State, this->Tmp, this);
	
	// This is wrong, try the replacement above ln 2 in for loop
	// R/cp * 1/p where p = rho R T gives:  R/ ( cp * rho * R * T) = 1/(cp rho T)

	for (int i = 0 ; i < this->NumGridPts; i ++ )
	{
		SSPtr[i]		= TmpPtr[i] / Data[i];

		//SSPtr[i]		= TmpPtr[i] * (RhoPtr[i] * this->GasConst )/(this->pb._p* CpPtr[i]);
	}

}

//
//Generates the heating Gaussian term
void myPb2::Set_VelocityDivergence_TempHeat(N_Vector State)
{
	realtype * VTempData 	= NV_DATA_S(this->SmallScrap);
	realtype * TmpD			= NV_DATA_S(this->Tmp);
	realtype * SD			= NV_DATA_S(State);
	if(this->HeatingOn==1)
	{
		this->RHS_Heat(0, State, this->Tmp, this);
		N_VScale(this->Power, this->Tmp, this->Tmp);			//Scale on/off
		for( int i =0 ; i < this->NumGridPts; i ++)
			VTempData[i] 	= TmpD[i]/SD[i];		
	}
	//N_VPrint_Serial(this->SmallScrap);
}

//==============================
//Generates the species term
//==============================
void myPb2::Set_VelocityDivergence_SpecDiff(N_Vector State)//This MUST be called first
{
	//Will write to small scrap
	//Intermediary vectors are large
	//These vectors are large, and will be picked apart later.
	//Get the Chemistry term
	this->RHS_Chem(this->t, State, this->Tmp, this); //Saves into Tmp
	//Get Diffusion
	this->RHS_Diff(this->t, State, this->Scrap, this);//into Scrap
	//Combine into Scrap
	N_VLinearSum(1.0, this->Tmp, 1.0, this->Scrap, this->Scrap);
	//Now Chem & Diffusion are on the same grid.
	//Get Cross Diffusion onto Tmp
	this->RHS_CrossDiff(this->t, State, this->Tmp, this);
	//Combine
	N_VLinearSum(1.0, this->Scrap, 1.0, this->Tmp, this->Scrap);
	//All terms combined
	this->Set_GasWeightBisetti(State);
	//Get Pointer to Diff-Chem terms
	realtype * SpecDiffPtr 	= NV_DATA_S(this->Scrap);

	//Pointer to the molar weight
	realtype * MassPtr 		= NV_DATA_S(this->MolarWeights); //Input
	//This is the gas weight
	realtype * GasWPtr 		= NV_DATA_S(GasWeight);
	//Pointer to output N_Vector
	realtype * TmpPtr		= NV_DATA_S(this->VelAve);//Output
	N_VConst(0.0, this->VelAve);

	int 		StateInd	= 0;
	for(int i = 0; i < this->NumGridPts; i ++)
	{
		for(int j = 1 ; j < this->num_equations; j++)
		{	
			//We jump to each SPECIES diffusion @ point i in the grid
			StateInd 	=   i + j * this->NumGridPts;
			//output   add    weight @ i  spec j mass  * Diffusion @ i on grid j 
			TmpPtr[i]  += SpecDiffPtr[StateInd] /MassPtr[j];
		}
		TmpPtr[i] = GasWPtr[i] * TmpPtr[i];
	}

}




//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//                                                                                          
//  @@@@@@@@  @@@@@@@    @@@@@@   @@@  @@@  @@@@@@@@  @@@ @@@   @@@@@@   @@@@@@@   @@@@@@@   
// @@@@@@@@@  @@@@@@@@  @@@@@@@@  @@@  @@@  @@@@@@@@  @@@ @@@  @@@@@@@@  @@@@@@@@  @@@@@@@@  
// !@@        @@!  @@@  @@!  @@@  @@!  @@@  @@!       @@! !@@  @@!  @@@  @@!  @@@  @@!  @@@  
// !@!        !@!  @!@  !@!  @!@  !@!  @!@  !@!       !@! @!!  !@!  @!@  !@!  @!@  !@!  @!@  
// !@! @!@!@  @!@!!@!   @!@!@!@!  @!@  !@!  @!!!:!     !@!@!   @!@!@!@!  @!@!!@!   @!@  !@!  
// !!! !!@!!  !!@!@!    !!!@!!!!  !@!  !!!  !!!!!:      @!!!   !!!@!!!!  !!@!@!    !@!  !!!  
// :!!   !!:  !!: :!!   !!:  !!!  :!:  !!:  !!:         !!:    !!:  !!!  !!: :!!   !!:  !!!  
// :!:   !::  :!:  !:!  :!:  !:!   ::!!:!   :!:         :!:    :!:  !:!  :!:  !:!  :!:  !:!  
//  ::: ::::  ::   :::  ::   :::    ::::     :: ::::     ::    ::   :::  ::   :::   :::: ::  
//  :: :: :    :   : :   :   : :     :      : :: ::      :      :   : :   :   : :  :: :  :   
//
//  Marked for removal upon release
//=======================================================================
//=======================================================================                                                                       
// 888888888888                             88                            
//      88                           ,d     ""                            
//      88                           88                                   
//      88   ,adPPYba,  ,adPPYba,  MM88MMM  88  8b,dPPYba,    ,adPPYb,d8  
//      88  a8P_____88  I8[    ""    88     88  88P'   `"8a  a8"    `Y88  
//      88  8PP"""""""   `"Y8ba,     88     88  88       88  8b       88  
//      88  "8b,   ,aa  aa    ]8I    88,    88  88       88  "8a,   ,d88  
//      88   `"Ybbd8"'  `"YbbdP"'    "Y888  88  88       88   `"YbbdP"Y8  
// 88888888ba   88                           88               aa,    ,88  
// 88      "8b  88                           88                "Y8bbdP"   
// 88      ,8P  88                           88                           
// 88aaaaaa8P'  88   ,adPPYba,    ,adPPYba,  88   ,d8                     
// 88""""""8b,  88  a8"     "8a  a8"     ""  88 ,a8"                      
// 88      `8b  88  8b       d8  8b          8888[                        
// 88      a8P  88  "8a,   ,a8"  "8a,   ,aa  88`"Yba,                     
// 88888888P"   88   `"YbbdP"'    `"Ybbd8"'  88   `Y8a        
//Test_VelocityDivergence
//=========================================================================
void myPb2::Test_VelocityDivergence(N_Vector TestState)
{
	//The function writes to this->SmallScrap
	realtype * SData 		= NV_DATA_S(TestState);
	realtype * CpPtr		= NV_DATA_S(this->CpGrid);
	realtype * RhoPtr		= NV_DATA_S(this->RhoGrid);
	realtype * DiffPtr		= NV_DATA_S(this->DiffGrid);
	realtype Err			= 0;
	std :: cout << "Testing the variable Velocity Divergence\n";
	N_VScale(0.0, this->Tmp, this->Tmp);
	N_VScale(0.0, this->SmallScrap, this->SmallScrap);
	//Run a simple state check and see if it returns a 0 vector
	this->SetTransportGrid(TestState);
	this->Set_ScalarGradient(TestState);
	this->Set_TransportGradient(TestState);
	
	//Run the test of the TempDiff term, see if 0 vector comes out.
	this->Set_VelocityDivergence_TempDiff(TestState);
	if(N_VDotProd(this->SmallScrap, this->SmallScrap) < 10*this->delx )
	{
		std :: cout << "The velocity temp diffusion Div passed the Const test!\n";
	}
	std :: cout << N_VDotProd(this->SmallScrap, this->SmallScrap) << "\n\n";

	
	std :: cout << "Testing the Species Term: \n\n";
	//N_VPrint_Serial(this->MolarWeights);
	N_VConst(0.0, this->VelAve);
	this->Set_VelocityDivergence_SpecDiff(TestState);
	if( isnan(N_VDotProd(this->VelAve,this->VelAve)))
	{
		std :: cout << "Species term generated NaN!\n";
	}
	else if(abs(sqrt(N_VDotProd(this->VelAve,this->VelAve))) < 10*this->delx)
	{
		std :: cout << "we passes the Species Diff-Reac test\n";
	}
	std :: cout << "running an integration check \n\n";
	//N_VPrint_Serial(this->VelAve);
	//Running and integration test to see if it runs.
	this->VelIntegrate(NV_DATA_S(this->VelAve), State, 0.0, 0.0);
	//This should create cx where c is  NV_DATA_S(this->Tmp)[0] )
	N_VConst(0.0, this->VelScrap);
	realtype c = 	NV_DATA_S(this->VelAve)[0];
	//std :: cout << "Integration scaling factor: " << c << std :: endl;
	//N_VPrint_Serial(this->Tmp);
	for( int i = 0 ; i < this->NumGridPts+1; i++)
	{
		NV_DATA_S(this->VelScrap)[i] = i*this->delx * c;
	}
	//N_VPrint_Serial(this->Vel);
	N_VLinearSum(1.0, this->Vel, -1.0, this->VelScrap, this->VelScrap);
	Err = N_VDotProd(this->VelScrap, this->VelScrap);
	if(sqrt(Err) < 10*this->delx*c)
	{
		std :: cout << "Error within tolerances: Passed test\n";
	}
	//N_VPrint_Serial(this->VelScrap);
	std :: cout << "Error from const integration for Species Diff Reac: " << sqrt(Err) << std :: endl;



	std :: cout << "Running a full test of Velocity update with a constant state\n";

	this->Set_VelocityDivergence(TestState);
	//this->VelIntegrate(NV_DATA_S(this->VelAve), State, 0.0, 0.0);
	this->VelIntegrate(NV_DATA_S(this->VelAve), State, this->VelAveLeftBnd, this->VelAveRightBnd);

	for( int i = 0 ; i < this->NumGridPts+1; i++)
	{
		NV_DATA_S(this->VelScrap)[i] = i*this->delx * c;
	}
	//N_VPrint_Serial(this->Vel);
	N_VLinearSum(1.0, this->Vel, -1.0, this->VelScrap, this->VelScrap);
	Err = N_VDotProd(this->VelScrap, this->VelScrap);
	if(sqrt(Err) < 10*this->delx*c)
	{
		std :: cout << "Error within tolerances: Passed test\n";
	}
	//N_VPrint_Serial(this->VelScrap);
	std :: cout << "Error from const integration for Full test: " << sqrt(Err) << std :: endl;
	N_VPrint_Serial(this->Vel);
	//Test T(x) = cos(2*pi*x), but manually decouple transport from T.
	//This should return the 4*pi^2 D_T
	SinWave(NV_DATA_S(TestState), this->NumGridPts, this->NumGridPts*this->num_equations, this->delx) ;
	N_VAddConst(TestState, -1.0, TestState);
	//Fix the ghosts for the derivative calcs
	realtype * GhostD	= NV_DATA_S(this->Ghost);
	GhostD[0] 			= SData[0];
	this->Set_VelocityDivergence_TempDiff(TestState);
	//Use VelAve as a proxy for the one's grid
	N_VConst(1.0, this->VelAve);
	//N_VPrint_Serial(this->SmallScrap);
	Err = abs(N_VDotProd(this->VelAve, this->SmallScrap) - this->NumGridPts * 4 * M_PI * M_PI* DiffPtr[0]/(RhoPtr[0] * CpPtr[0] ) );
	if(Err < 10 * this->delx)
	{
		std :: cout <<"Congrats, the Temp Diff passes the Cosine test!\n";
	}
	std :: cout << "Error: " << Err << std :: endl;

}
//  _____             _                                         _        _              _    ___    ___      __ _  _     __ ____   ____  
// /__   \  ___  ___ | |_  ___    ___   ___   _ __ ___   _ __  | |  ___ | |_   ___   __| |  / _ \  ( _ )    / // |/ |   / /|___ \ |___ \ 
//   / /\/ / _ \/ __|| __|/ __|  / __| / _ \ | '_ ` _ \ | '_ \ | | / _ \| __| / _ \ / _` | | | | | / _ \   / / | || |  / /   __) |  __) |
//  / /   |  __/\__ \| |_ \__ \ | (__ | (_) || | | | | || |_) || ||  __/| |_ |  __/| (_| | | |_| || (_) | / /  | || | / /   / __/  / __/ 
//  \/     \___||___/ \__||___/  \___| \___/ |_| |_| |_|| .__/ |_| \___| \__| \___| \__,_|  \___/  \___/ /_/   |_||_|/_/   |_____||_____|
//                                                      |_|                                                                              
//Diagnostic testing for the scalar gradient grid.
//Use the testing initial conditions to check if we can pass the const and sine tests
int myPb2::Test_TransportGradient(N_Vector TestState)
{
	//Const test
	//Fill the transport data and Ghost with const data
	std :: cout << "Running Const 1 test\n";
	realtype * TSD 	= NV_DATA_S(TestState);
	realtype * GP   = NV_DATA_S(this->Ghost);
	realtype * RSD  = NV_DATA_S(this->Scrap);
	int failFlag 	= 0;
	int failState	= 0;
	for( int i = 0 ; i < this->num_equations; i++)
	{
		GP[i]= 1;
	}
	N_VConst(1.0, TestState);
	//Set grids to 1
	N_VConst(1.0, this->CpGrid);
	N_VConst(1.0, this->RhoGrid);
	N_VConst(1.0, this->DiffGrid);
	//Change the tables to match the grid for ghost calculations
	N_VConst(1.0, this->CPPoly);
	N_VConst(1.0, this->RhoTable);
	N_VConst(1.0, this->DiffTable);
	this->Set_TransportGradient(TestState);
	if(N_VMaxNorm(this->CpGrad)+N_VMaxNorm(this->RhoGrad)+ N_VMaxNorm(this->DiffGrad)==0)
		std :: cout << "Const test passed!\n";
	else
	{
		std :: cout << "Const test failed!\n";
		std :: cout << N_VMaxNorm(this->CpGrad)+N_VMaxNorm(this->RhoGrad)+ N_VMaxNorm(this->DiffGrad) << std :: endl;
	}
	std :: cout << "completed Transport const test\n";
	

	//=================
	//Begin sine test.
	//=================
	//Cosine(2 pi x) + 1
	//=================
	std :: cout << "\n\n\n\n";
	std :: cout << "Running Cosine test\n";
	//============v omg this was annoying==========================================
	SinWave(NV_DATA_S(CpGrid), this->NumGridPts, this->NumGridPts, this->delx);
	N_VScale(1.0, CpGrid, RhoGrid);
	//Set the giant DiffTable
	SinWave(NV_DATA_S(DiffGrid), this->NumGridPts, this->NumGridPts*this->num_equations, this->delx);
	realtype * GhostSetter 	= NV_DATA_S(CpGrid);
	//std :: cout << "testing ghost point: " << GhostSetter[0] << std :: endl;
	for( int i = 0 ; i < this->num_equations; i++)
	{
		GP[i]= GhostSetter[0];//symmetry of cos
	}
	N_VConst(GP[0], this->CPPoly);
	N_VConst(GP[0], this->RhoTable);
	N_VConst(GP[0], this->DiffTable);
	//N_VPrint_Serial(CpGrid);
	
	this->Set_TransportGradient(TestState);


	//Now that the data has been set, measure the error.
	//Reference is pm 2*pi Other sine wave function. 
	//===================================================
	//Set a reference for a single grid of NumGridPts.
	//===================================================
	//Ref: -2 pi Sine(2 pi x)
	N_Vector Ref 	=	N_VClone(this->SmallScrap);
	realtype * RefPt=	NV_DATA_S(Ref);
	realtype * BRefPt=	NV_DATA_S(this->Scrap);
	SinWave2(NV_DATA_S(Ref), this->NumGridPts, this->NumGridPts, this->delx);
	SinWave2(NV_DATA_S(this->Scrap), this->NumGridPts, this->NumGridPts*this->num_equations, this->delx);
	N_VAddConst(Ref, -1.0, Ref);
	N_VAddConst(this->Scrap, -1.0, this->Scrap);
	//scale								vthis is our reference
	N_VScale(-2.0 * M_PI, Ref,Ref);
	N_VScale(-2.0 * M_PI, this->Scrap,this->Scrap); 

	//cosine cp
	N_VLinearSum(1.0, this->CpGrad, -1.0, Ref, this->SmallScrap );
	if(N_VDotProd(this->SmallScrap, this->SmallScrap)< 10*this->delx )
	{
		std ::cout << "Congrats: CpGrad passed the Cosine test!\n";
	}
	else
	{
		std :: cout <<"Error not within tolerance!\n";
		failFlag =	1;
		failState ++;
	}
	std :: cout << "2-norm error from solution: " << N_VDotProd(this->SmallScrap, this->SmallScrap) << std::endl;
	//cosine rho
	N_VLinearSum(1.0, this->RhoGrad, -1.0, Ref, this->SmallScrap );
	if(N_VDotProd(this->SmallScrap, this->SmallScrap)< 10*this->delx )
	{
		std ::cout << "Congrats: RhoGrad passed the Cosine test!\n";
	}
	else
	{
		std :: cout <<"Error not within tolerance!\n";
		failFlag =	1;
		failState ++;
	}
	std :: cout << "2-norm error from solution: " << N_VDotProd(this->SmallScrap, this->SmallScrap) << std::endl;

	
	//Diff
	N_VLinearSum(1.0, this->DiffGrad, -1.0, this->Scrap, this->Scrap );

	if(N_VDotProd(this->Scrap, this->Scrap)< 10*this->delx )
	{
		std ::cout << "Congrats: Diffusion coefficient gradients passed the Cosine test!\n";
	}
	else
	{
		std :: cout <<"Error not within tolerance!\n" << std :: endl;
		failFlag =	1;
		failState ++;
	}
	std :: cout << "2-norm error from solution: " << N_VDotProd(this->Scrap, this->Scrap) << std::endl;
	std :: cout << "completed Transport sine test\n";



	
	//Start Sine->Cosine test
	//=================
	//Sine(2 pi x) + 1
	//=================
	std :: cout << "\n\n\n\n";
	std :: cout << "Running Sine test\n";
	//============v omg this was annoying==========================================
	SinWave2(NV_DATA_S(CpGrid), this->NumGridPts, this->NumGridPts, this->delx);
	N_VScale(1.0, CpGrid, RhoGrid);

	//Set the giant DiffTable
	SinWave2(NV_DATA_S(DiffGrid), this->NumGridPts, this->NumGridPts*this->num_equations, this->delx);
	for( int i = 0 ; i < this->num_equations; i++)
	{
		GP[i]= 2-GhostSetter[0];//antisymmetry of sine
	}
	N_VConst(GP[0], this->CPPoly);
	N_VConst(GP[0], this->RhoTable);
	N_VConst(GP[0], this->DiffTable);

	//N_VPrint_Serial(CpGrid);
	//N_VPrint_Serial(this->DiffGrid);

	
	this->Set_TransportGradient(TestState);

	//manually fix the right boundary conditions double them to account for the error
	NV_DATA_S(this->CpGrad)[this->NumGridPts-1]= 2*NV_DATA_S(this->CpGrad)[this->NumGridPts-1];
	NV_DATA_S(this->RhoGrad)[this->NumGridPts-1]= 2*NV_DATA_S(this->RhoGrad)[this->NumGridPts-1];

	for (int i  = 0 ; i < this->NumGridPts; i++)
		NV_DATA_S(this->DiffGrad)[(i+1)*(this->NumGridPts)-1]= 2*NV_DATA_S(this->DiffGrad)[(i+1)*(this->NumGridPts)-1];

	//===================================================
	//Set a reference for a single grid of NumGridPts.
	//===================================================
	//Ref: 2 pi Cosine(2 pi x)
	SinWave(NV_DATA_S(Ref), this->NumGridPts, this->NumGridPts, this->delx);
	SinWave(NV_DATA_S(this->Scrap), this->NumGridPts, this->NumGridPts*this->num_equations, this->delx);
	N_VAddConst(Ref, -1.0, Ref);
	N_VAddConst(this->Scrap, -1.0, this->Scrap);
	//scale								vthis is our reference
	N_VScale(2.0 * M_PI, Ref,Ref);
	N_VScale(2.0 * M_PI, this->Scrap,this->Scrap); 

	N_VLinearSum(1.0, this->CpGrad, -1.0, Ref, this->SmallScrap );
	if(N_VDotProd(this->SmallScrap, this->SmallScrap)< 10*this->delx )
	{
		std ::cout << "Congrats: CpGrad passed the sine test!\n";
	}
	else
	{
		std :: cout <<"Error not within tolerance!\n";
		failFlag =	1;
		failState ++;
	}
	std :: cout << "2-norm error from solution: " << N_VDotProd(this->SmallScrap, this->SmallScrap) << std::endl;
	//cosine rho
	N_VLinearSum(1.0, this->RhoGrad, -1.0, Ref, this->SmallScrap );
	if(N_VDotProd(this->SmallScrap, this->SmallScrap)< 10*this->delx )
	{
		std ::cout << "Congrats: RhoGrad passed the sine test!\n";
	}
	else
	{
		std :: cout <<"Error not within tolerance!\n";
		failFlag =	1;
		failState ++;
	}
	std :: cout << "2-norm error from solution: " << N_VDotProd(this->SmallScrap, this->SmallScrap) << std::endl;
	//Diff
	N_VLinearSum(1.0, this->DiffGrad, -1.0, this->Scrap, this->Scrap );

	if(N_VDotProd(this->Scrap, this->Scrap)< 10*this->delx )
	{
		std ::cout << "Congrats: Diffusion coefficient gradients passed the Sine test!\n";
	}
	else
	{
		std :: cout <<"Error not within tolerance!\n" << std :: endl;
		failFlag =	1;
		failState ++;
	}
	std :: cout << "2-norm error from solution: " << N_VDotProd(this->Scrap, this->Scrap) << std::endl;
	std :: cout << "completed sine tests!\n";
	std :: cout << "completed all tests!\n";
	std :: cout << "Report:  Failed unit tests: " << failState << std :: endl;

	return 0;
}

//---------------------------------------------------------------------------------------------------------                                                                                          
//---------------------------------------------------------------------------------------------------------                                                                                          
//---------------------------------------------------------------------------------------------------------                                                                                          
//---------------------------------------------------------------------------------------------------------                                                                                          
//---------------------------------------------------------------------------------------------------------                                                                                          
//---------------------------------------------------------------------------------------------------------                                                                                          
//---------------------------------------------------------------------------------------------------------                                                                                          
//---------------------------------------------------------------------------------------------------------                                                                                          
//    ___                                            
//   / __\ _ __   ___   ___  ___                     
//  / /   | '__| / _ \ / __|/ __|                    
// / /___ | |   | (_) |\__ \\__ \                    
// \____/ |_|    \___/ |___/|___/                    
                                                  
//     ___  _   __   __              _               
//    /   \(_) / _| / _| _   _  ___ (_)  ___   _ __  
//   / /\ /| || |_ | |_ | | | |/ __|| | / _ \ | '_ \ 
//  / /_// | ||  _||  _|| |_| |\__ \| || (_) || | | |
// /___,'  |_||_|  |_|   \__,_||___/|_| \___/ |_| |_|
                                                  
//  _____             _                              
// /__   \  ___  ___ | |_                            
//   / /\/ / _ \/ __|| __|                           
//  / /   |  __/\__ \| |_                            
//  \/     \___||___/ \__|          
// ooooo     ooo              o8o      .      ooooooooooooo                        .            
// `888'     `8'              `"'    .o8      8'   888   `8                      .o8            
//  888       8  ooo. .oo.   oooo  .o888oo         888       .ooooo.   .oooo.o .o888oo  .oooo.o 
//  888       8  `888P"Y88b  `888    888           888      d88' `88b d88(  "8   888   d88(  "8 
//  888       8   888   888   888    888           888      888ooo888 `"Y88b.    888   `"Y88b.  
//  `88.    .8'   888   888   888    888 .         888      888    .o o.  )88b   888 . o.  )88b 
//    `YbodP'    o888o o888o o888o   "888"        o888o     `Y8bod8P' 8""888P'   "888" 8""888P' 
//   .oooooo.                                          oooo                .                    
//  d8P'  `Y8b                                         `888              .o8                    
// 888           .ooooo.  ooo. .oo.  .oo.   oo.ooooo.   888   .ooooo.  .o888oo  .ooooo.         
// 888          d88' `88b `888P"Y88bP"Y88b   888' `88b  888  d88' `88b   888   d88' `88b        
// 888          888   888  888   888   888   888   888  888  888ooo888   888   888ooo888        
// `88b    ooo  888   888  888   888   888   888   888  888  888    .o   888 . 888    .o        
//  `Y8bood8P'  `Y8bod8P' o888o o888o o888o  888bod8P' o888o `Y8bod8P'   "888" `Y8bod8P'        
//                                           888                                                
//                                          o888o                      
//=============================//
//Sets transport grids and grad//
//=============================//
void myPb2::SetTransportGrid(N_Vector State)//One velocity grid for derivatives
{
	int index				= 0;
	realtype * LookupTemp	= NV_DATA_S(this->TempTable);
	realtype * LookupCp		= N_VGetArrayPointer(this->CPPoly);
	realtype * LookupDiff	= N_VGetArrayPointer(this->DiffTable);
	realtype * LookupRho	= N_VGetArrayPointer(this->RhoTable);
	realtype * CpData		= NV_DATA_S(this->CpGrid);
	realtype * RhoData		= NV_DATA_S(this->RhoGrid);
	realtype * DiffData 	= NV_DATA_S(this->DiffGrid);
	realtype * data 		= NV_DATA_S(State);
	realtype * GhostP		= NV_DATA_S(this->Ghost);
	//Set the temperature based grids
	for(int i = 0 ; i < this->NumGridPts ; i ++ )
	{
		index 				= this->TempTableLookUp(data[i], this->TempTable);
		//std :: cout << "o";
		CpData[i]			= LookupCp[index];//Check for a bug on this line with pointer dissappearing.
		RhoData[i]			= LookupRho[index];
		//std :: cout << "x";
		for( int j = 0 ; j < this->num_equations; j ++)
		{
			DiffData[j*this->NumGridPts + i] = LookupDiff[j*500 + index];

			//DiffData[j*this->num_equations + i] = LookupDiff[j*500 + index];
		}
	}
}

void myPb2::Test_RHS_CrossDiff(N_Vector TestState)
{
	std :: cout << " Beginning Cross-Diffusion test\n";

	this->SetTransportGrid(TestState);						//Sets the transport with respect to the current point's temperature
	// this->Set_ScalarGradient(TestState);
	// this->Set_TransportGradient(TestState);
	std :: cout << "Checking RHS_CrossDiff\n";
	this->RHS_CrossDiff(0, TestState, this->Tmp, this);
	if( N_VDotProd(this->Tmp, this->Tmp)<1e-5)
		std :: cout << "Passed initial state test\n";
	else
		std :: cout << "Failed initial state test\n";
	std :: cout << N_VDotProd(this->Tmp, this->Tmp) << std :: endl;




	std :: cout << "running a Cosine test\n";
	//rho, cp, state held constant, but Diffusion changed to cosine waves.
	//this->Set_TransportGradient(TestState);
	SinWave(NV_DATA_S(DiffGrid), this->NumGridPts, this->NumGridPts*this->num_equations, this->delx);

	//Manually fix the diffGrad boundary via the lookup table
	realtype * GP 			= NV_DATA_S(this->Ghost);
	realtype * GhostSetter 	= NV_DATA_S(DiffGrid);
	N_VConst(GhostSetter[0], DiffTable);
	this->Set_TransportGradient(TestState);
	this->Set_ScalarGradient(TestState);
	N_VConst(1.0, this->RhoTable);
	N_VConst(1.0, this->CPPoly);

	N_VConst(1.0, this->RhoGrid);
	N_VConst(1.0, this->CpGrid);

	N_VConst(0.0, this->CpGrad);
	N_VConst(0.0, this->RhoGrad);
	N_VConst(1.0, this->ScalarGradient);
	

	//Prerun print
	//N_VPrint_Serial(this->RhoGrid);
	this->RHS_CrossDiff(0, TestState, this->Tmp, this);
	//N_VPrint_Serial(this->Tmp);
	N_Vector Ref = N_VClone(this->Tmp);
	SinWave2(NV_DATA_S(Ref), this->NumGridPts, this->NumGridPts*this->num_equations, this->delx);
	N_VAddConst(Ref, -1.0, Ref);
	N_VScale(-2.0 * M_PI, Ref,Ref);
	
	//N_VPrint_Serial(this->Tmp);
	
	N_VLinearSum(1.0, this->Tmp, -1.0, Ref, Ref);
	//Ref now has the error check the error
	realtype Err = N_VDotProd(Ref, Ref);
	if(Err>10*this->delx)
	{
		std:: cout << "You have failed the const rho, cp, state variable diffusion grid test\n";
	}
	else
		std:: cout << "Error: " << Err << std:: endl;

}

void myPb2::Test_JtV_CrossDiff(N_Vector TestState)
{
	std :: cout << " Beginning Cross-Diffusion jtv test\n";

	this->SetTransportGrid(TestState);						//Sets the transport with respect to the current point's temperature
	this->Set_ScalarGradient(TestState);
	this->Set_TransportGradient(TestState);
	std :: cout << "Checking JtV_CrossDiff\n";
	N_VConst(1.0, this->Tmp);
	N_VConst(0.0, this->DiffGrad);
	this->JtV_CrossDiff(TestState, this->Tmp, 0, TestState, TestState, this, this->Tmp);
	if( N_VDotProd(this->Tmp, this->Tmp)<1e-5)
		std :: cout << "Passed initial state test\n";
	else
		std :: cout << "Failed initial state test\n";
	std :: cout << N_VDotProd(this->Tmp, this->Tmp) << std :: endl;




	std :: cout << "running a Cosine test\n";
	//rho, cp, state held constant, but Diffusion changed to cosine waves.
	//this->Set_TransportGradient(TestState);
	//manually set Diffgrid
	SinWave(NV_DATA_S(TestState), this->NumGridPts, this->NumGridPts*this->num_equations, this->delx);
	N_VConst(1.0, this->DiffGrad);
	N_VConst(1.0, this->ScalarGradient);
	N_VConst(1.0, this->RhoGrad);
	N_VConst(1.0, this->CpGrad);
	N_VConst(1.0, this->RhoGrid);
	N_VConst(1.0, this->CpGrid);
	this->JtV_CrossDiff(TestState, this->Tmp, 0, TestState, TestState, this, this->Tmp);
	
	//N_VPrint_Serial(this->Tmp);
	//Appears to give the correct response

	//Set the reference
	N_Vector Ref = N_VClone(this->Tmp);
	SinWave2(NV_DATA_S(Ref), this->NumGridPts, this->NumGridPts*this->num_equations, this->delx);
	N_VAddConst(Ref, -1.0, Ref);
	N_VScale(-2.0 * M_PI, Ref,Ref);
	//Manually fix the left boudaries due to stenciling
	realtype 	denom		= 1.0/(2*this->delx);
	for (int i  = 0 ; i < this->NumGridPts; i++)
	NV_DATA_S(Ref)[(i)*(this->NumGridPts)]= NV_DATA_S(TestState)[(i)*(this->NumGridPts)+1]*denom;

	//N_VPrint_Serial(this->Tmp);

	//N_VPrint_Serial(Ref);
	
	N_VLinearSum(1.0, this->Tmp, -1.0, Ref, Ref);
	//Ref now has the error check the error
	realtype Err = N_VDotProd(Ref, Ref);
	if(Err>10*this->delx)
	{
		std:: cout << "You have failed the const rho, cp, state variable diffusion grid test\n";
	}
	else
	{
		std :: cout << "Congrats!  Passed sine test with variable state and fixed transport\n";
	}
		std:: cout << "Error: " << Err << std:: endl;
		//N_VPrint_Serial(Ref);
		std :: cout << "Max error: " << N_VMaxNorm(Ref)<< std :: endl;
}

//=====================================================================
//Passed diagnostic second order convergence via eyeball norm 08/09/22
//=====================================================================
//Diagnostic testing for the scalar gradient grid.
int myPb2::Test_ScalarGradient(N_Vector TestState)
{
	//======================
	//Const 1
	//======================
	std :: cout << "Running Const 1 test\n";
	realtype * TSD 	= NV_DATA_S(TestState);
	realtype * GP   = NV_DATA_S(this->Ghost);
	realtype * RSD  = NV_DATA_S(this->Scrap);
	//N_VConst(1.0, TestState);
	//Modify the ghost points
	for( int i = 0 ; i < this->num_equations; i++)
	{
		GP[i]= 1;
	}
	N_VConst(1.0, TestState);
	this->Set_ScalarGradient(TestState);
	//N_VPrint_Serial(this->ScalarGradient);
	std :: cout << "\n\n\n\n\n Const Error: " << N_VMaxNorm(this->ScalarGradient) << std :: endl;
	
	
	//=================
	//Cosine(2 pi x) + 1
	//=================
	std :: cout << "Running Cosine test\n";
	//============v omg this was annoying==========================================//08/09/22
	SinWave(TSD, this->NumGridPts, this->NumGridPts*this->num_equations, this->delx);
	for( int i = 0 ; i < this->num_equations; i++)
	{
		GP[i]= TSD[0];//symmetry of cos
	}
	//N_VPrint_Serial(TestState);
	std :: cout << "\n\n\n\n";
	this->Set_ScalarGradient(TestState);

	//=================
	//Cosine(2 pi x) + 1
	//=================
	std :: cout << "Running Sine test\n";
	//============v omg this was annoying==========================================//08/09/22
	SinWave2(TSD, this->NumGridPts, this->NumGridPts*this->num_equations, this->delx);
	for( int i = 0 ; i < this->num_equations; i++)
	{
		GP[i]= -1*TSD[0];//symmetry of sine
	}
	//N_VPrint_Serial(TestState);
	//Note the right most point will only be half its expected value due to 0 boundary conditions.
	std :: cout << "\n\n\n\n";
	this->Set_ScalarGradient(TestState);

	return 0;
}

//====================================================================
// 888     888          888                   d8b 888             
// 888     888          888                   Y8P 888             
// 888     888          888                       888             
// Y88b   d88P  .d88b.  888  .d88b.   .d8888b 888 888888 888  888 
//  Y88b d88P  d8P  Y8b 888 d88""88b d88P"    888 888    888  888 
//   Y88o88P   88888888 888 888  888 888      888 888    888  888 
//    Y888P    Y8b.     888 Y88..88P Y88b.    888 Y88b.  Y88b 888 
//     Y8P      "Y8888  888  "Y88P"   "Y8888P 888  "Y888  "Y88888 
// 888     888               888          888                 888 
// 888     888               888          888            Y8b d88P 
// 888     888               888          888             "Y88P"  
// 888     888 88888b.   .d88888  8888b.  888888 .d88b.           
// 888     888 888 "88b d88" 888     "88b 888   d8P  Y8b          
// 888     888 888  888 888  888 .d888888 888   88888888          
// Y88b. .d88P 888 d88P Y88b 888 888  888 Y88b. Y8b.              
//  "Y88888P"  88888P"   "Y88888 "Y888888  "Y888 "Y8888           
//             888                                                
//             888                                                
//             888                                                
//===================================================
//Called by the main function in the integration loop
//===================================================
//Awaiting final test clear
//DoTo:  Clean up extra memory allocations, use existing data
//===================================================
//Called by:  Main but commented out.
void myPb2::UpdateOneDVel(N_Vector State)
{
	//Declare and clean
	N_Vector VTemp 			= N_VClone(this->VelAve);			//Scalar size
	N_Vector TempGrad		= N_VClone(this->Vel);				//Velocity size.
	N_Vector Scratch		= N_VClone(this->SmallChem);		//Used for right ghost
	N_VScale(0.0, VTemp, VTemp);								//Clean
	N_VScale(0.0, Tmp, Tmp);									//Clean
	N_VScale(0.0, TempGrad, TempGrad);							//Clean
	//Set relevent data pointers.
	realtype * SD			= N_VGetArrayPointer(State);
	realtype * VTempData	= N_VGetArrayPointer(VTemp);
	realtype * TmpD			= N_VGetArrayPointer(Tmp);			//Small Tmp data
	realtype * SCP			= N_VGetArrayPointer(this->SmallChem);
	realtype * GD			= N_VGetArrayPointer(this->Ghost);
	realtype * TG			= N_VGetArrayPointer(TempGrad);
	realtype * ScratchData	= N_VGetArrayPointer(Scratch);
	//Lookup data tables pointers
	realtype * LookupTemp	= N_VGetArrayPointer(this->TempTable);
	realtype * LookupCp		= N_VGetArrayPointer(this->CPPoly);
	realtype * LookupDiff	= N_VGetArrayPointer(this->DiffTable);
	realtype * LookupRho	= N_VGetArrayPointer(this->RhoTable);

	realtype LeftTemp		= 0;
	realtype RightTemp		= 0;
	int End 				= this->NumGridPts-1;			//The final vector entry
	int TInd				= 0;
	//See Dr. Bisetti's book, chapter 17 for full details.
	
	//Set the chemistry information:  omegaT / (rho * cp * Temp)
	SUPER_2_VEC(End, ScratchData, SD, this->num_equations, this->NumGridPts);//copy right ghost state from far right
	this->GhostChem(0, Scratch, this->SmallChem, this);			//Right Floating Chem using Scratch
	this->RHS_Chem(0, State, this->Tmp, this);
	//SUPER_CHEM_RHS_TCHEM(0, State, this->Tmp, this);			//Set post integration Chem term

	//calculate OmegaT/(T rho c_p)= SUPER_CHEM_RHS_TCHEM / T
	for(int i = 0 ; i < this->NumGridPts; i ++)				//Put the temp int VTemp
		VTempData[i]	= TmpD[i]/(SD[i]); 					//Divide by temperature
	LeftTemp  = SCP[0]/(GD[0]);								//May need to be zero //Set left boundary data, fixed point
	RightTemp = VTempData[End];								//Set right boundary 0 N conditions
	
	N_VScale(0.0, this->Tmp, this->Tmp);					//Clean out Tmp Vec for later use


	//New Derivative and flame front position.
	//Set grad(T)(x)
	/*
	this->TempGradient(State, TempGrad);				//Call the method to set TempGrad
	N_VScale(this->Diff, TempGrad, TempGrad);			//Scale TempGrad by Diff setting, 0 or 1 
	for(int i = 0; i < this->NumGridPts; i++)			//Loop over Temp indices.
	{
		TInd 	= this->TempTableLookUp(SD[i], this->TempTable);	//Find Temp lookup value
		TG[i]	= LookupDiff[ TInd ] / (LookupRho[ TInd ] * LookupCp [ TInd]  ) * TG[i]/SD[i];//Divide by the temperature.
	}
	TG[End+1]	= TG[End];
	*/
	//Set the diffusion term.
	this->RHS_Diff(0, State, this->Tmp, this);
	//SUPER_RHS_DIFF_CP(0, State, this->Tmp, this);
	N_VScale(this->Diff, TempGrad, TempGrad);			//Scale TempGrad by Diff setting, 0 or 1 
	for(int i = 0; i < this->NumGridPts; i++)			//Loop over Temp indices.
		VTempData[i]	+=	TmpD[i]/SD[i];				//New, adds the Temp Integral
	RightTemp+=VTempData[End];							//New, possibly problematic

	//Heating component
	//This will turn on/off depending on Temperature max value or time.
	if(this->HeatingOn==1)
	{
		this->RHS_Heat(0, State, this->Tmp, this);
		//SUPER_RHS_HEATING(0, State, this->Tmp, this);			//Calculate heating
		N_VScale(this->Power, this->Tmp, this->Tmp);			//Scale on/off
		this->HeatingRightGhost	= this->HeatingOn *
			this->Power*this->HeatingRightGhost;	//Scale on/off the ghost point
		for( int i =0 ; i < this->NumGridPts; i ++)
			VTempData[i] 	+= TmpD[i]/SD[i];		
		//Set boundaries
		LeftTemp += 0;											//Should always be zero
		RightTemp+= this->HeatingRightGhost/SD[End];			//Call the end boundary.
	}
	/**/
	//======================
	//Finalize
	//======================
	this->VelIntegrate(VTempData, State, LeftTemp, RightTemp);
	//V(x) = Integral( omega/(rho T c_p) ,[0,x] ) + V(0);		See VelIntegrate function
	N_VLinearSum(1.0, this->Vel, 1.0, TempGrad, this->Vel);		//V+=Grad(T)(x)
	//N_VAddConst(this->Vel, -1.0*TG[0], this->Vel);				//V-=Grad(T)(0), when Grad(T) method is used.
	this->SetVelAve();											//Modify VelAve
	if(N_VMin(this->Vel)<0  &&  abs(N_VMin(this->Vel)>1e-1) )
		std :: cout << "Warning: negative vel @" << this->t <<"\n";
	//Destroy Temp Vectors
	N_VDestroy(VTemp);
	N_VDestroy(TempGrad);
	N_VDestroy(Scratch);
	//Exit
}

// ================================
// The dead
// ================================
// Officially depreciated 10/04/22                                                                                      
// void MatrixVectorProduct(int number_of_equations, realtype * JacD, N_Vector x, N_Vector tmp, realtype * JV)
// {
//     realtype * TMP  = NV_DATA_S(tmp);
// 	realtype * X 	= NV_DATA_S(x);
// 	for (int i=0; i< number_of_equations; i++)//for each state
// 	{
// 		for(int j=0; j<number_of_equations; j++)//marches across the column
// 		{
// 			TMP[j]=JacD[j+i*number_of_equations];//Stored row major
// 			//JV[i] += X[j] * JacD[i * number_of_equations + j];
// 		}
// 		JV[i]=N_VDotProd(tmp,x);
// 	}
// }																						  \
// Officially depreciated 10/04/22                                                                                      
//void myPb2::Set_GasWeight(N_Vector State)
// {
// 	realtype * GWPtr 	= NV_DATA_S(this->GasWeight);
// 	realtype * SDPtr 	= NV_DATA_S(State);
// 	realtype * SMPtr 	= NV_DATA_S(this->MolarWeights);
// 	N_VConst(0.0, this->GasWeight);
// 	//March through each grid point and grab each mass fraction and multiply by the weight of each molecule
// 	for(int i = 0; i < this->NumGridPts; i ++)
// 	{
// 		for(int j = 0; j < this->num_equations; j++)
// 		{	//Sigma W_i*Y_i
// 			//update    =  current + W_i    *  Y_i
// 			//GWPtr[i]	= GWPtr[i] + SMPtr[j] * SDPtr[i + j * this->num_equations];
// 			GWPtr[i]	+= SMPtr[j] * SDPtr[i + j * this->NumGridPts];

// 		}
// 		//std :: cout << GWPtr[i] << std :: endl;
// 	}
// 	// std :: cout << std :: endl;	
// 	// std :: cout << N_VMaxNorm(this->GasWeight) << std :: endl;

// 	//N_VPrint_Serial(this->GasWeight);  //I don't understand why this line prints 0's

// }
//void myPb2::Test_GasWeight(N_Vector TestState)
// {
// 	//grab the grid length scratch vector and fill it with ones.
// 	N_VConst(1.0, this->SmallScrap);

// 	std :: cout << "Running Gas weight tests\n";

// 	//Tests the constant default state
// 	this->Set_GasWeight(TestState);
// 	//This should return constant check this:
// 	if(abs(N_VDotProd(this->GasWeight, this->SmallScrap) - N_VGetLength(this->GasWeight)* N_VMaxNorm(this->GasWeight))<this->delx)
// 		std :: cout << "Congrats, you passed the constant default state test\n";

// 	std :: cout << abs(N_VDotProd(this->GasWeight, this->SmallScrap) - N_VGetLength(this->GasWeight)* N_VMaxNorm(this->GasWeight)) << "\n";
	
// 	//Set the states to a cosine wave.  the sum of the grid should be 0 in this case.
// 	SinWave(NV_DATA_S(TestState), this->NumGridPts, this->NumGridPts*this->num_equations, this->delx);
// 	N_VAddConst(TestState, -1.0, TestState);
// 	this->Set_GasWeight(TestState);
// 	//Mimic the integral with the non-delx integeral, just sum the indices.
// 	if( N_VDotProd(this->GasWeight, this->SmallScrap) < 10*this->delx)
// 	{
// 		std :: cout << "Congrats, you passed the cosine state test\n";
// 	}
// 	std :: cout << N_VDotProd(this->GasWeight, this->SmallScrap) << std :: endl;
	
// }