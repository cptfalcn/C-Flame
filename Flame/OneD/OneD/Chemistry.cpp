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
myPb2::myPb2(ordinal_type num_eqs, real_type_1d_view_type work, WORK kmcd, int GridPts,
			N_Vector y,realtype Delx)
{
	this->num_equations= num_eqs;
	this->pb;//create the problem
	this->pb._p 	= 101325;			//Standard (1 atm) pressure
	this->pb._work	= work;
	this->pb._kmcd 	= kmcd;
	this->NumGridPts= GridPts;
	this->vecLength = num_eqs * GridPts;
	sunindextype Block = this->vecLength * this->vecLength;
	real_type_1d_view_type fac("fac", 2*vecLength);
	this->pb._fac   = fac;
	this->Jac	= N_VNew_Serial(Block);
	this->Mat	= SUNDenseMatrix(num_eqs*GridPts,num_eqs*GridPts);
	this->Vel 	= N_VNew_Serial( (GridPts+1)); //This is staggered bigger
	this->Ghost 	= N_VNew_Serial(num_eqs);
	N_VScale(1.0 , y, this->Ghost);
	this->Tmp 	= N_VNew_Serial(vecLength);
	this->delx 	= Delx;
	this->VelAve	= N_VNew_Serial(GridPts);
	this->LeftDiff	= 0;
	this->SmallScrap= N_VNew_Serial(GridPts);
	this->VelScrap	= N_VNew_Serial(GridPts+1);
	this->SmallChem	= N_VNew_Serial(num_eqs);
	this->HeatingOn = 1;
}

void myPb2::SetAdvDiffReacPow(realtype adv, realtype diff, realtype reac, realtype pow, bool up)
{
	this->Adv=adv;
	this->Diff=diff;
	this->React=reac;
	this->Power=pow;
	this->VelUp=up;
}

void myPb2::SetGhost(N_Vector y)
{
	N_VScale(1.0, y, this->Ghost);
}

void myPb2::SetLeftDiff(N_Vector State)
{
	realtype * Data	= N_VGetArrayPointer(State);
	this->LeftDiff 	= 1.0/(this->delx * this->delx) * (Data[1]-Data[0]);	// i+1 - 2 * i + i-1 / del^2
}

//=================
//Destructor
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
}

//===============================
//Depreciated, marked for removal
//===============================
void myPb2::PrintGuts(void)
{
	std :: cout << BAR << "General Parameters" << BAR << std :: endl;
	std :: cout << "Number of grid points & ghosts:" << this->NumGridPts << "\t\t";
	std :: cout << "Jacobians stored:" << N_VGetLength(this->Jac)/(this->num_equations*
					this->num_equations) << std :: endl;
}

//==============================
//Scale the pressure
//==============================
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
void myPb2::SetVels(int velLength, int velOpt)
{
	realtype * VelData 	= NV_DATA_S(this->Vel);
	realtype * VelAveData	= NV_DATA_S(this->VelAve);
	N_VConst(velOpt, this->Vel);
	this->SetVelAve();
}

//===============================
//Set Vel Average
//===============================
void myPb2::SetVelAve()
{
	realtype * VelData 	= N_VGetArrayPointer(this->Vel);
	realtype * VelAveData	= N_VGetArrayPointer(this->VelAve);
	for(int i = 0 ; i < this->NumGridPts-1; i++)
		VelAveData[i]   = (VelData[i] + VelData[i+1])/2;
}
//===========================
//Get TempGradient
//===========================
//Note here we step up to the Velocity grid.
void myPb2::TempGradient(N_Vector State, N_Vector Gradient)
{
	//Declare
	realtype * SD		= N_VGetArrayPointer(State);
	realtype * GD		= N_VGetArrayPointer(this->Ghost);
	realtype * GradD	= N_VGetArrayPointer(Gradient);
	int	MaxPos		= 0;
	N_VScale(0.0, Gradient, Gradient);				//Zero out
	//Loop
	GradD[0]=(SD[0]-GD[0])/this->delx;				// Boundary
	for( int i =1 ; i < this->NumGridPts-1; i ++)			// Treat boundary separately
	{
		GradD[i]= (SD[i]-SD[i-1])/this->delx;			// (T[i]-T[i-1])/Delx
		if(GradD[i]>GradD[MaxPos])
			MaxPos=i;
	}
	this->FlameFrontLocation=MaxPos;
}

//===================================================
//Called by the main function in the integration loop
//===================================================
//Still debugging
//===================================================
void myPb2::UpdateOneDVel(N_Vector State)
{//This might be bugged.  Need further investigation.
	//Declare and clean
	N_Vector VTemp 		= N_VClone(this->VelAve);		//Scalar size
	N_Vector TempGrad	= N_VClone(this->Vel);			//Velocity size.
	N_VScale(0.0, VTemp, VTemp);					//Clean
	N_VScale(0.0, Tmp, Tmp);					//Clean
	realtype * SD		= N_VGetArrayPointer(State);
	realtype * VTempData	= N_VGetArrayPointer(VTemp);
	realtype * TmpD		= N_VGetArrayPointer(Tmp);
	realtype * SCP		= N_VGetArrayPointer(this->SmallChem);
	realtype * GD		= N_VGetArrayPointer(this->Ghost);
	realtype * TG		= N_VGetArrayPointer(TempGrad);
	realtype LeftTemp	= 0;
	realtype RightTemp	= 0;
	//realtype x		= 0;
	int End 		= this->NumGridPts-1;			//The final vector entry

	//See Dr. Bisetti's book, chapter 17 for full details.
	//Begin Chem term
	this->GhostChem(0, this->Ghost, this->SmallChem, this);		//Set left boundary Chem term
	SUPER_CHEM_RHS_TCHEM(0, State, this->Tmp, this);		//Set post integration Chem term
	//The grids are different sizes so we need to manually move the data.
	//calculate Omega/(T rho c_p)
	for(int i = 0 ; i < this->NumGridPts; i ++)			//Put the temp int VTemp
		VTempData[i]	= TmpD[i]/SD[i];			//Divide by temperature

	LeftTemp  = SCP[0]/GD[0];//May need to be zero			//Set left boundary data, fixed point
	RightTemp = VTempData[End];					//Set right boundary 0 N conditions
	N_VScale(0.0, this->Tmp, this->Tmp);				//Zero out the Tmp Vector

	//Start Diff Term
	/*
	SUPER_RHS_DIFF(0, State, this->Tmp, this);			//Set regular State Diffuision
	N_VScale(this->Diff, this->Tmp, this->Tmp);			//Scale by Diff term
									//Different sized grids, loop
	for(int i = 0 ; i < this->NumGridPts; i ++)			//Place Temp Diff into VTemp
		VTempData[i]	+= TmpData[i];				//Only use regular diffusion.
	//Empty Right boundary						//0 Diffusion boundary condition
	LeftTemp += this->Diff * (SD[0] - GD[0])/(pow(this->delx,2));	//Set Left Boundary condtion.
	*/
	//New Derivative and flame front position.
	//Set grad(T)(x)
	this->TempGradient(State, TempGrad);
	N_VScale(this->Diff, TempGrad, TempGrad);

	//Start Heat source term
	/**/
	//This will turn on/off depending on Temperature max value or time.
	SUPER_RHS_HEATING(0, State, this->Tmp, this);			//Calculate heating
	N_VScale(this->Power, this->Tmp, this->Tmp);			//Scale appropriately
	this->HeatingRightGhost= this->Power*this->HeatingRightGhost;	//Scale the ghost point
	for( int i =0 ; i < this->NumGridPts; i ++)
		VTempData[i] 	+= TmpD[i]/SD[i];

	LeftTemp += 0;							//Should always be zero
	RightTemp+= this->HeatingRightGhost/SD[End];			//Call the end boundary
	/**/
	this->VelIntegrate(VTempData, State, LeftTemp, RightTemp);
	//V(x) = Integral( omega/(rho T c_p) ,[0,x] ) + V(0);		See VelIntegrate function
	N_VLinearSum(1.0, this->Vel, 1.0, TempGrad, this->Vel);		//V+=Grad(T)(x)
	N_VAddConst(this->Vel, -1.0*TG[0], this->Vel);			//V-=Grad(T)(0)
	this->SetVelAve();						//Modify VelAve
	//Destroy Temp Vectors
	N_VDestroy(VTemp);
	N_VDestroy(TempGrad);
	//Exit
}

//=====================================================================
//Verified on both Sine and Cosine waves to show second order in space.
//=====================================================================
//Velocity integrate function
//Solves:	V(x) = integral( f(x), [0,x] ) + V(0);
//Pass the off centered values as an N_Vector and it will calculate a second order integral.
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

//=============================================
//Verification tests for the Simpson Integrator
//=============================================
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
void myPb2::SetGhostPVel(N_Vector y, int Experiment, int SampleNum, realtype VelVal)
{
	this->SetGhost(y);
	this->ScaleP(Experiment, SampleNum);
	if(SampleNum ==10)//AdvDiffCoS problem
		this->SetVels(this->NumGridPts+1, 1);
	else
		this->SetVels(this->NumGridPts+1,VelVal);
}



//=====================================
//End Class functions
//===========================================================
//||  <:::::>  <::::>  <:::::>     <:::>         <::::>    ||
//||     .*/   |::     |::  ::>   <>   <>        |::  `>   ||
//||    .*/    |::::>  |:::::<   <>     <>  <::> |::   *>  ||
//||   .*/     |::     |:: ':>    <>   <>        |::  .>   ||
//||  <:::::>  <::::>  |::   :>    <:::>         <::::>    ||
//===========================================================
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
//=============================
//  _  _  ____  ___  ___  ____
// | |/ /| <> || _ \| _ \| <> |
// |   < |    ||  _/|  _/|    |
// |_|\_\|_||_||_|  |_|  |_||_|
//=============================
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



//=================================================
//||      <:::>   :>   ::  <::::>      <::::>    ||
//||     <>   <>  |:>  ::  |::         |::  `>   ||
//||    <>     <> |: > ::  |::::> <::> |::   *>  ||
//||     <>   <>  |:  >::  |::         |::  .>   ||
//||      <:::>   |:   >:  <::::>      <::::>    ||
//=================================================
//  ___  _  _  ___  ___  ___
// / __/| || || _ \|  _|| _ \
// \__ \| || ||  _/|  _||   <
// /___/ \__/ |_|  |___||_|_|
//
//========================================================
//These use the second version of the problem class, myPb2.
//========================================================
//Setting React=0 is not the same as removing it
int SUPER_RHS(realtype t, N_Vector State, N_Vector StateDot, void * UserData)
{
	myPb2 * pb{static_cast<myPb2 *> (UserData)};//Recast
	realtype * TmpData 	= NV_DATA_S(pb->Tmp);
	realtype * SDD 		= NV_DATA_S(StateDot);
	N_VScale(0.0, StateDot, StateDot);
	N_VScale(0.0, pb->Tmp, pb->Tmp);
	if(pb->React>0){//Chemistry RHS
		SUPER_CHEM_RHS_TCHEM(t, State, StateDot, UserData);
		N_VScale(pb->React, StateDot, StateDot);
	}

	SUPER_RHS_HEATING(t, State, pb->Tmp, UserData); 		//Heating, should be handled correctly
	N_VLinearSum(pb->Power, pb->Tmp, 1.0, StateDot, StateDot);

	if(pb->NumGridPts > 1)						//Skip if not enough points
	{
		SUPER_RHS_DIFF_NL(t, State, pb->Tmp, UserData);
		N_VLinearSum(pb->Diff, pb->Tmp, 1.0, StateDot, StateDot);
		SUPER_RHS_ADV_VEL(t,State,pb->Tmp,UserData);
		//N_VMaxNorm(pb->Tmp);
		N_VLinearSum(pb->Adv, pb->Tmp, 1.0, StateDot, StateDot);
		N_VScale(0.0, pb->Tmp, pb->Tmp);
		//SUPER_RHS_HEATING(t, State, pb->Tmp, UserData);
		//N_VMaxNorm(pb->Tmp);
		//N_VLinearSum(pb->Power, pb->Tmp, 1.0, StateDot, StateDot);
	}
	N_VDotProd(StateDot,StateDot);
	return 0;
}

//=============================
//Same issue as above
//=============================
int SUPER_JTV(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void * pb, N_Vector tmp)
{
	myPb2 * pbPtr{static_cast<myPb2 *> (pb)};//Recast

        N_VScale(0.0, Jv, Jv);
        N_VScale(0.0, tmp, tmp);
	if(pbPtr->React>0){
      		SUPER_CHEM_JTV(v, Jv, t, u, fu, pb, tmp);
		N_VScale(pbPtr->React, Jv, Jv);
	}
	if(pbPtr->NumGridPts > 1 )//skip if not enough points
	{
		SUPER_ADV_VEL_JTV(v, tmp, t, v, fu, pb, tmp);
		N_VLinearSum(pbPtr->Adv, tmp , 1.0, Jv, Jv);
		SUPER_DIFF_NL_JTV(v,tmp,t,u, fu, pb, tmp);
        	N_VLinearSum(pbPtr->Diff, tmp, 1.0, Jv, Jv);
	}
        return 0;
}


//===========================
//RHS
//===========================
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
//CHEM Stuff
//===============================
//Super Jac Via TChem
//===============================
int SUPER_CHEM_JAC_TCHEM(realtype t, N_Vector State, N_Vector StateDot, SUNMatrix Jac, void * pb, N_Vector tmp1,
				N_Vector tmp2, N_Vector tmp3)
{
	myPb2 * pbPtr{static_cast<myPb2 *> (pb)};//Recast
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
	}
	N_VDestroy_Serial(SmallJac);
	N_VDestroy_Serial(SmallTmp);
	N_VDestroy_Serial(SmallJV);
	return 0;
}


/*
//============================
//Vector Shenanagins
//==================================================================================================
//Do the clever MatVec product//bugged
//Tested on two problems, Jac is block diagonal ones, and Jac is block diagonal 0,1,2,...,n in rows.
//==================================================================================================
int CleverMatVec(int i, int num_eqs, int Block, realtype * JV, realtype * V, realtype * Jac)
{
	for(int j = 0 ; j < num_eqs; j ++)//row
		for( int k = 0; k < num_eqs; k++)//column
			JV[i * num_eqs + j ] += V[i * num_eqs + k] * Jac[i*Block + j * num_eqs + k];
	return 0;
}
*/
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
		SUPERDATA[i + j * grid_sz] = VECDATA[j];//Working?
		//SUPERDATA[j + i * grid_sz] = VECDATA[j];
	return 0;
}

/*
//Jac into SuperJac||bugged
int Jac_2_SuperJac(int i, realtype * Jac, realtype * SuperJac, int num_eqs, int grid_sz)
{
	int Block = num_eqs * grid_sz;
	for( int j = 0 ; j < num_eqs; j ++)
		for(int k = 0; k < num_eqs; k++)
			SuperJac[i*Block + j * num_eqs + k] = Jac[j*num_eqs + k];
	return 0;
}
*/
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

int RescaleTemp(realtype scaling, realtype * StateData, int numPts)
{//Rescales the temperature points
	for( int i = 0 ; i < numPts; i++)
		StateData[i]= scaling*StateData[i];
	return 0;
}

//==================================
//ADV Stuff
//==================================
//Super Adv
//==================================
int SUPER_RHS_ADV(realtype t, N_Vector State, N_Vector StateDot, void * userData)
{
	myPb2 * problem{static_cast<myPb2 *> (userData)};//Recast
	realtype delx = 1.0 / (problem->NumGridPts);
	realtype divisor = -1.0 / (2 * delx);
	realtype * uData	= N_VGetArrayPointer(State);
	realtype * resultData 	= N_VGetArrayPointer(StateDot);
	realtype * Ghost = NV_DATA_S(problem->Ghost);
	int numPts = problem->NumGridPts;
	int vecLength = problem->num_equations * numPts;
	int Grid = 0;
	for (int i = 0; i < vecLength ; i++)
        {
		Grid = floor(i/problem->NumGridPts);
                if( i % numPts == 0) //left boundary
                        resultData[i] = divisor * (uData[i + 1]- Ghost[Grid]); //fix boundary condition here
                else if( i % numPts == (numPts-1) )//right boundary
                        resultData[i] = divisor * (uData[i] - uData[i - 1]);
                else
                        resultData[i] = divisor* (uData[i + 1] - uData[i - 1]);
        }
	return 0;
}

//==================================
//Super Adv Vel
//==================================
int SUPER_RHS_ADV_VEL(realtype t, N_Vector State, N_Vector StateDot, void * userData)
{
        myPb2 * problem{static_cast<myPb2 *> (userData)};//Recast
        realtype delx = 1.0 / (problem->NumGridPts);
        realtype divisor = -1.0 / (2 * delx);
        realtype * uData        = N_VGetArrayPointer(State);
        realtype * resultData   = N_VGetArrayPointer(StateDot);
        realtype * Ghost = NV_DATA_S(problem->Ghost);
	realtype * VAD		= N_VGetArrayPointer(problem->VelAve);
	problem->SetVelAve();
        int numPts = problem->NumGridPts;
        int vecLength = problem->num_equations * numPts;
        int Grid = 0;
	int FI 		= 0;
        for (int i = 0; i < vecLength ; i++)
        {
                Grid = floor(i/problem->NumGridPts);
		FI   =	i % numPts;
                if( i % numPts == 0) //left boundary
                        resultData[i] = divisor * VAD[FI]* (uData[i + 1]- Ghost[Grid]); //fix boundary condition here
                else if( i % numPts == (numPts-1) )//right boundary
                        resultData[i] = divisor * VAD[FI]* (uData[i] - uData[i - 1]);
                else
                        resultData[i] = divisor* VAD[FI]* (uData[i + 1] - uData[i - 1]);
        }
        return 0;
}

int SUPER_ADV_JTV(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void* userData, N_Vector tmp)
{
        myPb2 * problem{static_cast<myPb2 *> (userData)};//Recast
        realtype delx 		= problem->delx;
        realtype divisor 	= -1.0 / (2 * delx);
        realtype * uData        = NV_DATA_S(u);//Stuff comes from here.
        realtype * resultData   = NV_DATA_S(Jv);//stuff goes in here.
        realtype * Ghost 	= NV_DATA_S(problem->Ghost);
	realtype * VelAveD	= NV_DATA_S(problem->VelAve);
        int numPts 		= problem->NumGridPts;
        int vecLength = problem->num_equations * problem->NumGridPts;
        int Grid = 0;
        for (int i = 0; i < vecLength; i++)
        {
                Grid = floor(i/problem->NumGridPts);
                if( i % numPts == 0) //left boundary
                        resultData[i] = divisor * uData[i + 1]; //fix boundary condition here
                else if( i % numPts == (numPts-1) )//right boundary
                        resultData[i] = divisor * (uData[i] - uData[i - 1]);
                else
                        resultData[i] = divisor* (uData[i + 1] - uData[i - 1]);
        }
        return 0;
}

//===========================
//Super ADV Vel JTV
//===========================
int SUPER_ADV_VEL_JTV(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void* userData, N_Vector tmp)
{
        myPb2 * problem{static_cast<myPb2 *> (userData)};//Recast
        realtype delx           = problem->delx;
        realtype divisor        = -1.0 / (2 * delx);
        realtype * uData        = NV_DATA_S(u);//Stuff comes from here.
        realtype * resultData   = NV_DATA_S(Jv);//stuff goes in here.
        realtype * Ghost        = NV_DATA_S(problem->Ghost);
        realtype * VelAveD      = NV_DATA_S(problem->VelAve);
        realtype * VAD          = N_VGetArrayPointer(problem->VelAve);
        int numPts              = problem->NumGridPts;
        int vecLength = problem->num_equations * problem->NumGridPts;
        int Grid = 0;
	int FI 	= 0;
        for (int i = 0; i < vecLength; i++)
        {
                Grid = floor(i/problem->NumGridPts);
		FI = i % numPts;
                if( i % numPts == 0) //left boundary
                        resultData[i] = divisor * VAD[FI] * uData[i + 1]; //fix boundary condition here
                else if( i % numPts == (numPts-1) )//right boundary
                        resultData[i] = divisor * VAD[FI] * (uData[i] - uData[i - 1]);
                else
                        resultData[i] = divisor * VAD[FI] * (uData[i + 1] - uData[i - 1]);
        }
        return 0;
}

//=========================
//Super Diff
//=========================
int SUPER_RHS_DIFF(realtype t, N_Vector u, N_Vector uDot, void * userData)
{
	myPb2 * problem{static_cast<myPb2 *> (userData)};//Recast
        realtype * uData        = NV_DATA_S(u);//Stuff comes from here.
        realtype * resultData   = NV_DATA_S(uDot);//stuff goes in here.
        int numPts = problem->NumGridPts;
	realtype delx 		= problem->delx;
	realtype divisor 	= 1.0/(delx * delx);
	int grid = 0;
	int tempInd = 0;
	realtype T = 1.0;
	realtype * Ghost = NV_DATA_S(problem->Ghost);
        int vecLength = problem->num_equations * problem->NumGridPts;
	realtype c = 1.0 * divisor;
	for (int i = 0; i < vecLength; i++)
        {
		tempInd = i%numPts;//Change later
		grid = floor(i/problem->NumGridPts);
                if(i% numPts == 0)//left
                        resultData[i] = T * c * (Ghost[grid] - 2*uData[i] + uData[i+1]);
                else if (i % numPts == (numPts - 1) )//right 0 neumann
                        resultData[i] = T * c * ( uData[i-1] - uData[i] );
                else
                        resultData[i] = T * c * (uData[i-1] - 2*uData[i] + uData[i+1]);
        }
	return 0;
}

int SUPER_RHS_DIFF_NL(realtype t, N_Vector u, N_Vector uDot, void * userData)
{
        myPb2 * problem{static_cast<myPb2 *> (userData)};//Recast
        realtype * uData        = NV_DATA_S(u);//Stuff comes from here.
        realtype * resultData   = NV_DATA_S(uDot);//stuff goes in here.
        int numPts = problem->NumGridPts;
        realtype delx           = problem->delx;
        realtype divisor        = 1.0/(delx * delx);
        int grid = 0;
        int tempInd = 0;
        realtype T = 1;
        realtype * Ghost = NV_DATA_S(problem->Ghost);
        int vecLength = problem->num_equations * problem->NumGridPts;
        realtype c = 1.0 * divisor;
        for (int i = 0; i < vecLength; i++)
        {
                tempInd = i%numPts;
                T = uData[tempInd];//This adds the non-linearity we wanted
                grid = floor(i/problem->NumGridPts);
                if(i% numPts == 0)//left
                        resultData[i] = T * c * (Ghost[grid] - 2*uData[i] + uData[i+1]);
                else if (i % numPts == (numPts - 1) )//right 0 neumann
                        resultData[i] = T * c * ( uData[i-1] - uData[i] );
                else
                        resultData[i] = T * c * (uData[i-1] - 2*uData[i] + uData[i+1]);
        }
        return 0;
}

int SUPER_DIFF_JTV(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void* userData, N_Vector tmp)
{
        myPb2 * problem{static_cast<myPb2 *> (userData)};//Recast
        realtype * uData        = NV_DATA_S(u);//Stuff comes from here.
        realtype * resultData   = NV_DATA_S(Jv);//stuff goes in here.
        int numPts = problem->NumGridPts;
        realtype delx           = problem->delx;
        realtype divisor        = 1.0/(delx * delx);
        int tempInd = 0;
	int grid = 0;
        realtype * Ghost = NV_DATA_S(problem->Ghost);
        int vecLength = problem->vecLength;
        realtype c = 1.0 * divisor;
        for (int i = 0; i < vecLength; i++)
        {
		grid = floor(i/problem->NumGridPts);
                tempInd = i%numPts;//Change later
                if(i% numPts == 0)//left
                        resultData[i] = c * (-2*uData[i] + uData[i+1]);
                else if (i % numPts == (numPts - 1) )//right 0 neumann
                        resultData[i] = c * ( uData[i-1] - uData[i] );
                else
                        resultData[i] = c * (uData[i-1] - 2*uData[i] + uData[i+1]);
        }
        return 0;
}

int SUPER_DIFF_NL_JTV(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu,void* userData, N_Vector tmp)
{
        myPb2 * problem{static_cast<myPb2 *> (userData)};//Recast
        realtype * uData        = NV_DATA_S(u);//Stuff comes from here.
        realtype * resultData   = NV_DATA_S(Jv);//stuff goes in here.
	realtype * vData 	= NV_DATA_S(v);//Stuff also comes from here
        int numPts = problem->NumGridPts;
        realtype delx           = problem->delx;
        realtype divisor        = 1.0/(delx * delx);
        int tempInd 		= 0;
	int grid 		= 0;
	realtype * Ghost = NV_DATA_S(problem->Ghost);
	realtype T 		= 0;
        int vecLength 		= problem->vecLength;
        realtype c 		= 1.0 * divisor;
        for (int i = 0; i < vecLength; i++)
        {
		tempInd = i%numPts;
		grid = floor(i/problem->NumGridPts);
                T = uData[tempInd];//This adds the non-linearity we wanted
                if(i% numPts == 0)//left
		{
                        resultData[i] = c * T * (-2*vData[i] + vData[i+1]);
			resultData[i]+= c * vData[tempInd] * (Ghost[grid] - 2 *uData[0] + uData[1]);
		}
                else if (i % numPts == (numPts - 1) )//right 0 neumann
		{
                        resultData[i] = c * T * ( vData[i-1] - vData[i] );
			resultData[i]+= c * vData[tempInd] *  (uData[i-1] - uData[i]);
		}
                else
		{
                        resultData[i] = c * T * (vData[i-1] - 2*vData[i] + vData[i+1]);
			resultData[i]+= c * vData[ tempInd ] * (uData[ i-1 ]-2 * uData[ i ]+ uData[ i-1 ] );
		}
        }
        return 0;
}

//
int SUPER_RHS_HEATING(realtype t, N_Vector u, N_Vector uDot, void * userData)
{
        myPb2 * problem{static_cast<myPb2 *> (userData)};//Recast
        realtype * uData        = NV_DATA_S(u);//Stuff comes from here.
        realtype * resultData   = NV_DATA_S(uDot);//stuff goes in here.
        int numPts = problem->NumGridPts;
        realtype delx           = problem->delx;
        realtype * Ghost = NV_DATA_S(problem->Ghost);
        int vecLength = problem->NumGridPts;//Only march through temp
	realtype x = 0;
	realtype OffSet = 0;
	realtype xEnd= vecLength*delx+delx/2;
	if(t>1e-6 || N_VMaxNorm(u)>3500)//N_VMaxNorm(u)>2200 & problem->HeatingOn == 1)
		problem->HeatingOn=0;
	if(problem->HeatingOn==1 )//t<1e-5 && N_VMaxNorm(u) <  2200)//2000 originally
	{
		//OffSet = xEnd/2;
		OffSet = delx*(round(0.9*problem->NumGridPts) - 0.5);	//Set 1/10 to left
		//OffSet = delx*(problem->NumGridPts - 0.5);
        	for (int i = 0; i < numPts; i++) //Only Heat the Temp
		{
			x = i*delx + delx/2;
			resultData[i] = 5e10 * exp( -1e8 * pow(x - OffSet, 2) );
			//resultData[i] = 5e10 * exp( -1e8 * pow( delx * (i+1-problem->NumGridPts) , 2 ) );
		}
		problem->HeatingRightGhost=5e10 * exp( -1e8 * pow(x+delx/2 - OffSet, 2) );
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
	CHEM_RHS_TCHEM(t, y, ydot, pb);
}
