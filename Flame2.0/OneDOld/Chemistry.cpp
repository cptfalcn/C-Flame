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
//==================================
myPb2::myPb2(ordinal_type num_eqs, real_type_1d_view_type work, WORK kmcd, int GridPts,
			N_Vector y,realtype Delx)
{
	this->num_equations= num_eqs;
	//Create TChem Pb
	this->pb;								//create the problem
	this->pb._p 	= 101325;				//Standard (1 atm) pressure
	this->pb._work	= work;
	this->pb._kmcd 	= kmcd;
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
	this->ScalarGradient= N_VClone(this->Scrap);
	//Options
	this->HeatingOn = 1;
	this->dumpJac	= 0;								//Do we dump the Jac for visualization
	//Transport grids
	this->CpGrid	= N_VNew_Serial(GridPts);		//The value of Cp based on grid position.
	this->RhoGrid	= N_VNew_Serial(GridPts);		//The value of Rho based on grid position.
	this->DiffGrid	= N_VNew_Serial(vecLength);			//Diffusion coefficient based on grid position.
	this->MolarWeights = N_VNew_Serial(GridPts-1);		//Only grab the species molar wieghts;
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

	// !!!This kills the 0-D problem!!!
	// realtype * 			MolarWeightsPtr =	NV_DATA_S(MolarWeights);
	// real_type_1d_view_type SpeciesMolecularWeights(MolarWeightsPtr   ,   GridPts-1);
	// SpeciesMolecularWeights = Kokkos::create_mirror_view(this->pb._kmcd.sMass);
	// Kokkos::deep_copy(SpeciesMolecularWeights, this->pb._kmcd.sMass);
	
	// //std :: cout << kmcd.sMass << std :: endl;
	// //const auto SpeciesMolecularWeights = Kokkos::create_mirror_view(kmcd.sMass);
	// std::setprecision(17);
	// std :: cout << MolarWeightsPtr[0] << " zero-th mass entry?" << std :: endl;
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
	realtype * VelData 		= NV_DATA_S(this->Vel);
	realtype * VelAveData	= NV_DATA_S(this->VelAve);
	N_VConst(velOpt, this->Vel);
	this->SetVelAve();
}

//===============================
//Set Vel Average
//===============================
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

//===================================================
//Called by the main function in the integration loop
//===================================================
//Awaiting final test clear
//DoTo:  Clean up extra memory allocations, use existing data
//===================================================
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
	if(SampleNum ==10|| SampleNum == 12 )//AdvDiffCoS problem
		this->SetVels(this->NumGridPts+1, 1);
	else
		this->SetVels(this->NumGridPts+1,VelVal);
}

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



void myPb2::Set_RHS(CVRhsFn OneD_RHS_Chem, CVRhsFn OneD_RHS_Adv, CVRhsFn OneD_RHS_Diff, CVRhsFn OneD_RHS_Heat)
{
	this->RHS_Adv		=OneD_RHS_Adv;
	this->RHS_Diff		=OneD_RHS_Diff;
	this->RHS_Chem		=OneD_RHS_Chem;
	this->RHS_Heat 		=OneD_RHS_Heat;
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
//  __ ==============_===================___==============_================_====_===================                   
// / _\  ___   __ _ | |  __ _  _ __     /   \  ___  _ __ (_)__   __  __ _ | |_ (_)__   __  ___  ___ 
// \ \  / __| / _` || | / _` || '__|   / /\ / / _ \| '__|| |\ \ / / / _` || __|| |\ \ / / / _ \/ __|
// _\ \| (__ | (_| || || (_| || |     / /_// |  __/| |   | | \ V / | (_| || |_ | | \ V / |  __/\__ \
// \__/ \___| \__,_||_| \__,_||_|    /___,'   \___||_|   |_|  \_/   \__,_| \__||_|  \_/   \___||___/
//==================================================================================================
//==================================================================================================      
//===============================================
//Sets scalar Temp and mass fraction gradients.
//===============================================
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


//This derivative is on the scalar grid
//Should be called only once per step.
//Set Gradients of the transport properties.  Called by cross diffusion.
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
	realtype * GradT	= NV_DATA_S(this->TempGrad);
	realtype * GradRho	= NV_DATA_S(this->RhoGrad);
	realtype * GradCp	= NV_DATA_S(this->CpGrad);
	realtype * GradDiff	= NV_DATA_S(this->DiffGrad);


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
		CpData[i]			= LookupCp[index];
		RhoData[i]			= LookupRho[index];
		for( int j = 0 ; j < this->num_equations; j ++)
		{
			DiffData[j*this->num_equations + i] = LookupDiff[j*500 + index];
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

