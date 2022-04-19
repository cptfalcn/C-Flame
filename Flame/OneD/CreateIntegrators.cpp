#include "CreateIntegrators.h"
#include "Chemistry.h"
#include <cmath>

//  __    __   __    ______    ____
// |  >> |  >>|  >> /      >> |   _>>
// |  >> |   >|  >> \_   ..>> |  |_
// |  >> |       >>   |  >>   |   _>>
// |  >> |  |\   >>   |  >>   |  |_
// |__>> |__| \__>>   |__>>   |____>>
//

//=====================
//Create a EPI2 Pointer
//=====================
Epi2_KIOPS* CreateEPI2Integrator(CVRhsFn RHS, CVSpilsJacTimesVecFn JtV, void * pbptr, int MaxKrylovIters,
                                N_Vector y, const int number_of_equations, int UseJac)
{
                switch(UseJac)
                {
                case 0:{
			/*
                        Epi2_KIOPS *integrator = new Epi2_KIOPS(RHS, pbptr,
                                MaxKrylovIters, y, number_of_equations);
			*/
			/**/
			unique_ptr<Epi2_KIOPS> EPI2_PTR (	new Epi2_KIOPS(RHS, pbptr,
					MaxKrylovIters, y, number_of_equations));//Smart Pointer
			/**/
                        cout<< BAR <<"\tEPI2 w/o JtV\t\t"<< BAR << endl;
                        //return integrator;
			return EPI2_PTR.release();
                        break;}
                case 1:{
                        Epi2_KIOPS *integrator = new Epi2_KIOPS(RHS, JtV, pbptr,
                                MaxKrylovIters, y, number_of_equations);

                        cout<< BAR << "\tEPI2 JtV\t\t"<< BAR << endl;
                        return integrator;
                        break;}
                default:
                        cout<< BAR <<"\t Invalid EPI2 Jacobian Option\t"<< BAR << endl;;
                        exit(EXIT_FAILURE);
                }
        return 0;
}

//======================
//Create an EPI3 Pointer
//======================
Epi3_KIOPS* CreateEPI3Integrator(CVRhsFn RHS, CVSpilsJacTimesVecFn JtV, void * pbptr,
		int MaxKrylovIters, N_Vector y, const int number_of_equations, int UseJac)
{
                switch(UseJac)
                {
                case 0:{
                        Epi3_KIOPS *integrator = new Epi3_KIOPS(RHS, pbptr,
                                MaxKrylovIters, y, number_of_equations);

                        cout<< BAR <<"\tEPI3 w/o JtV\t\t"<< BAR << endl;
                        return integrator;
                        break;}
                case 1:{
                        Epi3_KIOPS *integrator = new Epi3_KIOPS(RHS, JtV, pbptr,
				MaxKrylovIters, y, number_of_equations);

                        cout<< BAR <<"\tEPI3 JtV\t\t"<< BAR << endl;
                        return integrator;
                        break;}
                default:
                        cout<<"===============EPI3 Jacobian  Option Invalid============\n";
                        exit(EXIT_FAILURE);
                }
        return 0;
}


//==========================
//Create an EPIRK4SC Pointer
//==========================
EpiRK4SC_KIOPS * CreateEPIRK4SCIntegrator(CVRhsFn RHS, CVSpilsJacTimesVecFn JtV, void * pbptr,
		 int MaxKrylovIters, N_Vector y, const int number_of_equations, int UseJac)
{
                switch(UseJac)
                {
                case 0:{
                        EpiRK4SC_KIOPS *integrator = new EpiRK4SC_KIOPS(RHS,pbptr,
                                MaxKrylovIters, y, number_of_equations);

                        cout<< BAR <<"\tEPIRK4 w/o JtV\t\t"<< BAR << endl;
                        return integrator;
                        break;}
                case 1:{
                        EpiRK4SC_KIOPS *integrator = new EpiRK4SC_KIOPS(RHS, JtV, pbptr,
                                MaxKrylovIters, y, number_of_equations);

                        cout<< BAR <<"\tEPIRK4 JtV\t\t"<< BAR << endl;
                        return integrator;
                        break;}
                default:
                        cout<< BAR << "EPIRK4 Jacobian  Option Invalid" << BAR << endl;
                        exit(EXIT_FAILURE);
                }
        return 0;
}

//======================================================
//General Cvode with general RHS and Jacobian functions.
//======================================================


void * CreateCVODE(CVRhsFn RHS, CVSpilsJacTimesVecFn JtV, CVLsJacFn JAC,  void * pbptr, SUNMatrix A,
                        SUNLinearSolver LS, SUNNonlinearSolver NLS , int num_eqs, N_Vector State,
                        realtype relTol, realtype absTol, realtype StepSize, int UseJac)
{
        void * cvode_mem;
        int retVal = 0;
        //realtype tret=0;
        N_Vector AbsTol= N_VNew_Serial(num_eqs);
        for ( int i = 0 ; i < num_eqs ; i++)
                NV_Ith_S(AbsTol,i)=absTol;

        cvode_mem = CVodeCreate (CV_BDF);
        //Give CVODE the user problem
        retVal = CVodeSetUserData(cvode_mem, pbptr);
        //Set intial state
        retVal = CVodeInit(cvode_mem, RHS, 0, State);
        //Set intial stepsize guess
        retVal = CVodeSetInitStep(cvode_mem, StepSize);
        //Set max stepsize
        retVal = CVodeSetMaxStep(cvode_mem, StepSize);
	//Set min stepsize
  	//retVal = CVodeSetMinStep(cvode_mem, StepSize/1e-4);
        //Set tolerances
        retVal = CVodeSVtolerances(cvode_mem, relTol, AbsTol);
        //Set linear solver
        retVal = CVodeSetLinearSolver(cvode_mem, LS, A);//A might be null, try that
        retVal = CVodeSetNonlinearSolver(cvode_mem, NLS);
        if (UseJac==1)
        {
		retVal = CVodeSetJacTimes(cvode_mem, NULL, JtV);
		retVal = CVodeSetLSetupFrequency(cvode_mem, 1); //Remove because also stupid.
		retVal = CVodeSetJacEvalFrequency(cvode_mem, 1);//Remove becuase this is stupid.
                retVal = CVodeSetJacFn(cvode_mem, JAC);		//Error if removed .
                retVal = CVodeSetMaxNumSteps(cvode_mem, 500);
        }
        N_VDestroy_Serial(AbsTol);
        return cvode_mem;
}


//===============================
//Set-up an Adaptive CVODE solver
//High max steps
//===============================
void * CreateCVODEAdp(CVRhsFn RHS, CVSpilsJacTimesVecFn JtV, CVLsJacFn JAC,  void * pbptr, SUNMatrix A,
                        SUNLinearSolver LS, SUNNonlinearSolver NLS , int num_eqs, N_Vector State,
                        realtype relTol, realtype absTol, realtype StepSize)
{
        void * cvode_mem;
        int retVal = 0;
        N_Vector AbsTol= N_VNew_Serial(num_eqs);
        for ( int i = 0 ; i < num_eqs ; i++)
                NV_Ith_S(AbsTol,i)=absTol;

        cvode_mem = CVodeCreate (CV_BDF);
        //Give CVODE the user problem
        retVal = CVodeSetUserData(cvode_mem, pbptr);
        //Set intial state, 	(mem, RHS, t_0, y_0)
        retVal = CVodeInit(cvode_mem, RHS, 0, State);
        //Set intial stepsize guess
        retVal = CVodeSetInitStep(cvode_mem, StepSize);
        //Set max stepsize
        retVal = CVodeSetMaxStep(cvode_mem, 100*StepSize);
        //Set min stepsize
        //retVal = CVodeSetMinStep(cvode_mem, StepSize/1e-4);
        //Set tolerances
        retVal = CVodeSVtolerances(cvode_mem, relTol, AbsTol);

	//Set linear solver
        retVal = CVodeSetLinearSolver(cvode_mem, LS, A);//A might be null, try that
        retVal = CVodeSetNonlinearSolver(cvode_mem, NLS);
        retVal = CVodeSetJacTimes(cvode_mem, NULL, JtV);
        retVal = CVodeSetLSetupFrequency(cvode_mem, 1); //Remove because also stupid.
        retVal = CVodeSetJacEvalFrequency(cvode_mem, 1);//Remove becuase this is stupid.
        retVal = CVodeSetJacFn(cvode_mem, JAC);         //Error if removed .
        retVal = CVodeSetMaxNumSteps(cvode_mem, 1e6);	//Max steps allowed
        N_VDestroy_Serial(AbsTol);
        return cvode_mem;
}




//=======================
//Create CVODE BDF Krylov
//=======================

void * CreateCVODEKrylov(CVRhsFn RHS, CVSpilsJacTimesVecFn JtV, void * pbptr, SUNMatrix A,
			SUNLinearSolver LS, SUNNonlinearSolver NLS , int num_eqs, N_Vector State,
			realtype relTol, realtype absTol, realtype StepSize, int UseJac)
{
	void * cvode_mem;
	int retVal = 0;
	realtype tret=0;
	N_Vector AbsTol= N_VNew_Serial(num_eqs);
	for ( int i = 0 ; i < num_eqs ; i++)
			NV_Ith_S(AbsTol,i)=absTol;

	cvode_mem = CVodeCreate (CV_BDF);
	//Give CVODE the user problem
	retVal = CVodeSetUserData(cvode_mem, pbptr);
	//Set intial state
	retVal = CVodeInit(cvode_mem, RHS, 0, State);
	//Set intial stepsize guess
	retVal = CVodeSetInitStep(cvode_mem, StepSize);
	//Set max stepsize
	retVal = CVodeSetMaxStep(cvode_mem, StepSize);
	//retVal = CVodeSetMinStep(cvode_mem, StepSize);
	//Set tolerances
	retVal = CVodeSVtolerances(cvode_mem, relTol, AbsTol);
	//Set linear solver
	retVal = CVodeSetLinearSolver(cvode_mem, LS, A);
	retVal = CVodeSetNonlinearSolver(cvode_mem, NLS);
	retVal = CVodeSetJacTimes(cvode_mem, NULL, JtV);
	if (UseJac==1)
	{
	retVal = CVodeSetJacFn(cvode_mem, SUPER_CHEM_JAC_TCHEM);//Updated Version
	//retVal = CVodeSetJacTimes(cvode_mem, NULL, JtV);
			retVal = CVodeSetMaxNumSteps(cvode_mem, max(500,num_eqs*num_eqs*num_eqs));
	}
	N_VDestroy_Serial(AbsTol);
	return cvode_mem;
}

//========================
//Integrate: Marked for removal
//========================
IntegratorStats* Integrate(int UseJac, realtype StepSize, realtype TNow, realtype TNext, int NumBands,
				N_Vector y, realtype KrylovTol, int startingBasisSizes[], void * pbptr,
				void * cvode_mem, string Method, Epi2_KIOPS * integrator,
				Epi3_KIOPS * integrator2, EpiRK4SC_KIOPS * integrator3,
				IntegratorStats *  IntStats)
{
	if(UseJac==1 && Method != "CVODEKry")
		CHEM_COMP_JAC(y,pbptr);
	if(Method=="EPI2")
	{
		IntStats = integrator->Integrate(StepSize, TNow, TNext, NumBands, y,
				KrylovTol, startingBasisSizes);
	}else if(Method == "EPI3")
        {
		IntStats = integrator2->Integrate(StepSize, TNow, TNext, NumBands, y, KrylovTol,
				startingBasisSizes);
	}else if(Method =="EPIRK4")
        {
		IntStats = integrator3->Integrate(StepSize, TNow, TNext, NumBands, y,
				KrylovTol, startingBasisSizes);
     	}else if(Method=="CVODEKry")
	{
		CVode(cvode_mem, TNow, y, &TNext, CV_NORMAL);
		//CVode(cvode_mem, TNow, y, &TNext, CV_ONE_STEP);
	}
	return IntStats;
}

//=======================================================
//  ___  _  _  ___  ___  ___
// / __/| || || _ \|  _|| _ \
// \__ \| || ||  _/|  _||   <
// /___/ \__/ |_|  |___||_|_|
//
//========================================================
//Create the void integrator: Marked for removal
//=====================================
void *		SelectIntegrator(CVRhsFn RHS, CVSpilsJacTimesVecFn JtV, CVLsJacFn Jac, void * pbptr,
				void * integrator, int MaxKrylovIters, N_Vector y, const int num_eqs,
				int UseJac, string Method, SUNMatrix A, SUNLinearSolver LS,
				SUNNonlinearSolver NLS, realtype absTol, realtype relTol, realtype StepSize)
{
		realtype TNext= 1e-6;
		int startingBasisSizes[] ={3,3};
		IntegratorStats * IntStats = NULL;
		if(Method=="EPI2")
		{
		integrator = reinterpret_cast<void *>
				(CreateEPI2Integrator(RHS, JtV, pbptr, MaxKrylovIters, y, num_eqs, UseJac));
		}
		else if(Method == "EPI3")
		{
		integrator = reinterpret_cast<void *>
				( CreateEPI3Integrator(RHS, JtV, pbptr, MaxKrylovIters, y, num_eqs, UseJac));
		}
		else if(Method == "EPIRK4")
		{
		integrator = reinterpret_cast<void *>
			( CreateEPIRK4SCIntegrator(RHS, JtV, pbptr, MaxKrylovIters, y, num_eqs, UseJac));
		}else//create CVODE
		{
			integrator = CreateCVODE(RHS, JtV, Jac, pbptr, A, LS, NLS , num_eqs, y,
						relTol, absTol, StepSize, UseJac);
			//integrator = CreateCVODEKrylov(RHS, JtV, pbptr, A, LS, NLS,
			//			num_eqs,y,relTol, absTol, StepSize, UseJac);
		}
		return integrator;
}
//Marked for removal
IntegratorStats * IntegrateWrapper(int UseJac, realtype StepSize, realtype TNow, realtype TNext,
				int NumBands, N_Vector y, realtype KrylovTol, int startingBasisSizes[],
				void * pbptr, void * ProtoInt, void * cvode_mem,string Method,
				IntegratorStats *  IntStats, CVLsJacFn JacFn)
{
	N_Vector StateDot = N_VClone(y);
	int Length = N_VGetLength(y);
	myPb2 * pbPtr{static_cast<myPb2 *> (pbptr)};
	realtype * STATEDATA = NV_DATA_S(y);
	if(UseJac==1 && Method != "CVODEKry")
		JacFn(TNow, y, StateDot, pbPtr->Mat, pbPtr, StateDot, StateDot, StateDot);
        if(Method=="EPI2")
        {
		Epi2_KIOPS * integrator{static_cast<Epi2_KIOPS *> (ProtoInt)};
               	IntStats = integrator->Integrate(StepSize, TNow, TNext, NumBands, y, KrylovTol,	startingBasisSizes);
        }else if(Method == "EPI3")
        {
		Epi3_KIOPS * integrator{static_cast<Epi3_KIOPS *> (ProtoInt)};
               	IntStats = integrator->Integrate(StepSize, TNow, TNext, NumBands, y, KrylovTol, startingBasisSizes);
        }else if(Method =="EPIRK4")
        {
		EpiRK4SC_KIOPS * integrator{static_cast<EpiRK4SC_KIOPS *> (ProtoInt)};
               	IntStats = integrator->Integrate(StepSize, TNow, TNext, NumBands, y, KrylovTol, startingBasisSizes);
        }else if(Method=="CVODEKry")
        {
               	CVode(cvode_mem, TNow, y, &TNext, CV_NORMAL);
               	//CVode(cvode_mem, TNow, y, &TNext, CV_ONE_STEP);
        }
	return  IntStats;
}


//======================================================
//  __    __   __    ______    ____     ____    ____
// |  >> |  >>|  >> /      >> |   _>>  /  _ >> |    >>
// |  >> |   >|  >> \_   ..>> |  |_   |  / \_>>|  x  >>
// |  >> |       >>   |  >>   |   _>> | (  ___ |    >>
// |  >> |  |\   >>   |  >>   |  |_   |  \_/ >>| |\  >>
// |__>> |__| \__>>   |__>>   |____>>  \____>> |_| \__>>
//
//======================================================
// ___  ___  ____  ____
//| __|| _ \|_  _||__  |
//| __||  _/ _||_ |__  |
//|___||_|  |____||____|
//=============================
//New TCHEM adaptive integrator
//=============================

Epi3VChem_KIOPS::Epi3VChem_KIOPS(CVRhsFn f, CVSpilsJacTimesVecFn jtv, void *userData, int maxKrylovIters,
					N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    	this->f = f;
    	this->jtv = jtv;
    	this->userData = userData;

    	krylov = new Kiops(MaxPhiOrder1 + 2, maxKrylovIters, tmpl, NEQ);//Originally MaxPhiOrder1+1
    	integratorStats = new IntegratorStats(NumProjectionsPerStep);


        realtype hMax=0;
    	fy = N_VClone(tmpl);
        Y1 = N_VClone(tmpl);
        fY1 = N_VClone(tmpl);
    	hfy = N_VClone(tmpl);
    	hb1 = N_VClone(tmpl);
    	hb2 = N_VClone(tmpl);
    	r1 = N_VClone(tmpl);
        r2=N_VClone(tmpl);
        r3=N_VClone(tmpl);
        Remainder = N_VClone(tmpl);
        Scratch1  = N_VClone(tmpl);
    	tmpVec = N_VClone(tmpl);
    	zeroVec = N_VClone(tmpl);
    	// Zero out zeroVec by using the trusted vector tmpl.
    	// If zeroVec has random data in it, it may have NAN values.
    	// zeroVec = 0.0 * zeroVec may accidentally be the calculation
    	// 0.0 * NAN = NAN, and zeroVec will not end up zero.
    	// Instead zeroVec = 0.0 * templ, since templ is trusted
    	// to have no NAN values in it.
    	N_VConst(0, fy);
    	N_VConst(0, Y1);
    	N_VConst(0, fY1);
    	N_VConst(0, hfy);
    	N_VConst(0, hb1);
    	N_VConst(0, hb2);
    	N_VConst(0, r1);
    	N_VConst(0, zeroVec);
    	N_VConst(0, r2);
    	N_VConst(0, r3);
        N_VConst(0, Scratch1);
}


Epi3VChem_KIOPS::Epi3VChem_KIOPS(CVRhsFn f, void *userData, int maxKrylovIters, N_Vector tmpl,
				const int vecLength) : NEQ(vecLength)
{
    	this->f = f;
    	this->userData = userData;
    	this->jtv = nullptr;
    	this->delta = nullptr;

    	krylov = new Kiops(MaxPhiOrder1 + 2, maxKrylovIters, tmpl, NEQ);//Originally MaxPhiOrder1+1
    	integratorStats = new IntegratorStats(NumProjectionsPerStep);

        realtype hMax=0;
    	fy = N_VClone(tmpl);
    	Y1 = N_VClone(tmpl);
    	fY1 = N_VClone(tmpl);
    	hfy = N_VClone(tmpl);
    	hb1 = N_VClone(tmpl);
    	hb2 = N_VClone(tmpl);
    	r1 = N_VClone(tmpl);
    	r2 = N_VClone(tmpl);
        r3=N_VClone(tmpl);
        Remainder = N_VClone(tmpl);
    	tmpVec = N_VClone(tmpl);
    	zeroVec = N_VClone(tmpl);
    	// Zero out zeroVec by using the trusted vector tmpl.
    	// If zeroVec has random data in it, it may have NAN values.
    	// zeroVec = 0.0 * zeroVec may accidentally be the calculation
    	// 0.0 * NAN = NAN, and zeroVec will not end up zero.
    	// Instead zeroVec = 0.0 * templ, since templ is trusted
    	// to have no NAN values in it.
    	N_VConst(0, fy);
    	N_VConst(0, Y1);
    	N_VConst(0, hfy);
    	N_VConst(0, hb1);
    	N_VConst(0, hb2);
    	N_VConst(0, r1);
    	N_VConst(0, zeroVec);
    	N_VConst(0, r2);
    	N_VConst(0, r3);

}

//==========
//Destructor
//==========
Epi3VChem_KIOPS::~Epi3VChem_KIOPS()
{
    	delete krylov;
    	delete integratorStats;

    	N_VDestroy(fy);
    	N_VDestroy(Y1);
    	N_VDestroy(fY1);
    	N_VDestroy(hfy);
    	N_VDestroy(hb1);
    	N_VDestroy(hb2);
    	N_VDestroy(r1);
    	N_VDestroy(r2);
    	N_VDestroy(r3);
    	N_VDestroy(Remainder);
    	N_VDestroy(tmpVec);
    	N_VDestroy(zeroVec);
        N_VDestroy(Scratch1);
}


IntegratorStats *Epi3VChem_KIOPS::Integrate(const realtype hStart, const realtype hMax, const realtype absTol,
                const realtype relTol, const realtype t0, const realtype tFinal,
                const int numBands, int basisSizes[], N_Vector y)
{
	realtype * 	data	= N_VGetArrayPointer(y);		//Query the state in debugger with this
	realtype fac		= pow(0.25, .5);
	realtype PercentDone	= 0;
	int PercentDots		= 0;
	int ProgressDots	= 0;
	myPb2 * pb{static_cast<myPb2 *> (userData)};    //Recast
	ofstream myfile;
	pb->MaxStepTaken		= 0;
	pb->MinStepTaken		= 1;
    myfile.open(pb->dumpJacFile, std:: ios_base::app);
	ofstream datafile;
	datafile.open("EPI3VData.txt", std::ios_base::app);
	int IgnDelay		= 0;
	//=======================
	//Standard error checking
	//=======================
	if( hStart < ZERO )
	{
		printf("Time step h is to small. \n");
		exit(EXIT_FAILURE);
	}
	if( tFinal < t0 )
	{
		printf("Starting time is larger the end time. \n");
		exit(EXIT_FAILURE);
	}
	realtype t = t0, hNew= hStart, h=hStart;                //Use this an initial guess
	int count 	= 0;
	int IntSteps	= 0;
        if(hMax < hStart)
        {
			hNew=hMax;
			h=hMax;
        }
        //End Error checking block.
        bool finalStep= false;                  //Set this necessary EOF flag
        N_Vector yTemp = N_VClone(y);
        N_VConst(1.0, yTemp);
        realtype sqrtN = EPICRSqrt(N_VDotProd(yTemp,yTemp));
        realtype krylovTol= 0.1 * sqrtN * absTol;
        //Main integration loop
        while(t<tFinal)
        {
			realtype Err=5.0;
			realtype ErrEst=0;
			count = 0;
			IntSteps ++;
			while(Err > 1 )//was 1
			{
				count ++;
				h=hNew;//h = k;
				// f is RHS function from the problem; y = u_n
				f(t, y, fy, userData); 							// f(t, y) = fy
				N_VScale(h, fy, hfy); 							//Scale f(y);
				SUPER_CHEM_JAC_TCHEM(t, y, y,     pb->Mat, this->userData, y,y,y);
				JTimesV jtimesv(jtv, f, delta, t, y, fy, userData, tmpVec);

				//Mayya's new method//
				// Y1= y_n + phi_1(6/8 hJ) h f(y_n)
				// R(z)= f(z)-f(y_n) - J*(z-y_n)
				// y(n+1)= y_n + phi_1(hJ) h f(y_n) + 2 phi_3(hj) h r(Y1)
				N_Vector stage1InputVecs[] = {zeroVec, hfy}; //Set the b vector
				N_Vector stage1OutputVecs[] = {r1,r2}; //Set output vectors
				N_Vector stage2OutputVecs[]= {r3};
				realtype timePts[] = {6.0/8.0, 1.0};
				realtype timePts2[]= {1.0};

				krylov->Compute(2, stage1InputVecs, timePts, 2, stage1OutputVecs, &jtimesv,
					h, krylovTol, basisSizes[0], &integratorStats->krylovStats[0]);
				//Set  r1= 6/8 phi_1(6/8 hJ), r2=phi_1(hJ)
				this->CheckNanBlowUp(r1, N_VGetLength(r1));
				this->CheckNanBlowUp(r2, N_VGetLength(r2));
				//std :: cout << "Second order cleared\n";

				N_VScale(8.0/6.0, r1, r1);              //Set r1=phi_1 (6/8 hJ) h f_n //PBFlag
				N_VLinearSum(1.0, y, 1.0, r1, Y1);      //Set Y1= y_n + phi_1(6/8 hJ) h f(y_n)= y_n + r1
				f(t,Y1, fY1, userData);                 //f(t, Y1)= fY1
				jtimesv.ComputeJv(r1,Remainder);        //Set Remainder = J (Y1-y_n)= J (r1)
				N_VLinearSum(-1.0, fy, -1.0, Remainder, Remainder);//Remainder = J(r1) - f(y_n)
				N_VLinearSum(1.0, Remainder, 1.0, fY1, Remainder);//Remainder = R(Y1) = J(r1) - f(y_n) + f(Y1)
				N_VScale(h,Remainder,Remainder);        //set Remainder= h R(Y1)
				N_Vector stage2InputVecs[]= {zeroVec, zeroVec, zeroVec, Remainder}; //[0,0,0,hR(Y1)]


				krylov->Compute(4, stage2InputVecs, timePts2, 1, stage2OutputVecs, &jtimesv,
						h, krylovTol, basisSizes[0], &integratorStats->krylovStats[0]);
				N_VScale(2.0, r3, r3);                       	//R3 is also the error estimate
				//Needs editing below
				N_VAbs(y, Scratch1);							//Scratch1 sp to be high order
				N_VScale(relTol, Scratch1, Scratch1);			//relTol*|y|->Scratch1
				N_VAddConst(Scratch1, absTol, Scratch1);		//relTol*|y|+absTol ->Scratch1
				N_VDiv(r3, Scratch1, Scratch1);					//ErrEst/(relTol*|y| + absTol)->Scratch1
				ErrEst = N_VDotProd(Scratch1, Scratch1);		//dot(ErrEst)
				ErrEst = ErrEst/ N_VGetLength(r3);				//Normalize ErrEst
				ErrEst = EPICRSqrt(ErrEst);						//sqrt ErrEst
				Err = ErrEst;                               	//Finalize Error Estimate

				//============================
				//Past this point errors arise
				//============================
				hNew = 0.9 * h * pow(ErrEst, -1.0/ 2.0);                  //Create New Step, usually .9
				if(hNew>100*h)//we increase the time
					hNew= 2*h;			//throttle
				if(1000*hNew<h)
					hNew= 0.01*h;		//throttle
				if( hNew>hMax)
					hNew=hMax;
				if( hNew <= ZERO)
				{
					//Perform Data dump
					printf("There is a possible singularity in the solution\n");
					std :: cout << "time stamp: " << t << std :: endl;
					std :: cout << "hNew: " << hNew<< std :: endl;
					std :: cout << "ErrEst: " << Err << std :: endl;
					std :: cout << "y: \n";
					N_VPrint_Serial(y);
					exit(EXIT_FAILURE);
				}
			}//Exit Adaptive loop
			//Step accepted, Do the following update
			N_VLinearSum(1.0, y, 1.0, r2, y); 		// Second Order = y + r2 = y + phi( ) * hfy
			N_VLinearSum(1.0 ,y, 1.0, r3, y); 		// Add the third order component.
			integratorStats->Step();                //Recompute integrator statistics

			this->Clean(y, N_VGetLength(y));		//Clean the data
			//this->CheckNanBlowUp(y, N_VGetLength(y));		//Check for errors
			t = t + h;                              		//Advance the time
			myPb2 * pb{static_cast<myPb2 *> (userData)};	//Recast

			//Check curret stepsize against mins and maxes
			if(h > pb->MaxStepTaken)
				pb->MaxStepTaken=h;
			if(h < pb->MinStepTaken)
				pb->MinStepTaken=h;

			if(data[0] >1500 && IgnDelay ==0)
			{
				cout << "Ign Delay between: ";
				cout << t-h << "and " << t << endl;
				IgnDelay=1;
				pb->ignTime= t - h/2.0;
			}
			pb->t=t;					//Set time for pass back;
			if(pb->Movie==1)
			{
				realtype * jacPtr	= N_VGetArrayPointer(pb->Jac);
				for(int i = 0 ; i < N_VGetLength(pb->Jac); i ++)
				{
					myfile << jacPtr[i];
					myfile.flush();
					myfile << "\n";
				}
				myfile << h << "\n";
				//Print y data to file
				for(int i = 0 ; i < N_VGetLength(y); i++)
					datafile << data[i] << endl;

				datafile << h << endl;
				datafile << Err << endl;
				datafile << IgnDelay << endl;
				datafile << t << endl;
				datafile << relTol << endl;
				datafile << absTol << endl;
			}
		ProgressDots = TrackProgress(tFinal, t, PercentDone, ProgressDots);
		//Check exit conditions
		if(finalStep)
			break;						//Yes?  	Exit
		if(t+hNew>=tFinal)				//Check Overstepping final time
		{
			hNew=tFinal-t;				//Set next step
			if(hNew < 1e-10)			//Close enough?
				break;					//Yes?  	Exit
			finalStep = true;			//No?  		One more step
		}
	}//Loop
	std :: cout << std :: endl;
	if(pb->Movie==1)
	{
		myfile.close();
		datafile.close();
	}
	return integratorStats;
}

//Cleaning function


void Epi3VChem_KIOPS::Clean(N_Vector y, int length)
{
        realtype * yD = N_VGetArrayPointer(y);
        for(int i = 0;   i<length;  i++)
                if(yD[i]<0)
		{
                        yD[i]=0;
		}

}


//is nan blowup check

void Epi3VChem_KIOPS::CheckNanBlowUp(N_Vector State, int vecLength)
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


//Does the matrix multiplication for the method
void Epi3VChem_KIOPS::MatMult(N_Vector R, N_Vector A, N_Vector B, int row, int col)
{
							//Remember this is row oriented data
	realtype * rPtr = N_VGetArrayPointer(R);
	realtype * aPtr	= N_VGetArrayPointer(A);
	realtype * bPtr	= N_VGetArrayPointer(B);
	int length	= N_VGetLength(A);
	int curCol	= 0;
	int curRow	= 0;

	for( int i = 0 ; i < length; i ++) 		//Over the entirety of the matrix
	{
		curRow = floor(i/row);			//Find your row location: floor(index/rowlength)
		curCol = i%col;				//Find your col:	  remainder(i/colLength)
		for( int j = 0; j < row; j ++)		//Work across the column (i,j)
			rPtr[i] +=aPtr[curRow*row + j]*	//R(i,j)= sum A(i,k)*B(k,j)
				bPtr[curCol + j*row];	//Set the result pointer
	}


}

void Epi3VChem_KIOPS::TestMatMult()
{
	int len		= 9;
	N_Vector R	= N_VNew_Serial(len);
	N_Vector A	= N_VClone(R);
	N_Vector B	= N_VClone(A);

	realtype * aPtr	= N_VGetArrayPointer(A);
	realtype * bPtr = N_VGetArrayPointer(B);
	realtype * rPtr	= N_VGetArrayPointer(R);

	for( int i = 0; i < len; i ++)
	{
		aPtr[i] = i;
		bPtr[i] = 1;
	}

	this->MatMult(R, A, B, 3, 3);
	std :: cout << "Printint test Vector\n";
	for( int i = 0; i < len; i++)
		std :: cout << rPtr[i] << "\t\t";
	std :: cout << "\n";
}

void Epi3VChem_KIOPS::Phi2(N_Vector Jac, realtype h, N_Vector Result)
{
	N_Vector Scratch= N_VClone(Jac);
	int len		= N_VGetLength(Jac);
	len 		= sqrt(len);
	this->MatMult(Scratch, Jac, Jac, len, len);		//Gets J^2
	N_VScale(1.0/6.0*h*h, Scratch, Scratch);		//Gets( hJ)^2/6
	N_VLinearSum(h/2.0, Jac, 1.0, Scratch, Scratch);	//Gets hJ/2 + (hJ)^2/6

	N_VDestroy_Serial(Scratch);

}

int Epi3VChem_KIOPS::TrackProgress(realtype FinalTime, realtype TNext, realtype PercentDone, int ProgressDots)
{
        PercentDone=floor(TNext/FinalTime*100);
        for (int k=0; k<PercentDone-ProgressDots;k++)
        {
                cout<<":";
                cout.flush();
        }
        return PercentDone;
}
