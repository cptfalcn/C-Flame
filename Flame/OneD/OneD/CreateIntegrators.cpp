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

//=======================
//Create CVODE BDF Direct
//This is depreciated and no longer works.
//=======================
void * CreateCVODE(CVRhsFn RHS, CVSpilsJacTimesVecFn JtV, void * pbptr, SUNMatrix A, SUNLinearSolver LS,
			int num_eqs, N_Vector State, realtype relTol, realtype absTol, realtype StepSize,
			int UseJac)
{
	void * cvode_mem;
	int retVal = 0;
	realtype tret = 0;
	N_Vector AbsTol= N_VNew_Serial(num_eqs);
 	for ( int i = 0 ; i < num_eqs ; i++)
		NV_Ith_S(AbsTol,i)=absTol;

	cvode_mem = CVodeCreate (CV_BDF);
	//Set intial state
	retVal = CVodeInit(cvode_mem, RHS, 0, State);
	//Give CVODE the user problem
	retVal = CVodeSetUserData(cvode_mem, pbptr);
	//Set intial stepsize guess
	retVal = CVodeSetInitStep(cvode_mem, StepSize);//Modify this line
	//Set max stepsize
	retVal = CVodeSetMaxStep(cvode_mem, StepSize);
	//Set tolerances
	retVal = CVodeSVtolerances(cvode_mem, 1e-8, AbsTol);
	//Set linear solver
	retVal = CVodeSetLinearSolver(cvode_mem, LS, A);
	retVal = CVodeSetJacTimes(cvode_mem, NULL, JtV);
	if (UseJac==1)
   	{
        	retVal = CVodeSetJacFn(cvode_mem, CHEM_COMP_JAC_CVODE);//likely bugged AF
       		retVal = CVodeSetMaxNumSteps(cvode_mem, max(500,num_eqs*num_eqs*num_eqs));
      	}
	N_VDestroy_Serial(AbsTol);
	return cvode_mem;
}

//=======================
//Create CVODE BDF Krylov
//=======================

void * CreateCVODEKrylov(CVRhsFn RHS, CVSpilsJacTimesVecFn JtV, void * pbptr, SUNMatrix A, SUNLinearSolver LS,
                        SUNNonlinearSolver NLS , int num_eqs, N_Vector State, realtype relTol,
			realtype absTol, realtype StepSize, int UseJac)
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
	retVal = CVodeSetJacTimes(cvode_mem, NULL, JtV);//Edit in and out
        if (UseJac==1)
        {
                //retVal = CVodeSetJacFn(cvode_mem, CHEM_COMP_JAC);
		retVal = CVodeSetJacFn(cvode_mem, CHEM_COMP_JAC_CVODE);//likely bugged AF
                retVal = CVodeSetMaxNumSteps(cvode_mem, max(500,num_eqs*num_eqs*num_eqs));
        }
        N_VDestroy_Serial(AbsTol);
        return cvode_mem;
}

//========================
//Integrate
//========================
IntegratorStats* Integrate(int UseJac, realtype StepSize, realtype TNow, realtype TNext, int NumBands,
				N_Vector y, realtype KrylovTol, int startingBasisSizes[], void * pbptr,
				void * cvode_mem, string Method, Epi2_KIOPS * integrator,
				Epi3_KIOPS * integrator2, EpiRK4SC_KIOPS * integrator3,
				IntegratorStats *  IntStats)
{
	if(UseJac==1 && Method != "CVODE" && Method != "CVODEKry" && Method !="CVODEBiCGS")
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
     	}else if(Method =="CVODE"||Method=="CVODEKry"||Method=="CVODEBiCGS")
	{
		CVode(cvode_mem, TNow, y, &TNext, CV_NORMAL);
		//CVode(cvode_mem, TNow, y, &TNext, CV_ONE_STEP);
	}
	return IntStats;
}

void *		SelectIntegrator(CVRhsFn RHS, CVSpilsJacTimesVecFn JtV, void * pbptr, int MaxKrylovIters,
                                	N_Vector y, const int num_eqs, int UseJac, string Method,
					SUNMatrix A, SUNLinearSolver LS, realtype absTol,
					realtype relTol, realtype StepSize)
{
		unique_ptr<Epi2_KIOPS> DUMMY;
		if(Method=="EPI2")
		{
		unique_ptr<Epi2_KIOPS> EPI2_PTR;
		EPI2_PTR.reset( CreateEPI2Integrator(RHS, JtV, pbptr, MaxKrylovIters, y, num_eqs, UseJac));
		return reinterpret_cast<void *> (EPI2_PTR.get() );
		}
		else if(Method == "EPI3")
		{
		unique_ptr<Epi3_KIOPS> EPI3_PTR;
		EPI3_PTR.reset( CreateEPI3Integrator(RHS, JtV, pbptr, MaxKrylovIters, y, num_eqs, UseJac));
		return reinterpret_cast<void *> (EPI3_PTR.get() );
		}
		else if(Method == "EPIRK4")
		{
		unique_ptr<EpiRK4SC_KIOPS> EPIRK4_PTR;
		EPIRK4_PTR.reset( CreateEPIRK4SCIntegrator(RHS, JtV, pbptr, MaxKrylovIters, y, num_eqs,
					UseJac));
                return reinterpret_cast<void *> (EPIRK4_PTR.get() );
		}else//create CVODE
		{
		return CreateCVODE(RHS, JtV, pbptr, A, LS, num_eqs, y, relTol, absTol, StepSize, UseJac);
		}
}
