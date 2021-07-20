/*
 * ===========================================================================================
 * This is serial implementation.
 * ===========================================================================================
 */

#include <stdio.h>
#include <math.h>
#include "Epic.h"
#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TChem_CommandLineParser.hpp"
#include "TChem_Util.hpp"
#include "TChem_KineticModelData.hpp"
#include "TChem_Impl_IgnitionZeroD_Problem.hpp" // here is where Ignition Zero D problem is implemented
#include "TChem_Impl_NewtonSolver.hpp"
#include "TChem_Impl_TrBDF2.hpp"
#include <chrono>
#include "Epi3_KIOPS.h"


//#define TCHEMPB TChem::Impl::IgnitionZeroD_Problem<TChem::KineticModelConstData<Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> >>
#define TCHEMPB TChem::Impl::IgnitionZeroD_Problem	<TChem::KineticModelConstData	<Kokkos::Device	<Kokkos::Serial, Kokkos::HostSpace>  >	>
//Change to Serial, this can be changed by altering the TChem master build profile to include OPENMP on or off

//=====================
//Prototypes & Classes
//=====================
//Add an additional class elements to pass to epic.
class myPb : public TCHEMPB{
	public:
	//members
	ordinal_type  num_equations;
	N_Vector Jac;
	//fucntions
	ordinal_type	get_num_equations(void)
	{
		return this->num_equations;
	}
};

//Initial conditions
N_Vector InitialConditionsKappa(int, realtype);//number of equations and eps
N_Vector InitialConditionsH(int,int);
N_Vector InitialConditionsGri(int);//Need to add other samples later

//RHS functions
int RHS_KAPPA(realtype , N_Vector , N_Vector, void *);
int RHS_TCHEM(realtype, N_Vector, N_Vector, void *);

//JtV functions
int Jtv_KAPPA(N_Vector , N_Vector , realtype , N_Vector , N_Vector , void * , N_Vector);
int Jtv_TCHEM(N_Vector , N_Vector ,realtype, N_Vector, N_Vector , void* , N_Vector );

//Misc functions
int CheckStep(realtype, realtype);
void PrintFromPtr(realtype *,  int);
void ErrorCheck(ofstream &, N_Vector, realtype *, int, realtype);
void PrintDataToFile(ofstream &, realtype *,int, realtype);
//Creating the integrators
Epi2_KIOPS *	CreateEPI2Integrator(CVRhsFn, CVSpilsJacTimesVecFn, void *, int, N_Vector, const int, int);
Epi3_KIOPS *	CreateEPI3Integrator(CVRhsFn, CVSpilsJacTimesVecFn, void *, int, N_Vector, const int, int);
EpiRK4SC_KIOPS *    CreateEPIRK4SCIntegrator(CVRhsFn, CVSpilsJacTimesVecFn, void *, int, N_Vector, const int, int);

//Used in Jtv
void MatrixVectorProduct(int, N_Vector, realtype, N_Vector, realtype *);
//More misc Functions
void TestMatrixVectorProduct();
realtype LocalErrorEstimate(realtype, realtype, N_Vector, N_Vector, N_Vector, void *);

//=====================
//Namespaces and globals
//=====================
using namespace std;

int JacCnt=0;
int RHSCnt=0;

double JacTime=0;
double RHSTime=0;
double MatVecTime=0;
Kokkos::Impl::Timer g_timer;

//====================
//Main
//====================
int main(int argc, char* argv[])
{
	//====================
	// Intial Declarations
	//====================
	static realtype FinalTime = 0;//1.0e-4;//1e-4 seems to be the max
	static const int NumBands = 3;//Epic stuff, default is 3.
	realtype StepSize=0;
	realtype KTime=0;
	realtype PrintTime=0;
	int ProgressDots=0; //From 0 to 100; only care about percentage
	int StepCount = 0;
	realtype PercentDone=0;//void *userData=nullptr; //Empty problem pointer
	realtype TNow=0;
	realtype TNext=0;
	string MyFile="Default.txt";
	string Method="EPI2";//This is the default method
	static realtype KrylovTol=1e-14;
	int UseJac=1; //will we use the Jacobian or not
	int SampleNum=0;
	int Experiment=1; //default is hydrogen, set to 0 for kapila
	int Profiling=0;//default to no profiling
	int number_of_equations=0;

	//=====================================================
	//Kokkos Input Parser
	//=====================================================
	std::string chemFile("chem.inp");
 	std::string thermFile("therm.dat");
 	/// parse command line arguments --chemfile=user-chemfile.inp --thermfile=user-thermfile.dat
	/// --Stepsize=1e-8  --FinalTime=1e-2  --MyFile="Filename"  --KrylovTole=1e-14
	//  --UseJac=0 --SampleNum={1,2,3}, --Method="Integrator name"
	//  --Experiment={0,1,2}
  	/// with --help, the code list the available options.
 	TChem::CommandLineParser opts("This example computes reaction rates with a given state vector");
 	opts.set_option<std::string>("chemfile", "Chem file name e.g., chem.inp", &chemFile);
 	opts.set_option<std::string>("thermfile", "Therm file name e.g., therm.dat", &thermFile);
	opts.set_option<realtype>("StepSize", "StepSize desired", &StepSize);//New
	opts.set_option<realtype>("FinalTime", "The final simulation time", &FinalTime);//New
	opts.set_option<std::string>("MyFile","Where we output data", &MyFile);//New
	opts.set_option<realtype>("KrylovTol", "KrylovTolerance", &KrylovTol);
	opts.set_option<int>("UseJac","Will we use the Jacobian", &UseJac);
	opts.set_option<int>("SampleNum", "Sample used", &SampleNum);
	opts.set_option<std::string>("Method", "The time integration method", & Method);
	opts.set_option<int>("Experiment", "Experiment chosen", &Experiment);
	opts.set_option<int>("Profiling", "Output profiling data or not", &Profiling);

	const bool r_parse = opts.parse(argc, argv);
	if (r_parse)
		return 0; // print help return


	//======================================================
	//Check inputs and output to console for review.
	//======================================================
	ofstream myfile(MyFile, std::ios_base::app);
	cout<<"\n\n\n";//Give the program some space.
	int Steps=CheckStep(FinalTime, StepSize);//Checks the number of steps, Moved
	cout<<"=====================Sim Parameters=======================\n";
	cout<<setprecision(17)<<fixed;
	cout<<"Final Time: " << FinalTime <<"\t\t Step Size: " <<StepSize;
	cout<<"\t\t Krylov Tolerance: " <<KrylovTol<< "\t\t Writing to file: "<<MyFile <<endl;
	switch (Experiment)
	{
	case 0:
		cout<<"===============Kapila experiment selected================\n";
		break;
	case 1:
		cout <<"==============Hydrogen Experiment selected==============\n";
		break;
	case 2:
		cout <<"==================Gri Experiment selected===============\n";
		break;
	}
	//=====================================================
	//Kokkos block: working, pruning old code
	//=====================================================

	Kokkos::initialize(argc, argv);
	/// all Kokkos varialbes are reference counted objects. they are deallocated within this local scope.
	{//begin local scope
		//=====================================
		//Kokkos Sub-declarations
		//=====================================
		/// scalar type and ordinal type
		using real_type = double;
		using ordinal_type = int;

		/// Kokkos environments - host device type and multi dimensional arrays
		/// note that the 2d view use row major layout while most matrix format uses column major layout.
		/// to make the 2d view compatible with other codes, we need to transpose it.
		using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
		using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
		//using real_type_2d_view_type = Tines::value_type_2d_view<real_type,host_device_type>;

	    	/// construct TChem's kinect model and read reaction mechanism
    		TChem::KineticModelData kmd(chemFile, thermFile);

	    	/// construct const (read-only) data and move the data to device.
	    	/// for host device, it is soft copy (pointer assignmnet).
	    	auto kmcd = kmd.createConstData<host_device_type>();

    		/// use problem object for ignition zero D problem, interface to source term and jacobian
    		/// other available problem objects in TChem: PFR, CSTR, and CV ignition
    		using problem_type = TChem::Impl::IgnitionZeroD_Problem<decltype(kmcd)>;

	    	/// state vector - Temperature, Y_0, Y_1, ... Y_{n-1}), n is # of species.
		//Declare my variables locally

		//Set number of equations
		if(Experiment==0)
			number_of_equations=3;
		else
			number_of_equations=problem_type::getNumberOfEquations(kmcd);

		//Declare the state variable locally
		N_Vector y = N_VNew_Serial(number_of_equations); //The data y
		//Adding another N_Vector declaration here creates a bug.  Reason unknown.

		//Set Initial conditions based on experiment#
		switch(Experiment)
		{
		case 0:
			N_VScale(1., InitialConditionsKappa(number_of_equations,1e-2),y);
			break;
		case 1:
        		N_VScale(1., InitialConditionsH(number_of_equations,SampleNum), y);
        		break;
		case 2:
			N_VScale(1., InitialConditionsGri(number_of_equations), y);
			break;
		}

		//set the pointer to the state
		realtype *data = NV_DATA_S(y);


    		/// TChem does not allocate any workspace internally.workspace should be explicitly given from users.
	    	/// you can create the work space using NVector, std::vector or real_type_1d_view_type (Kokkos view)
	    	/// here we use kokkos view
		//Set up an ignition problem and the workspace.
		//Note: Kapila (Experiment 0) can run this, but it will be unused
		const ordinal_type problem_workspace_size = problem_type::getWorkSpaceSize(kmcd);
		real_type_1d_view_type work("workspace", problem_workspace_size);
	    	/// set problem
		myPb problem;
		void *pbptr= &problem;

      		/// initialize problem
		const real_type pressure(101325);//Constant pressure
		real_type_1d_view_type fac("fac", 2*number_of_equations);
      		problem._p = pressure; // pressure
		problem._fac= fac;
      		problem._work = work;  // problem workspace array
     		problem._kmcd = kmcd;  // kinetic model


		problem.num_equations=number_of_equations;
		problem.Jac=N_VNew_Serial(number_of_equations*number_of_equations);//Make the Jacobian
		const int MaxKrylovIters = 500;//500


		//==============================================
		//Create integrator
		//==============================================
		Epi2_KIOPS *integrator = NULL;
		Epi3_KIOPS *integrator2 = NULL;
		EpiRK4SC_KIOPS *integrator3 = NULL;

		//Parse the experiment cases
		if(Experiment!=0)
		{
			if(Method=="EPI2")
				integrator = CreateEPI2Integrator(RHS_TCHEM, Jtv_TCHEM, pbptr,
					MaxKrylovIters, y,number_of_equations, UseJac);
			else if(Method=="EPI3")
				integrator2 = CreateEPI3Integrator(RHS_TCHEM, Jtv_TCHEM, pbptr,
					MaxKrylovIters, y, number_of_equations, UseJac);
			else if(Method=="EPIRK4")
				integrator3 = CreateEPIRK4SCIntegrator(RHS_TCHEM, Jtv_TCHEM, pbptr,
					 MaxKrylovIters, y, number_of_equations, UseJac);
		}else{
                        if(Method=="EPI2")
                                integrator = CreateEPI2Integrator(RHS_KAPPA, Jtv_KAPPA, pbptr,
                                        MaxKrylovIters, y,number_of_equations, UseJac);
                        else if(Method=="EPI3")
                                integrator2 = CreateEPI3Integrator(RHS_KAPPA, Jtv_KAPPA, pbptr,
                                        MaxKrylovIters, y, number_of_equations, UseJac);
                        else if(Method=="EPIRK4")
                                integrator3 = CreateEPIRK4SCIntegrator(RHS_KAPPA, Jtv_KAPPA, pbptr,
                                         MaxKrylovIters, y, number_of_equations, UseJac);
                }


		//========================
        	// Set integrator parameters
        	//========================
        	int startingBasisSizes[] = {3, 3};

        	//=================================
        	// Run the time integrator loop
        	//=================================
        	cout<<"\t\tIntegrator progress:\n[";
		cout.flush();
		auto Begin=std::chrono::high_resolution_clock::now();
		while(StepCount<Steps)
        	//while(TNext<FinalTime)
        	{
                	TNow= StepCount*StepSize;
                	TNext=(StepCount+1)*StepSize;
			auto Start=std::chrono::high_resolution_clock::now();//Time integrator

			if(	Method=="EPI2")
			{
                		integrator->Integrate(StepSize, TNow, TNext, NumBands, y,
						KrylovTol, startingBasisSizes);
			}else if(Method == "EPI3")
			{
				integrator2->Integrate(StepSize, TNow, TNext, NumBands, y,
						KrylovTol, startingBasisSizes);
			}else if(Method =="EPIRK4")
			{
				integrator3->Integrate(StepSize, TNow, TNext, NumBands, y,
                                                 KrylovTol, startingBasisSizes);
			}

			//Clock the time spent in the integrator
			auto Stop=std::chrono::high_resolution_clock::now();
                        auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
			KTime+=Pass.count()/1e9;
			//integrator->Integrate(StepSize/10, StepSize, KrylovTol, KrylovTol, TNow, TNext, NumBands, startingBasisSizes, y);//EpiRKxSV
			//=======================================
			//Error checking
			//=======================================
			if(Experiment!=0) //Error checking invalid for Experiment 0.
				ErrorCheck(myfile, y, data, number_of_equations, TNext);

			//clean
                	for (int j=0; j<number_of_equations; j++)
             		{
				if(data[j]<0)
                        	        data[j]=0;
                	}

			//=======================
			//Track the progress
			//=======================
                	PercentDone=floor(TNext/FinalTime*100);
                	for (int k=0; k<PercentDone-ProgressDots;k++)
                	{
                	        cout<<":";
                	        cout.flush();
                	}

			/*//Use if we want to do timeseries plots
			if(ProgressDots!=PercentDone)
				PrintDataToFile(myfile, data, number_of_equations,TNext);
			*/
                	ProgressDots=PercentDone;
                	StepCount++;

        	}//End integration loop
		auto End=std::chrono::high_resolution_clock::now();
    		auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(End - Begin);
		double IntTime = elapsed.count()/1e9;
        	TNow=TNext;
        	cout << "]100%\n\n";
		//Delete integrators
		delete integrator;
		delete integrator2;
		delete integrator3;
		cout<<"====================Integration complete===============\n";
	        //=======================================
        	//Console Output
        	//=======================================
		if(Experiment!=0)
		{
		//Hydrogen layout, should also be valid for GRI 3.0
 		//Temp H2 O2 O OH H2O H HO2 H2O2
        		cout << "===========================Data========================\n";
        		cout << "Temp=" << data[0] << "\t\t H2=" << data[1] <<"\t\t O2=" << data[2]<<endl;
        		cout << "Total Mass Fractions: " <<N_VL1NormLocal(y)-data[0];
        		cout << "\t\t Mass Fraction error: "<<abs( N_VL1NormLocal(y)-data[0]-1.0)<<endl;
		}else
			cout << "y=" << data[0] << "\t\t z=" << data[1] <<"\t\t Temp=" <<data[2] <<endl;

		//Case agnostic output
        	cout << "Simulation ended at time: " << TNow << "\t\t Steps taken: " << StepCount <<endl;

		if (Profiling ==1){//Invalid for experiment 0
			cout << "=======================Performance data================\n";
			cout << "Integration loop CPU time: "<<IntTime<< " seconds" <<endl;
			cout << "Time Integrator CPU time: "<<KTime<<" seconds\n";
			cout << "Time integrator percentage of loop: " << KTime/IntTime*100 << "%\n";
			if (UseJac==1)
			{
				cout << "\tJacobian Calls: " << JacCnt;
				cout << "\t\tJacobian Time: " << JacTime << " seconds\n";
				cout << "\t\t" << JacTime*100/KTime<<"%\n";
				cout << "\t\tJac Time cost: " <<JacTime/JacCnt << " seconds\n";
				cout << "\tMatVec Time: " << MatVecTime << " seconds\n";
				cout << "\t\t" << MatVecTime*100/KTime <<"%\n\n";
			}
			cout << "\tRHS Function calls: " << RHSCnt << "\t\tRHS time: ";
			cout << RHSTime <<" seconds\n";
			cout << "\t\t" << RHSTime/KTime*100<< "%\n";
			cout << "\t\tRHS Time cost: " << RHSTime/RHSCnt << " seconds\n\n";
			if(UseJac==1)
				cout << "\tRelative cost ration Jac/RHS: ";
				cout << JacTime/RHSTime/JacCnt*RHSCnt << endl;

		}
		cout << "=====================Printing data to file=================\n";
                PrintDataToFile(myfile,data,number_of_equations,StepSize);//Only print here for conv studies
		myfile << "\t\t" << IntTime <<endl;
        	myfile.close();
        	cout << "======================Exiting without error=============\n";
  	}//end local kokkos scope.

 	/// Kokkos finalize checks any memory leak that are not properly deallocated.
  	Kokkos::finalize();
}



/*
//=========================================================
//  //====  ==   ==  ==     ==    _.==_    ========  _===_
//  ||      ||   ||  ||\\   ||  ./     \\     ||    //   \\
//  ||====  ||   ||  || \\  ||  ||            ||    \\___
//  ||      ||   ||  ||  \\ ||	~\            ||         \\
//  ||      \\___//  ||   \\||   ~=___-//     ||    \\___//
//=========================================================
*/


//===========================================
//Check the number of steps
//FinalTime
//StepSize
//===========================================
int CheckStep(realtype FinalTime, realtype StepSize)
{
	cout<<"======================Steps=========================\n";
	cout<<"Checking the proposed number of steps...\n";
	cout<<std::setprecision(17);
	cout<<FinalTime/StepSize << " steps proposed...";

        if(floor( FinalTime/StepSize ) == FinalTime/StepSize )
	{
                static const int Steps= FinalTime/StepSize;
                cout << " accepted...\n";
		return Steps;
        }else if (abs(round(FinalTime/StepSize))-FinalTime/StepSize <1e-6 ){
		static const int Steps= round (FinalTime/StepSize);
		cout << Steps << " steps has been approximated...\n";
		cout << "We use this modified Step size "<<fixed<<setprecision(9)<< (FinalTime/Steps) <<endl;
		return Steps;
	}else{
                cout<<"Cannot perform non-integer number of steps!!\n";
                exit(1);
        }
}

	//======================================================
	//Initial condition functions
	//======================================================
//====================================================
//Kapila problem
//====================================================
N_Vector InitialConditionsKappa(int number_of_equations, realtype EPS)
{
        N_Vector y0 = N_VNew_Serial(number_of_equations);
        realtype *data = NV_DATA_S(y0);
        data[0]=1.0-.5*EPS;
        data[1]=.5*EPS;
        data[2]=1.0;
        return y0;
}


//====================================================
//Hydrogen problem intial conditions
//====================================================
N_Vector InitialConditionsH(int number_of_equations, int SampleNum)
{
	//Temp H2 O2 O OH H2O H HO2 H2O2
	N_Vector y0 = N_VNew_Serial(number_of_equations);
        realtype *data = NV_DATA_S(y0);
        data[0]=900;//1032.0; //They put 1200 for some reason
	//First sample	: 900	2.7431213550e-01        7.2568786450e-01
	//Second sample	: 900	5.926665910131902193e-02 9.407333408986811030e-01
	//Third sample	: 900	1.438365693957185942e-01 8.561634306042813503e-01
	//Forth sample	: 1000	5.926665910131902193e-02 9.407333408986811030e-01
	switch(SampleNum)
	{
	case 1:
         	data[1]=2.7431213550e-01;
         	data[2]=1-data[1];
		cout<<"=======================Sample #1 Selected=====================\n";
		break;
	case 2:
		data[1]=5.926665910131902193e-02;
		data[2]=9.407333408986811030e-01;
		cout<<"=======================Sample #2 Selected=====================\n";
		break;
	case 3:
		data[1]=1.438365693957185942e-01;
		data[2]=8.561634306042813503e-01;
		cout<<"=======================Sample #3 Selected======================\n";
		break;
	case 4: //1.000000000000000000e+03 1.013250000000000000e+05 5.926665910131902193e-02 9.407333408986811030e-01
		//Run to 1e-5 time;
		data[0] = 1000.0;
		data[1] = 5.926665910131902193e-02;
		data[2] = 9.407333408986811030e-01;
		cout<<"=======================Sample #4 Selected======================\n";
		break;
	default:
		cout<<"=========================Unknown Sample========================\n";
		exit(EXIT_FAILURE);
	}
        return y0;
}



//=======================================================
//Gri Initial Conditions
//=======================================================
N_Vector InitialConditionsGri(int number_of_equations)
{
	//1.000000000000000000e+03 1.013250000000000000e+05
	// 5.483858845025559731e-02 2.187578062376045740e-01 7.137587863547695255e-01 1.264481895737025810e-02
	N_Vector y0= N_VNew_Serial(number_of_equations);
	realtype *data = NV_DATA_S(y0);
	data[0] = 1000.0;//use 1000 for base
	data[1] =  5.483858845025559731e-02;
	data[2] =  2.187578062376045740e-01;
	data[3] =  7.137587863547695255e-01;
	data[4] =  1.264481895737025810e-02;
	return y0;
}


	//====================================================
	//RHS functions
	//====================================================

/*
 * ===========================================================================================
 *
 * Function RHS
 *
 * If y' = f(t, y) then RHS is function f.
 * Function RHS is written is such way that it can handle arrays or arbitrarly length.
 * It is important that all data is in form of an array.
 *
 * Inputs:
 * t          time
 * u          input vector
 * udot       result
 * userData
 *
 * ===========================================================================================
 */

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


/*
//============================================================
//RHS that uses the TCHEM rhs functions
//============================================================
*/
int RHS_TCHEM(realtype t, N_Vector u, N_Vector udot, void * pb)
{
	//TCHEMPB *pbPtr{ static_cast<TCHEMPB*>(pb)} ;//recast the type here.
	myPb * pbPtr{static_cast<myPb *> (pb)};//Recast
	ordinal_type number_of_equations=pbPtr->get_num_equations();
	realtype *y= NV_DATA_S(u);
	realtype *dy=NV_DATA_S(udot);
	/// we use std vector mimicing users' interface
	using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
	using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
        // output rhs vector and x via a kokkos wrapper.
	real_type_1d_view_type x(y,	number_of_equations);
	real_type_1d_view_type f(dy,	number_of_equations);
        //=================================
        // Compute right hand side vector
        //=================================
	auto member =  Tines::HostSerialTeamMember();
	g_timer.reset();
	pbPtr->computeFunction(member, x ,f);
	RHSTime+=g_timer.seconds();
	RHSCnt++;
	return 0;
}

/*
/============================================================
//Estimate the local error for adaptive step sizes
//At each step we need to estimate the LTE
//Attempt to use a smaller step as a reference
//===========================================================
*/

realtype LocalErrorEstimate(realtype StepSize, realtype t, N_Vector y, N_Vector dy, N_Vector scrap, void * pb)
{
	realtype LTE =0;
	return LTE;

}





/*
//=============================================================
//Matrix Vector product
//
//
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



/*
 * ===========================================================================================
 *
 * Function Jtv
 *
 * This function computes Jacobian matrix times some vector v. In Epirk Jacobian is never
 * stored and it is computed every time we need it and every time we compute it we acctualy
 * compute it we compute its product with some vector v.
 * Function Jtv is written is such way that it can handle arrays or arbitrarly length.
 * It is important that all data is in form of an array.
 * It is used in JTimesv class.
 *
 * Inputs:
 * v          vector being multiplied with Jacobian Matrix
 * Jv         result
 * t          time
 * u          vector used ot compute Jacobian Matrix
 * fu         f(u), i.e., RHS of u
 * userData
 * tmp
 * ===========================================================================================
 */

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
||Jtv using TCHEM
||N_Vector v: 	the v in Jv
||N_Vector Jv:  Jv
||realtype t:  	time, needed by KIOPS
||N_Vector u:	the state vector
||N_Vector fu:	rhs evaluated at u
||void * pb:	pointer to my modified problem class
||N_Vector tmp:	temporary vector used here to grab rows of the Jacobian
\\=============================================================
*/
int Jtv_TCHEM(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void* pb, N_Vector tmp)
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
	g_timer.reset();
        pbPtr->computeJacobian(member, x ,J);
	JacTime+=g_timer.seconds();
	JacCnt++;
	//===================
	//Compute JV
	//===================

	g_timer.reset();
	MatrixVectorProduct(number_of_equations, JacD, v, tmp, JV);
	MatVecTime+=g_timer.seconds();
	return 0;
}




//==============================================
//Print from pointer
//realtype* ptr 	object pointer
//int num_eqs		number of equations
//==============================================
void PrintFromPtr(realtype * ptr, int num_eqs)
{
	cout<<endl;
	for(int i=0; i<num_eqs; i++)
	{
		cout<<ptr[i]<<endl;
	}
	cout<<endl;
}

//===================================================
//Print data to output file
//myfile:	pointer to the filestream
//data:		pointer to the data
//number_of_equations
//t:		time/Stepsize
//===================================================


void PrintDataToFile(ofstream & myfile, realtype * data, int number_of_equations, realtype t)
{
	for (int i=0; i<number_of_equations; i++)
	{
		myfile<<setprecision(20)<<fixed <<data[i] <<"\t\t";
	}
	myfile << "\t\t " << t;
	myfile.flush();

}

//=======================================================
//RHS Error Check
//realtype* y 		y Data ptr
//realtype* rhs		rhs Data ptr
//int num_eqs		number of equations in the system
//========================================================
void ErrorCheck(ofstream & myfile, N_Vector y, realtype * data, int number_of_equations, realtype TNext)
{
	realtype MassError=abs( N_VL1NormLocal(y)-data[0]-1.0);
	if(MassError>.1 || abs(data[0])>1e5 )
	{
		PrintDataToFile(myfile, data, number_of_equations, TNext);
                cout<<"\n===============!!!Critical error!!!================\n";
                                if(MassError>.1)
                                {
                                        cout<<"Severe mass loss detected. Check output file for details.\n";
                                        myfile<<"Mass fraction error : " << abs( N_VL1NormLocal(y)-data[0]-1.0)<< endl;
                                }
                                if(abs(data[0])>1e5)
                                        cout<<"Temperature instability detected. Check output file for details.\n";

                                cout<<"\nFailed after step: "<< TNext<<endl;
                                exit(EXIT_FAILURE);
	}

}



		//===================================================
		//	========  ==	==   ========	.====.
		//	   ||     ||\   ||	||     //    \\
		//	   ||     || \  ||	||     \\_____
		//	   ||     ||  \ ||	||     _     \\
		//	========  ||   \||	||     \\____//
		//===================================================

//===================================================================
//Create a pointer to an EPI2_KIOPS integrator
//RHS
//JtV
//MaxKrylovIters		Max Krylov Iterations
//y				the reference template vector
//number_of_equations
//UseJac			boolean if we use the Jacobian or not
//===================================================================
Epi2_KIOPS* CreateEPI2Integrator(CVRhsFn RHS, CVSpilsJacTimesVecFn JtV, void * pbptr, int MaxKrylovIters,
				N_Vector y, const int number_of_equations, int UseJac)
{
                switch(UseJac)
                {
		case 0:{
                //EpiRK4SC_KIOPS *integrator = new EpiRK4SC_KIOPS(RHS_TCHEM,
                //EpiRK4SV *integrator = new EpiRK4SV(RHS_TCHEM,
                        Epi2_KIOPS *integrator = new Epi2_KIOPS(RHS,
                                pbptr,
                                MaxKrylovIters,
                                y,
                                number_of_equations);

                        cout<<"=================EPI2 w/o JtV created====================\n";
                        return integrator;
			break;}
                case 1:{
                        Epi2_KIOPS *integrator = new Epi2_KIOPS(RHS,
                                JtV,
                                pbptr,
                                MaxKrylovIters,
                                y,
                                number_of_equations);//This line, NSE will cause issues later.
                        cout<<"===================EPI2 JtV created=====================\n";
			return integrator;
                        break;}
		default:
			cout<<"===============EPI2 Jacobian  Option Invalid============\n";
			exit(EXIT_FAILURE);
		}
	return 0;
}

//===================================================================
//Create a pointer to an EPI3_KIOPS integrator
//RHS
//JtV
//MaxKrylovIters                Max Krylov Iterations
//y                             the reference template vector
//number_of_equations
//UseJac                        boolean if we use the Jacobian or not
//===================================================================
Epi3_KIOPS* CreateEPI3Integrator(CVRhsFn RHS, CVSpilsJacTimesVecFn JtV, void * pbptr, int MaxKrylovIters, N_Vector y, const int number_of_equations, int UseJac)
{
                switch(UseJac)
                {
                case 0:{
                //EpiRK4SC_KIOPS *integrator = new EpiRK4SC_KIOPS(RHS_TCHEM,
                //EpiRK4SV *integrator = new EpiRK4SV(RHS_TCHEM,
                        Epi3_KIOPS *integrator = new Epi3_KIOPS(RHS,
                                pbptr,
                                MaxKrylovIters,
                                y,
                                number_of_equations);

                        cout<<"=================EPI3 w/o JtV created====================\n";
                        return integrator;
                        break;}
                case 1:{
                        Epi3_KIOPS *integrator = new Epi3_KIOPS(RHS,
                                JtV,
                                pbptr,
                                MaxKrylovIters,
                                y,
                                number_of_equations);//This line, NSE will cause issues later.
                        cout<<"===================EPI3 JtV created=====================\n";
                        return integrator;
                        break;}
                default:
                        cout<<"===============EPI3 Jacobian  Option Invalid============\n";
                        exit(EXIT_FAILURE);
                }
        return 0;
}


//===================================================================
//Create a pointer to an EPIRK4SC_KIOPS integrator
//RHS
//JtV
//MaxKrylovIters                Max Krylov Iterations
//y                             the reference template vector
//number_of_equations
//UseJac                        boolean if we use the Jacobian or not
//===================================================================
EpiRK4SC_KIOPS * CreateEPIRK4SCIntegrator(CVRhsFn RHS, CVSpilsJacTimesVecFn JtV, void * pbptr, int MaxKrylovIters, N_Vector y, const int number_of_equations, int UseJac)
{
                switch(UseJac)
                {
                case 0:{
                //EpiRK4SC_KIOPS *integrator = new EpiRK4SC_KIOPS(RHS_TCHEM,
                //EpiRK4SV *integrator = new EpiRK4SV(RHS_TCHEM,
                        EpiRK4SC_KIOPS *integrator = new EpiRK4SC_KIOPS(RHS,
                                pbptr,
                                MaxKrylovIters,
                                y,
                                number_of_equations);

                        cout<<"=================EPIRK4SC_K w/o JtV created====================\n";
                        return integrator;
                        break;}
                case 1:{
                        EpiRK4SC_KIOPS *integrator = new EpiRK4SC_KIOPS(RHS,
                                JtV,
                                pbptr,
                                MaxKrylovIters,
                                y,
                                number_of_equations);//This line, NSE will cause issues later.
                        cout<<"===================EPIRK4SC_K JtV created=====================\n";
                        return integrator;
                        break;}
                default:
                        cout<<"===============EPIRK4SC_K Jacobian  Option Invalid============\n";
                        exit(EXIT_FAILURE);
                }
        return 0;
}



//========================================================
//This function runs simple tests to see if the matrix
//vector product function is correct or not.
//Appears to be working for the three test cases.
//========================================================
void TestMatrixVectorProduct()
{
	//initialize test
	int num_equations=2;
	N_Vector Test_X=N_VNew_Serial(num_equations);
	N_Vector Test_Jac=N_VNew_Serial(num_equations*num_equations);
	N_Vector tmp=N_VNew_Serial(num_equations);
	N_Vector Jv=N_VNew_Serial(num_equations);
	realtype * JV=NV_DATA_S(Jv);
	realtype * JacD=NV_DATA_S(Test_Jac);
	realtype * XPtr= NV_DATA_S(Test_X);
	int TestCase=3;
	//set test case;
	switch(TestCase)
	{
	case 1: 	// [  0, 1  // 1, 0 ]  x [1 // 0] = [0 //1 ]
		JacD[1]=1;
		JacD[2]=1;
		XPtr[0]=1;
		break;
	case 2: 	//[1 , 1 // 1,1] x [1,1] = [2// 2]
		JacD[0]=1;
		JacD[1]=1;
		JacD[2]=1;
		JacD[3]=1;
		XPtr[0]=1;
		XPtr[1]=1;
		break;
	case 3:		//[1, -1 \\ -1 ,1] x [1 // 1] = [0 // 0]
		JacD[0]=1;
		JacD[1]=-1;
		JacD[2]=-1;
		JacD[3]=1;
		XPtr[0]=1;
		XPtr[1]=1;
		break;
	default:
		cout<<"Test not run..."<<endl;
		break;
	}
	MatrixVectorProduct(num_equations, JacD, Test_X, tmp, JV);
	PrintFromPtr(JV, num_equations);

}
