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
#include <omp.h>

#define TCHEMPB TChem::Impl::IgnitionZeroD_Problem<TChem::KineticModelConstData<Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> >>





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

N_Vector InitialConditionsH(int);
N_Vector InitialConditionsGri(int);
int RHS_TCHEM(realtype, N_Vector, N_Vector, void *);
int Jtv_TCHEM(N_Vector , N_Vector ,realtype, N_Vector, N_Vector , void* , N_Vector );
int CheckStep(realtype, realtype);
void PrintFromPtr(realtype *,  int);
void ErrorCheck(realtype *, int, string);
void PrintDataToFile(ofstream &, realtype *,int, realtype);


//=====================
//Namespaces
//=====================
using namespace std;



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
	int ProgressDots=0; //From 0 to 100; only care about percentage
	int StepCount = 0;
	realtype PercentDone=0;//void *userData=nullptr; //Empty problem pointer
	realtype TNow=0;
	realtype TNext=0;
	string MyFile="Default.txt";
	static realtype KrylovTol=1e-14;


	//=====================================================
	//Kokkos Input block
	//=====================================================
	std::string chemFile("chem.inp");
 	std::string thermFile("therm.dat");
 	/// parse command line arguments --chemfile=user-chemfile.inp --thermfile=user-thermfile.dat
	/// --Stepsize=1e-8  --FinalTime=1e-2  --MyFile="Filename"  --KrylovTole=1e-14
  	/// with --help, the code list the available options.
 	TChem::CommandLineParser opts("This example computes reaction rates with a given state vector");
 	opts.set_option<std::string>("chemfile", "Chem file name e.g., chem.inp", &chemFile);
 	opts.set_option<std::string>("thermfile", "Therm file name e.g., therm.dat", &thermFile);
	opts.set_option<realtype>("StepSize", "StepSize desired", &StepSize);//New
	opts.set_option<realtype>("FinalTime", "The final simulation time", &FinalTime);//New
	opts.set_option<std::string>("MyFile","Where we output data", &MyFile);//New
	opts.set_option<realtype>("KrylovTol", "KrylovTolerance", &KrylovTol);
	const bool r_parse = opts.parse(argc, argv);
	if (r_parse)
		return 0; // print help return


	//======================================================
	//Check inputs and output to console for review.
	//======================================================
	ofstream myfile(MyFile, std::ios_base::app);
	int Steps=CheckStep(FinalTime, StepSize);//Checks the number of steps, Moved
	cout<<"=====================Sim Parameters=======================\n";
	cout<<setprecision(17)<<fixed;
	cout<<"Final Time: " << FinalTime <<"\t\t Step Size: " <<StepSize;
	cout<<"\t\t Krylov Tolerance: " <<KrylovTol<< "\t\t Writing to file: "<<MyFile << endl;


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
    		const ordinal_type number_of_equations = problem_type::getNumberOfEquations(kmcd);
		N_Vector y = N_VNew_Serial(number_of_equations); //The data y
        	N_VScale(1., InitialConditionsH(number_of_equations), y);
		N_Vector yOld= N_VNew_Serial(number_of_equations);
		N_Vector yDiff=N_VNew_Serial(number_of_equations);
        	realtype *data = NV_DATA_S(y);
		realtype *dataOld=NV_DATA_S(yOld);


    		/// TChem does not allocate any workspace internally.workspace should be explicitly given from users.
	    	/// you can create the work space using NVector, std::vector or real_type_1d_view_type (Kokkos view)
	    	/// here we use kokkos view
		const ordinal_type problem_workspace_size = problem_type::getWorkSpaceSize(kmcd);
		real_type_1d_view_type work("workspace", problem_workspace_size);
	    	/// set problem
		myPb problem;
		void *pbptr= &problem;

      		/// initialize problem
		const real_type pressure(101325);//Constant pressure
      		problem._p = pressure; // pressure
      		problem._work = work;  // problem workspace array
     		problem._kmcd = kmcd;  // kinetic model
		problem.num_equations=number_of_equations;
		problem.Jac=N_VNew_Serial(number_of_equations*number_of_equations);//Make the Jacobian
		const int MaxKrylovIters = 500;//500


		//==============================================
		//Create integrator
		//==============================================
		//TChem::Example::TestProblemTrBDF problem;
		//EpiRK4SC_KIOPS *integrator = new EpiRK4SC_KIOPS(RHS_TCHEM,
		//EpiRK4SV *integrator = new EpiRK4SV(RHS_TCHEM,
		Epi2_KIOPS *integrator = new Epi2_KIOPS(RHS_TCHEM,
				//Jtv_TCHEM,//Turn this on and off
				pbptr,
				MaxKrylovIters,
				y,
				number_of_equations);//This line, NSE will cause issues later.
		cout<<"=====================Integrator created====================\n";
       		//========================
        	// Set integrator parameters
        	//========================
        	//const realtype KrylovTol = RCONST(1.0e-10);//1e-14
        	int startingBasisSizes[] = {1, number_of_equations};

        	//=================================
        	// Run the time integrator loop
        	//=================================
        	cout<<"Integrator progress:[";
		cout.flush();
		while(StepCount<Steps)
        	//while(TNext<FinalTime)
        	{
                	TNow= StepCount*StepSize;
                	TNext=(StepCount+1)*StepSize;
                	//Integrate
                	integrator->Integrate(StepSize,TNow, TNext, NumBands, y, KrylovTol, startingBasisSizes);
			//integrator->Integrate(StepSize/10, StepSize, KrylovTol, KrylovTol, TNow, TNext, NumBands, startingBasisSizes, y);//EpiRKxSV

			//=======================================
			//Error checking
			//=======================================
			if(abs( N_VL1NormLocal(y)-data[0]-1.0)>.1 || abs(data[0])>1e5 )
			{
				PrintDataToFile(myfile, data, number_of_equations, TNext);
				cout<<"\n===============!!!Critical error!!!================\n";
				if(abs( N_VL1NormLocal(y)-data[0]-1.0)>.1)
				{
					cout<<"Severe mass loss detected. Check TChemHydrogenData.txt for details\n";
					myfile<<"Mass fraction error : " << abs( N_VL1NormLocal(y)-data[0]-1.0)<< endl;
				}
				if(abs(data[0])>1e5)
					cout<<"Temperature instability detected. Check TChemHydrogenData.txt for details \n";

				cout<<"\nFinal time step: "<< TNow;
				cout<<"\t\t We completed " << PercentDone <<"%\n";
				exit(EXIT_FAILURE);
			}



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
        	TNow=TNext;
        	cout << "]100%\n\n";
		delete integrator;
		//If doing convergence studies only print the final data.
		cout<<"====================Integration complete===============\n";
	        //=======================================
        	//Console Output
        	//=======================================
        	//Temp H2 O2 O OH H2O H HO2 H2O2
        	cout<<"===========================Data========================\n";
        	cout <<"Temp=" << data[0] << "\t\t H2=" << data[1] <<"\t\t O2=" << data[2]<<endl;
        	cout << "Total Mass Fractions: " <<N_VL1NormLocal(y)-data[0];
        	cout << "\t\t Mass Fraction error: "<<abs( N_VL1NormLocal(y)-data[0]-1.0)<<endl;
        	cout << "Simulation ended at time: " << TNow;
        	cout << "\t\t Steps taken: " << StepCount <<endl;
                PrintDataToFile(myfile,data,number_of_equations,StepSize);

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



//====================================================
//Hydrogen problem intial conditions
//====================================================
N_Vector InitialConditionsH(int number_of_equations)
{
	//Temp H2 O2 O OH H2O H HO2 H2O2
         N_Vector y0 = N_VNew_Serial(number_of_equations);
         realtype *data = NV_DATA_S(y0);
         data[0]=900;//1032.0; //They put 1200 for some reason
	 /*2.7431213550e-01        7.2568786450e-01*/
         data[1]=2.7431213550e-01;
         data[2]=1-data[1];
         return y0;
}


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
	pbPtr->computeFunction(member, x ,f);
	//===============================
	//Check for NaN
	//===============================
	//ErrorCheck(y,number_of_equations, "y");
	//ErrorCheck(dy, number_of_equations, "rhs");
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
        pbPtr->computeJacobian(member, x ,J);

	//===================
	//Compute JV
	//===================
	//#pragma omp parallel for private(TMP)
	for (int i=0; i< number_of_equations; i++)//for each state
	{

		//#pragma omp parallel for private(JV)
		for(int j=0; j<number_of_equations; j++)//marches across the column
		{
			TMP[j]=JacD[j+i*number_of_equations];//Stored row major
			if (j+1==number_of_equations)
				JV[i]=N_VDotProd(tmp,v);
		}
		//JV[i]=N_VDotProd(tmp,v);
	}
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
	myfile << "\t\t " << t << "\n";
	myfile.flush();

}

//=======================================================
//RHS Error Check
//realtype* y 		y Data ptr
//realtype* rhs		rhs Data ptr
//int num_eqs		number of equations in the system
//========================================================

void ErrorCheck(realtype * data, int num_eqs, string name)
{
	for (int i=0; i<num_eqs; i++)
        {
		if (isnan(data[i])==1)
		{
			cout<< endl << "NaN in "<< name <<" detected\n";
			PrintFromPtr(data,num_eqs);
                       	exit(EXIT_FAILURE);
		}
	}

}
