/*
 * ===========================================================================================
 * 
 * In this example we show how Epic package can be used with test problems which are not 
 * written using objects. 
 * 
 * -------------------------------------------------------------------------------------------
 * -------------------------------------------------------------------------------------------
 * 
 * This is serial implementation. 
 * 
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

#define NEQ 3 //number of equations
#define NSE 9 //Number of species, 9 for Hydrogen
#define EPS 1e-2
#define TCHEMPB TChem::Impl::IgnitionZeroD_Problem<TChem::KineticModelConstData<Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> >>

//====================
//Prototypes
//====================
//auto kmcd = kmd.createConstData<host_device_type>();
//using problem_type = TChem::Impl::IgnitionZeroD_Problem<decltype(kmcd)>;
//template<typename KineticModelConstDataType>;

N_Vector InitialConditions();
int RHS(realtype, N_Vector, N_Vector, void *);
int RHS_TCHEM(TCHEMPB , realtype, N_Vector, N_Vector, void *, ordinal_type);
int Jtv(N_Vector, N_Vector, realtype, N_Vector, N_Vector , void *, N_Vector);
int CheckStep(realtype, realtype);

//=====================
//Namespaces
//=====================


using namespace std;

int main(int argc, char* argv[])
{
	//====================
	// Declarations
	//====================
	static const realtype FinalTime = 10.0;
	static const int NumBands = 3;//Epic stuff
	realtype StepSize=1e-2;
	CheckStep(FinalTime, StepSize);//Checks the number of steps
	int ProgressDots=0; //From 0 to 100; only care about percentage
	int StepCount = 0;
	realtype PercentDone=0;
	void *userData=nullptr; //Empty problem pointer
	realtype TNow=0;
	realtype TNext=0;
	N_Vector y = N_VNew_Serial(NEQ); //The data y
	N_VScale(1., InitialConditions(), y);
	realtype *data = NV_DATA_S(y);
	bool LastStepAdjusted = 0;

	//=====================================
	//This block is for testing.
	//Reactivate at your own risk,
	//The entire block is a mess.
	//=====================================
	//N_Vector Is = N_VNew_Serial(NEQ);
	//N_Vector Jv = N_VNew_Serial(NEQ);
	N_Vector rhs= N_VNew_Serial(NSE);
	//N_Vector y0 = N_VNew_Serial(NEQ);
        //N_Vector y0 = N_VNew_Serial(NEQ);
	//N_VScale(1., InitialConditions(), y0);
        //realtype *Y0= NV_DATA_S(y0);
	//realtype *IS = NV_DATA_S(Is);
	//realtype *JV = NV_DATA_S(Jv);
	//realtype *RHSy=NV_DATA_S(rhs);
	//realtype *Y0= NV_DATA_S(y0);
	//N_VScale(1., InitialConditions(), y0);
	//for (int i=0; i<NEQ; i++){
	//	IS[i]=1.0;}
	//RHS(0, y, rhs, userData);
	//Jtv(Is, Jv, 1, y , rhs, userData, Is);
	//cout << "=================Initial Data===================\n";
	//cout << "=========================y======================\n";
	//cout << "y[0]=" << data[0] <<"\t\t y[1]=" << data[1] << "\t\t y[2]=" << data[2] << endl;
	//cout << "======================rhs========================\n";
	//cout << "rhs[0]=" << RHSy[0] << "\t\t rhs[1]=" << RHSy[1] << "\t\t rhs[2]="<<RHSy[2] <<endl;
	//cout << "========================JtV=====================\n";
	//cout <<  "Jtv[0] =" <<JV[0] <<"\t\t Jtv[1]="<< JV[1] <<"\t\t Jtv[2]=" << JV[2]  <<endl;
	//end test

	//=====================================================
	//Kokkos Input block: Working
	//=====================================================
	/**/
	std::string chemFile("chem.inp");
 	std::string thermFile("therm.dat");
        std::cout<<"Starting\n";
 	/// parse command line arguments --chemfile=user-chemfile.inp --thermfile=user-thermfile.dat
  	/// with --help, the code list the available options.
 	TChem::CommandLineParser opts("This example computes reaction rates with a given state vector");
 	opts.set_option<std::string>("chemfile", "Chem file name e.g., chem.inp", &chemFile);
 	opts.set_option<std::string>("thermfile", "Therm file name e.g., therm.dat", &thermFile);
	const bool r_parse = opts.parse(argc, argv);
	if (r_parse)
	return 0; // print help return
	cout<< "=================Chemistry files==========================\n";
	cout<< chemFile <<"\t" << thermFile << endl;
	/**/
	//=====================================================
	//Testing Kokkos block: working
	//=====================================================

	Kokkos::initialize(argc, argv);
	{
	///
    	/// 3. Type definitions
	///

	/// scalar type and ordinal type
		using real_type = double;
		using ordinal_type = int;

	/// Kokkos environments - host device type and multi dimensional arrays
	/// note that the 2d view use row major layout while most matrix format uses column major layout.
	/// to make the 2d view compatible with other codes, we need to transpose it.
		using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
		using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
		using real_type_2d_view_type = Tines::value_type_2d_view<real_type,host_device_type>;

       	//
	//4. Construction of a kinetic model
	//

    	/// construct TChem's kinect model and read reaction mechanism
    		TChem::KineticModelData kmd(chemFile, thermFile);

    	/// construct const (read-only) data and move the data to device.
    	/// for host device, it is soft copy (pointer assignmnet).
	    	auto kmcd = kmd.createConstData<host_device_type>();

    	///
    	/// 5. Problem setup
    	///

    	/// use problem object for ignition zero D problem providing interface to source term and jacobian
    	/// other available problem objects in TChem: PFR, CSTR, and CV ignition
    		using problem_type = TChem::Impl::IgnitionZeroD_Problem<decltype(kmcd)>;

    	/// state vector - Temperature, Y_0, Y_1, ... Y_{n-1}), n is # of species.
    		const ordinal_type number_of_equations = problem_type::getNumberOfEquations(kmcd);

    	///
    	/// TChem does not allocate any workspace internally.
    	/// workspace should be explicitly given from users.
    	/// you can create the work space using NVector, std::vector or real_type_1d_view_type (Kokkos view)
    	/// here we use kokkos view
		const ordinal_type problem_workspace_size = problem_type::getWorkSpaceSize(kmcd);
		real_type_1d_view_type work("workspace", problem_workspace_size);

	///
    	/// 6. Input state vector
    	///

    	/// set input vector
    		const real_type pressure(101325);

    	/// we use std vector mimicing users' interface
    		std::vector<real_type> x_std(number_of_equations, real_type(0));
    		x_std[0] = 1200; // temperature
    		x_std[1] = 0.1;  // mass fraction species 1
    		x_std[2] = 0.9;  // mass fraction species 2

    	/// output rhs vector and J matrix
    		std::vector<real_type> rhs_std(number_of_equations);
    		std::vector<real_type> J_std(number_of_equations * number_of_equations);

    	/// get points for std vectors
    		real_type * ptr_x = x_std.data();
    		real_type * ptr_rhs = rhs_std.data();
    		real_type * ptr_J = J_std.data();

    	///
    	/// 7. Compute right hand side vector and Jacobian
    	///
		//{
      			problem_type problem;

      	/// initialize problem
      			problem._p = pressure; // pressure
      			problem._work = work;  // problem workspace array
     			problem._kmcd = kmcd;  // kinetic model

      	/// create view wrapping the pointers
      			//real_type_1d_view_type x(ptr_x,   number_of_equations);
      			//real_type_1d_view_type f(ptr_rhs, number_of_equations);
      			//real_type_2d_view_type J(ptr_J,   number_of_equations, number_of_equations);

      	/// create a fake team member
      			const auto member = Tines::HostSerialTeamMember();

      	/// compute rhs and Jacobian
	/*
			//=====================================
 			//Here the  function for  calculating the right hand side goes
	  		//=====================================
      			problem.computeFunction(member, x, f);
      			problem.computeJacobian(member, x, J);

      	/// change the layout from row major to column major
      			for (ordinal_type j=0;j<number_of_equations;++j)
        			for (ordinal_type i=0;i<j;++i)
                			std::swap(J(i,j), J(j,i));
	*/
		//}
	///
    	/// 8. Check the rhs and Jacobian
    	///
	/*
	printf("RHS std vector \n" );
    		for (ordinal_type i=0;i<number_of_equations;++i) 
		{
      			printf("%e\n", rhs_std[i]);
        		std::cout<<"\n";
    		}

		printf("Jacobian std vector \n" );
	    	for (ordinal_type i=0;i<number_of_equations;++i)
		{
      			for (ordinal_type j=0;j<number_of_equations;++j)
			{
        			printf("%e ", J_std[i+j*number_of_equations] );
      			}
      			printf("\n" );
   	 	}

	*/
    	/// all Kokkos varialbes are reference counted objects. they are deallocated
    	/// within this local scope.
	//Try EPIC here.
	 const int MaxKrylovIters = 500;
        /* 
	Epi2_KIOPS *integrator = new Epi2_KIOPS(RHS_TCHEM,
                                       //Jtv,
                                       userData,
                                       MaxKrylovIters,
                                       y,
                                       NEQ);
	delete integrator;
	*/
	RHS_TCHEM(problem,1,y,rhs,userData, number_of_equations);
  	}

 	 /// Kokkos finalize checks any memory leak that are not properly deallocated.
  	Kokkos::finalize();


	//===================================================
	// Create the integrator
	// In this case we use EpiRK5C, for other integrators just use their name.
    	//===================================================
	const int MaxKrylovIters = 500;
    	Epi2_KIOPS *integrator = new Epi2_KIOPS(RHS,
                                      Jtv,
                                      userData,
                                      MaxKrylovIters,
                                      y,
                                      NEQ);

	//========================
	// Set integrator parameters
	//========================
	const realtype KrylovTol = RCONST(1.0e-14);//1e-14
	int startingBasisSizes[] = {3, 3};



	//=================================
    	// Run the time integrator loop
	//=================================
	cout<<"Integrator progress:[";
	while(TNext<FinalTime)
	{
		TNow= StepCount*StepSize;
		TNext=(StepCount+1)*StepSize;
		if(TNext>FinalTime&&(FinalTime != TNow))
		{
			StepSize=FinalTime-TNow;
			LastStepAdjusted=1;
		}
		//Integrate
		integrator->Integrate(
        	StepSize,TNow, TNext, NumBands, y,
		KrylovTol, startingBasisSizes);
		//If the data goes negative, which is not physical, correct it.
		//Removing this loop will cause stalls.
		for (int j=0; j<NEQ; j++)
		{
			if(data[j]<0)
			{
				data[j]=0;
			}
		}
		//Track the progress
		PercentDone=floor(TNext/FinalTime*100);
		for (int k=0; k<PercentDone-ProgressDots;k++)
		{
			cout<<":";
			cout.flush();
		}
		ProgressDots=PercentDone;
		StepCount++;
	}
	TNow=TNext;
	cout << "]100%\n\n";
	//IntegratorStats->PrintStats();


	//=====================================================
	//Need to figure out why this  doesn't work as intended
	//=====================================================
    	//IntegratorStats *integratorStats = integrator->Integrate(
	//StepSize,InitTime, FinalTime, NumBands, y, KrylovTol, startingBasisSizes);
    	//Print the statistics
    	//printf("Run stats:\n");
    	//integratorStats->PrintStats();



	//=======================================
	//Console Output
	//=======================================
	cout<<"===========================Data========================\n";
	cout <<"y=" << data[0] << "\t\t z=" << data[1] <<"\t\t Temp=" << data[2]<<endl;
	if (LastStepAdjusted == 1)
	{
		cout<<"We altered the last time step to not go over the final time\n";
	}
	cout << "Simulation ended at time: " << TNow << endl;
	cout << "Steps taken: " << StepCount <<endl;
	//cout <<"Sum of y=" <<data[0]+ data [1] + data [2]<<"\t\t";
	//cout <<"Error in sum=" << data[0]+data[1]+data[2]-Y0[0]-Y0[1]-Y0[2]<<endl;
	cout <<"===================Printing data to file==============\n";


	//=======================================
	//Data file output
	//=======================================
	ofstream myfile("TChemData.txt", std::ios_base::app);
	myfile << setprecision(17) << fixed << data[0] << "\t\t" << data[1] << "\t\t";
	myfile << data[2] << "\t\t" <<StepSize <<endl;
	myfile.close();
	cout << "======================Simulation complete=============\n";

   	 // Clean up the integrator
    	delete integrator;



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
		return Steps;
	}else{
                cout<<"Cannot perform non-integer number of steps!!\n";
                exit(1);
        }
}

/*
 * ===========================================================================================
 * 
 * Function InitialConditions
 * 
 * Computes initial condition. 
 * 
 * Output:
 * y0   data
 * 
 * ===========================================================================================
 */

N_Vector InitialConditions()
{
        N_Vector y0 = N_VNew_Serial(NEQ);
        realtype *data = NV_DATA_S(y0);
	data[0]=1.0;
	data[1]=.5*EPS;
	data[2]=1.0;
        return y0;
}


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

int RHS(realtype t, N_Vector u, N_Vector udot, void *userData)
{
	realtype *y = NV_DATA_S(u);
	realtype *dy = NV_DATA_S(udot);
	double Omega= .5 * y[0] * y[1] * exp ( (y[2]-1.0 ) /  (EPS * y[2]  ) );
	dy[0] = -Omega;
	dy[1] = Omega-y[1];
	dy[2] = y[1];
//	cout << dy[0] <<"\t\t" << dy[1] << "\t\t" << dy[2] <<endl;
        return 0;
}
/*
//============================================================
//RHS that uses the TCHEM rhs functions
//============================================================
*/
int RHS_TCHEM(TCHEMPB pb, realtype t, N_Vector u, N_Vector udot, void *userData, ordinal_type number_of_equations)
{
	//problem_type problem;
	realtype *y= NV_DATA_S(u);
	realtype *dy=NV_DATA_S(udot);
	/// we use std vector mimicing users' interface
	using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
	using real_type_1d_view_type = Tines::value_type_1d_view<real_type,host_device_type>;
	using real_type_2d_view_type = Tines::value_type_2d_view<real_type,host_device_type>;
	std::vector<real_type> x_std(number_of_equations, real_type(0));
	x_std[0] = 1200; // temperature
	x_std[1] = 0.1;  // mass fraction species 1
	x_std[2] = 0.9;  // mass fraction species 2

         /// output rhs vector and J matrix
	std::vector<real_type> rhs_std(number_of_equations);
	std::vector<real_type> J_std(number_of_equations * number_of_equations);

	real_type * ptr_x = x_std.data();
	real_type * ptr_rhs = rhs_std.data();
	//real_type * ptr_J = J_std.data();

	real_type_1d_view_type x(ptr_x,   number_of_equations);
	real_type_1d_view_type f(dy, number_of_equations);


        ///
        /// Compute right hand side vector
        ///
	auto member =  Tines::HostSerialTeamMember();
	pb.computeFunction(member, x ,f);
	cout<<"We entered the new RHS function!\n";
	printf("RHS std vector \n" );
	for (ordinal_type i=0;i<number_of_equations;++i)
	{
		printf("%e\n", f[i]);
                std::cout<<"\n";
        }
	return 0;
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
 * 
 * ===========================================================================================
 */

int Jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, void *userData, N_Vector tmp)
{

	realtype *JV = NV_DATA_S(Jv);
	realtype *Y = NV_DATA_S(u);
	realtype *V = NV_DATA_S(v);
        double EXP= exp ( (Y[2] - 1.0 ) / (EPS * Y[2]));
        double DEXP=1/(EPS * Y[2]*Y[2]);
	JV[0] = -.5*EXP* (V[0]*Y[1]+V[1]*Y[0]+V[2]*Y[0]*Y[1]*DEXP);
	JV[1] =-JV[0]-V[1];
	JV[2] = V[1];
	//cout << "====================Data==================\n";
	//cout << Y[0] << "\t\t" << Y[1] << "\t\t" << Y[2] <<endl;
	/*
	for (int i=0; i<NEQ; i++){
		if (isnan(JV[i])==1){
			cout <<"-------------BOOOOOMMM!---------\n";
			cout <<"-----------Dumping data---------\n";
			cout <<"Jtv[0]=" << JV[0]<< "\t\t" << "Jtv[1]=" <<JV[1] <<"\t\t Jtv[2]="<< JV[2]<<endl;
			cout <<"y0="<<Y[0] <<"\t\ty1=" <<Y[1] << "\t\ty2=" << Y[2]<<endl;
			exit(1);
			}
		}*/
        return 0;
}
