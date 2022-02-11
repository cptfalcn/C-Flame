/*
 * -------------------------------------------------------------------------------------------
 * 
 * This is serial implementation. 
 * 
 * ===========================================================================================
 */

#include <stdio.h>
#include <math.h>
#include "Epic.h"
#include "Epi3_KIOPS.h"
#include "Epi3SC_KIOPS.h"
#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <chrono>
// Number of points in the grid
#define NEQ 3
#define EPS 1e-2
N_Vector InitialConditions();
int RHS(realtype, N_Vector, N_Vector, void *);
int Jtv(N_Vector, N_Vector, realtype, N_Vector, N_Vector , void *, N_Vector);
int CheckStep(realtype, realtype);

using namespace std;

int main(int argc, char* argv[])
{
	//====================
	// Declarations
	//====================
	static const realtype FinalTime = 1e1;
	static const int NumBands = 3;
	string MyFile = "OutPutFile.txt";
	realtype StepSize=1e-2;
	int ProgressDots=0; //From 0 to 100; only care about percentage
	int StepCount = 0;
	realtype PercentDone=0;
	void *userData=nullptr;
	realtype TNow=0;
	realtype TNext=0;
	N_Vector y = N_VNew_Serial(NEQ); //The data y
	N_VScale(1.0, InitialConditions(), y);
	realtype *data = NV_DATA_S(y);
	string opt = "EPI2" ;
	bool LastStepAdjusted = 0;
	realtype KTime = 0;
	realtype PrintStepSize = 0;
	realtype KrylovTol = 1e-14;
	//Input Parsing
	if(argc == 1)
	{
		//CheckStep(FinalTime, StepSize);
	}
	else if(argc== 2)
	{
		StepSize=stof(argv[1]);
		std :: cout << "Step-size: " << StepSize << endl;
		//CheckStep(FinalTime, StepSize);
		opt = "EPI2";
	}
	else if(argc == 3)
	{
		StepSize=stof(argv[1]);
		std :: cout << "Step-size: " << StepSize << endl;
		opt = argv[2];
		cout << "Method: " << opt << endl;
	}
	else if(argc == 4 )
	{
		StepSize = stof(argv[1]);
		std :: cout << "Step-size: " << StepSize << endl;
		opt = argv[2];
		cout << "Method: " << opt << endl;
		MyFile = argv[3];
		cout << "Output file: " << MyFile << endl;
	}
	else if(argc == 5)
	{
		StepSize = stof(argv[1]);
		opt = argv[2];
		MyFile = argv[3];
		KrylovTol = stof(argv[4]);
	}
	else
	{
		cout << "!!!Failure:  Incorrect number of args!!!" << endl;
		exit(EXIT_FAILURE);
	}
	PrintStepSize=StepSize;
	//===================================================
	// Create the integrator
	// In this case we use EpiRK5C, for other integrators just use their name.
    	//===================================================
	const int MaxKrylovIters = 500;
	//EpiRK4SV *integrator = new EpiRK4SV(RHS,
	//EpiRK4SC_KIOPS *integrator = new EpiRK4SC_KIOPS(RHS,
	//EpiRK5C *integrator = new EpiRK5C(RHS,
	IntegratorStats *Stats = NULL;
    	Epi2_KIOPS *integrator 	= new Epi2_KIOPS(RHS, Jtv, userData, MaxKrylovIters, y, NEQ);
	Epi3_KIOPS * Epi3 	= new Epi3_KIOPS(RHS, Jtv, userData, MaxKrylovIters, y, NEQ);
	Epi3SC_KIOPS * Epi3SC  	= new Epi3SC_KIOPS(RHS, Jtv, userData, MaxKrylovIters, y, NEQ);
	//========================
	// Set integrator parameters
	//========================
	//const realtype KrylovTol = RCONST(1.0e-14);//1e-14
	int startingBasis[] = {3, 3};



	//=================================
    	// Run the time integrators outside of the loop
	//=================================
	/**/
	if(opt == "EPI3SC")
	{
		auto Start=std::chrono::high_resolution_clock::now();//Time integrator
		Stats = Epi3SC->Integrate(StepSize, StepSize*100, 1e-5, 1e-8, 0, FinalTime,
						NumBands, startingBasis, y);
		auto Stop=std::chrono::high_resolution_clock::now();
                auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
                KTime+=Pass.count()/1e9;
	}
	/**/
	TNext = StepSize;
	if(TNext>FinalTime)
		StepSize=FinalTime;

	cout<<"[";
	while(TNow<FinalTime)
	{
		//Integrate
		if(opt=="EPI2")
		{
			auto Start=std::chrono::high_resolution_clock::now();
			Stats= integrator->Integrate(StepSize,TNow, TNext, NumBands, y,
						KrylovTol, startingBasis);
			auto Stop=std::chrono::high_resolution_clock::now();
                        auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
                        KTime+=Pass.count()/1e9;

		}
		else if(opt=="EPI3")
		{
			auto Start=std::chrono::high_resolution_clock::now();
			Stats= Epi3->Integrate(StepSize, TNow, TNext, NumBands, y,
						KrylovTol, startingBasis);
			auto Stop=std::chrono::high_resolution_clock::now();
                        auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
                        KTime+=Pass.count()/1e9;
		}
		else if(opt=="EPI3SC"){
			break;
		}
		else
		{
			std :: cout <<"!!!Error:  Invalid Integrator selected\n";
			exit(EXIT_FAILURE);
		}

		//Removing this cleaning loop will cause stalls.
		for (int j=0; j<NEQ; j++)
			if(data[j]<0)
				data[j]=0;
		//Track the progress
		PercentDone=floor(TNext/FinalTime*100);
		for (int k=0; k<PercentDone-ProgressDots;k++)
		{
			cout<<":";
			cout.flush();
		}
		ProgressDots=PercentDone;

		//Set next iterate information
		TNow 	+= StepSize;
		TNext 	=  TNow+StepSize;
		if(TNext>FinalTime)
		{
			StepSize=FinalTime-TNow;
			LastStepAdjusted=1;
			TNext = TNow+StepSize;
		}
		StepCount ++;
		//End loop
	}
	cout << "]100%\n\n";
	/**/

	printf("Run stats:\n");
    	Stats->PrintStats();



	//=======================================
	//Console Output
	//=======================================
	cout<<"===========================Data========================\n";
	cout << setprecision(17);
	cout <<"y=" << data[0] << "\t\t z=" << data[1] <<"\t\t Temp=" << data[2]<<endl;
	if (LastStepAdjusted == 1)
		cout<<"Altered final ";
	cout << "Simulation ended at time: " << setprecision(16) << TNow << endl;
	cout << "Steps taken: " << StepCount <<endl;
	cout << "Integration time: " <<KTime << endl;
	cout << "StepSize: " << StepSize << endl;
	//cout <<"Sum of y=" <<data[0]+ data [1] + data [2]<<"\t\t";
	//cout <<"Error in sum=" << data[0]+data[1]+data[2]-Y0[0]-Y0[1]-Y0[2]<<endl;
	cout <<"===================Printing data to " <<MyFile <<"==============\n";


	//=======================================
	//Data file output
	//=======================================
	ofstream myfile(MyFile, std::ios_base::app);
	myfile << setprecision(17) << fixed << data[0] << "\t\t" << data[1] << "\t\t";
	myfile << data[2] << "\t\t" << PrintStepSize << "\t\t" << KTime;
	myfile << "\t\t" << KrylovTol;
	myfile << "\t\t" << Stats->numTimeSteps << endl;
	myfile.close();
	cout << "======================Simulation complete=============\n";

   	 // Clean up the integrator
    	delete integrator;
	delete Epi3SC;



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
        }else if (abs(round(FinalTime/StepSize))-FinalTime/StepSize <1e-3 ){
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
	data[0]=1.0-.5*EPS;
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
