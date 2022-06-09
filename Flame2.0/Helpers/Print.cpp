#include "Print.h"

//====== ____===========================//
//      |  _ \     _                    //
//      | (_) )   (_)         _         //
//      |  __/`._  _  _ __  _| |_       //
//      | |  |  _)| || `  \(_   _)      //
//      | |  | |  | || |\ |  | |        //
//      |_|  |_|  |_||_||_|  |_|        //
//========================================

void PrintExpData(realtype * data, int Experiment, realtype Norm)
{
        if(Experiment==0)
                std :: cout << "y=" << data[0] << "\t\t z=" << data[1] <<"\t\t Temp=" <<data[2] << std :: endl;
        else if(Experiment==1)
        {
                //Temp H2 O2 O OH H2O H HO2 H2O2
                std :: cout << "Temp=" << data[0] << "\t\t H2=" << data[1] <<"\t\t O2=" << data[2]<<std :: endl;
                std :: cout << "Total Mass Fractions: " <<Norm-data[0];
                std :: cout << "\t\t Mass Fraction error: "<<abs( Norm -data[0]-1.0)<<std :: endl;
        }else if (Experiment==2){//Gri
                std :: cout << "Temp=" << data[0] << "\t\t CH4=" << data[14] <<"\t\t O2=" << data[4];
                std :: cout << std :: endl << "Total Mass Fractions: " <<Norm -data[0];
                std :: cout << "\t\t Mass Fraction error: "<<abs( Norm -data[0]-1.0)<< std :: endl;
        }

}

void PrintFromPtr(realtype * ptr, int num_eqs)
{
        cout<<endl;
        for(int i=0; i<num_eqs; i++)
        {
                cout<<ptr[i]<<endl;
        }
        cout<<endl;
}
//====================================
//Print the state to a file
//=====================================
void PrintDataToFile(ofstream & myfile, realtype * data, int number_of_equations, realtype t, string BAR,
			string MyFile, realtype KTime)
{
	//cout << BAR <<"Printing data to "<< MyFile << BAR << endl;
        for (int i=0; i<number_of_equations; i++)
        {
                myfile<<setprecision(20)<<fixed <<data[i] <<"\t\t";
        }
        myfile << "\t\t " << t << "\t\t" << KTime <<endl;//Print the integrator time
        //myfile.flush();

}

//========================
//Print the Parameters
//========================
void PrintExpParam(realtype FinalTime, realtype TNow, realtype StepSize, realtype StepCount,
			realtype KrylovTol, realtype AbsTol, realtype RelTol, realtype KTime, string BAR)
{
	//cout << BAR << "\tSim Parameters\t\t" << BAR << endl;
	cout << "Exact Final Time: " << FinalTime << endl;
	cout << "Simulation Final Time: " << TNow;
	//cout << "Step Size: " << StepSize << "\t\tNumber of Steps: "<<StepCount;
	//cout << "\t\tKrylov Tolerance: " <<KrylovTol<<  endl;
	//cout << "Absolute tolerance: " << AbsTol << "\tRelative tolerance: " << RelTol;
	cout << "\t\tIntegration time: " << KTime <<endl;
}

void PrintPreRun(realtype StepSize, realtype Delx, int Steps, realtype KryTol, realtype AbsTol,
			realtype RelTol, string Method, int num_pts, string Bar)
{
        cout << "Delt: " << StepSize << "\t\tNumber of Steps: "<<Steps;
	cout << "\t\tFinal Time: " << StepSize*Steps << endl << "Delx: " << Delx << "\t\t";
	cout << "Grid Points: " << num_pts;
	if(Delx!=0)
		cout << "\t\t\tX = [0, " << num_pts*Delx << "]";
	cout << endl;
	if(Method == "CVODEKry")
		cout << "Absolute tolerance: " << AbsTol << "\tRelative tolerance: " << RelTol<<endl;
	else
        	cout << "Krylov Tolerance: " <<KryTol << endl;
	PrintMethod(Method, Bar);
}

//================================
//Print Super Data
//================================
void PrintSuperVector(realtype * data, int Experiment, int num_pts, string BAR)
{
	//cout << BAR << "Printing SuperVector" << BAR << endl;
	for(int i = 0; i < num_pts; i++)
	{
		if(Experiment==0)
		{
                	std :: cout << "y=" << data[0+i] << "\t\t z=" << data[num_pts+i];
			std :: cout <<"\t\t Temp=" <<data[2*num_pts + i];
		}
		else if(Experiment==1)
        	{
				//Temp H2 O2 O OH H2O H HO2 H2O2
				std :: cout << "Temp=" << data[i] << "\t\t H2=" << data[num_pts+i];
				std :: cout <<"\t\t O2=" << data[2*num_pts + i];
        	}else if (Experiment==2){//Gri
				std :: cout << "Temp=" << data[i] << "\t\t CH4=" << data[14*num_pts + i];
				std :: cout <<"\t\t O2=" << data[4*num_pts + i];
        	}else if (Experiment==3){//Iso-Oct
				std :: cout << "Temp=" << data[i] << "\t\t Iso-Oct=" << data[742*num_pts + i];
				std :: cout <<"\t\t O2=" << data[4*num_pts + i] << "\t\t CH4=" << data[14*num_pts + i];
			}
		std :: cout << std :: endl;
	}
}

// void PrintProfiling(IntegratorStats* integratorStats, int Profiling, string Method, string Bar, void* cvode_mem)
// {
// 	if (Profiling ==1 && Method!="CVODEKry")
// 	{
// 		cout << Bar << "\tPerformance data\t" << Bar << endl;
//          	integratorStats->PrintStats();
// 	}else if(Profiling==1 && Method=="CVODEKry")
// 	{
// 		long int steps 		= 0;
// 		long int fEval 		= 0;
// 		long int linSetups 	= 0;
// 		long int errTestFails	= 0;
// 		int lastOrd		= 0;
// 		int nextOrd		= 0;
// 		realtype realHInit	= 0;
// 		realtype hLast		= 0;
// 		realtype hCurr		= 0;
// 		realtype tCurr		= 0;
// 		long int nonLinIters	= 0;
// 		long int Projections	= 0;
// 		CVodeGetNumSteps(cvode_mem, &steps);
// 		CVodeGetIntegratorStats(cvode_mem, &steps, &fEval, &linSetups, &errTestFails, &lastOrd,
// 					&nextOrd, &realHInit, &hLast, &hCurr, &tCurr);
// 		CVodeGetNumNonlinSolvIters(cvode_mem, &nonLinIters);
// 		///CVodeGetNumProjEvals(cvode_mem, &Projections);
// 		cout << Bar << "\tPerformance data \t" << Bar << endl;
// 		//cout << "Sorry, no data available right now" << endl;
// 		cout << "CVODE steps:" << steps << "\t\t\n";
// 		cout << "Number of linear up steps: " << linSetups << "\n";
// 		cout << "Last order: " << lastOrd << endl;
// 		cout << "Next order: " << nextOrd << endl;
// 		cout << "Real initial h: " << realHInit << endl;
// 		cout << "Last h: " << hLast << endl;
// 		cout << "Current h: " << hCurr << endl;
// 		cout << "Current t: " << tCurr << endl;
// 		cout << "NonLineraIters: "<< nonLinIters << endl;
// 		//cout << "Projections: " << Projections << endl;
// 	}
// }

void PrintMethod(string Method, string Bar)
{
	if(Method == "EPI2")
		cout << Bar << "\t EPI2 Selected \t" << Bar << endl;
	if(Method == "EPI3")
		cout << Bar << "\t EPI3 Selected \t" << Bar << endl;
	if(Method == "CVODEKry")
		cout << Bar << "\tCVODE Selected \t" << Bar << endl;
	if(Method == "EPIRK4")
		cout << Bar << "\tEPIRK4 Selected\t" << Bar << endl;
}
