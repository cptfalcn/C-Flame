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
void PrintDataToFile(ofstream & myfile, realtype * data, int number_of_equations, realtype t)
{
        for (int i=0; i<number_of_equations; i++)
        {
                myfile<<setprecision(20)<<fixed <<data[i] <<"\t\t";
        }
        myfile << "\t\t " << t;
        //myfile.flush();

}

//========================
//Print the Parameters
//========================
void PrintExpParam(realtype FinalTime, realtype TNow, realtype StepSize, realtype StepCount,
			realtype KrylovTol, realtype AbsTol, realtype RelTol, realtype KTime)
{
	cout << "Exact Final Time: " << FinalTime << "\t\tSimulation Final Time: " <<TNow<<endl;
	cout << "Step Size: " << StepSize << "\t\tNumber of Steps: "<<StepCount<<endl;
	cout << "Krylov Tolerance: " <<KrylovTol<<  endl;
	cout << "Absolute tolerance: " << AbsTol << "\t\tRelative tolerance: " << RelTol << endl;
	cout << "Time integration time: " << KTime <<endl;
}
