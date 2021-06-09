/*
 * ===========================================================================================
 * 
 * In this example we show how Epic package can be used with test problems which are not 
 * written using objects. 
 * 
 * -------------------------------------------------------------------------------------------
 * 
 * This is implementaion of Burgers equation for Epic,
 * 
 *                                    u_t + u u_x = Nu u_xx
 * 
 * with Dirichlet boundary conditions and initial condition is given by
 * 
 *                               u = (sin(3 pi x))^3 (1-x)^(3/2)
 * 
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

// Number of points in the grid
#define NEQ 3
#define EPS 1e-2
N_Vector InitialConditions();
void Lu(N_Vector, N_Vector);
void Mu(realtype, N_Vector, N_Vector);
int RHS(realtype, N_Vector, N_Vector, void *);
int Jtv(N_Vector, N_Vector, realtype, N_Vector, N_Vector , void *, N_Vector);
int CheckStep(realtype, realtype);

using namespace std;

int main()
{
    // Declaration of start and finish time
    //static const realtype InitTime = 0.0;
    static const realtype FinalTime = 10.0;
    static const int NumBands = 3;
    static const realtype StepSize=1e-1;
	//check if we can make the stepsize
	static const int Steps=CheckStep(FinalTime, StepSize);
	void *userData=nullptr;
	realtype TNow=0;
	realtype TNext=0;
    // Setup the initial condition (and solution) y
    N_Vector y = N_VNew_Serial(NEQ);
	N_VScale(1., InitialConditions(), y);
	realtype *data = NV_DATA_S(y);
	//=====================================
	//This block is testing, remove later
	//actually, we need some of this.
	//=====================================
	//N_Vector Is = N_VNew_Serial(NEQ);
	//N_Vector Jv = N_VNew_Serial(NEQ);
	//N_Vector rhs= N_VNew_Serial(NEQ);
	N_Vector y0 = N_VNew_Serial(NEQ);
	//realtype *IS = NV_DATA_S(Is);
	//realtype *JV = NV_DATA_S(Jv);
	//realtype *RHSy=NV_DATA_S(rhs);
	realtype *Y0= NV_DATA_S(y0);
	N_VScale(1., InitialConditions(), y0);
	//for (int i=0; i<NEQ; i++){
	//	IS[i]=1.0;}
	//RHS(0, y, rhs, userData);
	//Jtv(Is, Jv, 1, y , rhs, userData, Is);
	//cout << "=================Initial Data===================\n";
	cout << "=========================y======================\n";
	cout << "y[0]=" << data[0] <<"\t\t y[1]=" << data[1] << "\t\t y[2]=" << data[2] << endl;
	//cout << "======================rhs========================\n";
	//cout << "rhs[0]=" << RHSy[0] << "\t\t rhs[1]=" << RHSy[1] << "\t\t rhs[2]="<<RHSy[2] <<endl;
	//cout << "========================JtV=====================\n";
	//cout <<  "Jtv[0] =" <<JV[0] <<"\t\t Jtv[1]="<< JV[1] <<"\t\t Jtv[2]=" << JV[2]  <<endl;
	//end test 

    // Create the integrator
    // In this case we use EpiRK5C, for other integrators just use their name like we did 
    const int MaxKrylovIters = 500;
    Epi2_KIOPS *integrator = new Epi2_KIOPS(RHS,
                                      Jtv,
                                      userData,
                                      MaxKrylovIters,
                                      y,
                                      NEQ);

    // Set run parameters
    //const realtype StepSize = RCONST(4e-5);
    const realtype KrylovTol = RCONST(1.0e-14);//1e-14
    int startingBasisSizes[] = {3, 3};
    
    // Run the integrator loop
	for (int i=0; i<Steps; i++){
		TNow= i*StepSize;
		TNext=(i+1)*StepSize;
		//cout<< "Completed step " << i+1 << endl;
		//cout<< "Stepped from " <<TNow << " to " <<TNext <<endl;
		integrator->Integrate(
        	StepSize,TNow, TNext, NumBands, y,
		KrylovTol, startingBasisSizes);
		//If the data goes negative, which is not physical, correct it.
		//Removing this loop will cause stalls.
		for (int j=0; j<NEQ; j++){
			if(data[j]<0){
				data[j]=0;
			}
		}

	}



    	//IntegratorStats *integratorStats = integrator->Integrate(
	//StepSize,InitTime, FinalTime, NumBands, y, KrylovTol, startingBasisSizes);
    	// Print the statistics
    	//printf("Run stats:\n");
    	//integratorStats->PrintStats();




    printf("\n\n");
	cout<<"===========================Data========================\n";
	cout <<"y=" << data[0] << "\t\t z=" << data[1] <<"\t\t Temp=" << data[2]<<endl;
	cout <<"Sum of y=" <<data[0]+ data [1] + data [2]<<endl;
	cout <<"Error in sum=" << data[0]+data[1]+data[2]-Y0[0]-Y0[1]-Y0[2]<<endl;
    // Clean up the integrator
    delete integrator;

}



int CheckStep(realtype FinalTime, realtype StepSize){
        if(fmod(FinalTime,1/StepSize) ==0){
                static const int Steps= FinalTime/StepSize;
                cout << "We will take " <<Steps << " steps.\n";
		return Steps;
        }else{
                cout<<"Cannot perform non-integer number of steps!!\n";
                cout<<"fmod(FinalTime, StepSize)="<<fmod(FinalTime,1/StepSize)<<endl;
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
