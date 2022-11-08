#include "Epi3VCntrl.h"
#include <cmath>

//Set the necessary userdata class used in zeroD.  THIS MUST MATCH MAIN.
using value_type = realtype;
#define TCHEMPB TChem::Impl::IgnitionZeroD_Problem<value_type, Tines::UseThisDevice<TChem::host_exec_space>::type > 
using namespace std;
class myPb : public TCHEMPB{
	public:
	//members
	ordinal_type  			num_equations;
	N_Vector 				Jac;
	realtype 				t;
	realtype 				MaxStepTaken;
	realtype 				MinStepTaken;
	realtype 				ignTime;
	realtype 				KiopsTime;
	SUNMatrix				Mat;
	int 					Movie;
	std :: string			dumpJacFile;
	int						InternalSteps;
	int						BadErrSteps;
	int						BlowupSteps;
	int						KiopsBlowups;
	std :: string			stepRatioFile;
	realtype 				stepRatio;
	realtype				ProjectTime;
	realtype				OrthogTime;
	//functions
	ordinal_type			get_num_equations(void)
	{
		return this->num_equations;
	}
};

// //Minifactorial function
// int factorial3(int n) { return (n == 1 || n == 0) ? 1 : factorial3(n - 1) * n; }



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
Epi3VCntrl::Epi3VCntrl(CVRhsFn f, CVSpilsJacTimesVecFn jtv, CVLsJacFn jacf, void *userData, int maxKrylovIters,
					N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    	this->f = f;
    	this->jtv = jtv;
		this->jacf = jacf;
    	this->userData = userData;

    	krylov = new Kiops(MaxPhiOrder1 + 2, maxKrylovIters, tmpl, NEQ);//Originally MaxPhiOrder1+1
		NewKrylov= new KiopsRetCode(MaxPhiOrder1 + 2, maxKrylovIters, tmpl, NEQ);//Originally MaxPhiOrder1+1

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
//==========
//Destructor
//==========
Epi3VCntrl::~Epi3VCntrl()
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



//Cleaning function
void Epi3VCntrl::Clean(N_Vector y, int length)
{
	realtype * yD = N_VGetArrayPointer(y);
	for(int i = 0;   i<length;  i++)
		if(yD[i]<0)
		{
			yD[i]=0;
		}

}


//is nan blowup check

void Epi3VCntrl::CheckNanBlowUp(N_Vector State, int vecLength)
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



int Epi3VCntrl::TrackProgress(realtype FinalTime, realtype TNext, realtype PercentDone, int ProgressDots)
{
        PercentDone=floor(TNext/FinalTime*100);
        for (int k=0; k<PercentDone-ProgressDots;k++)
        {
                cout<<":";
                cout.flush();
        }
        return PercentDone;
}




//Playing with a new cleaned up integration with a modified kiops
IntegratorStats *Epi3VCntrl::NewIntegrate(const realtype hStart, const realtype hMax, const realtype absTol,
			const realtype relTol, const realtype t0, const realtype tFinal,
			int basisSizes[], N_Vector y)
{
	realtype * 	data	= N_VGetArrayPointer(y);		//Query the state in debugger with this
	realtype fac		= pow(0.25, .5);
	realtype PercentDone= 0;
	realtype hOld		= hStart;
	realtype hOldOld	= hStart;
	realtype ErrOld		= 1;
	realtype ErrOldOld	= 1;
	realtype ErrNorm	= 0;
	realtype ConfR		= .9;
	int PercentDots		= 0;
	int ProgressDots	= 0;
	myPb * pb{static_cast<myPb *> (userData)};    		//Recast
	ofstream myfile;
	pb->MaxStepTaken		= 0;
	pb->MinStepTaken		= 1;

	ofstream datafile;

	//Debugging output files
	//Need RHS(y), hJac, y, and Remainder.
	ofstream	StepFile;
	StepFile.open(pb->stepRatioFile, std::ios_base::app);
	int IgnDelay		= 0;
	int KiopsErrors		= 0;		//When Kiops throws an error
	int NaNPhiCount		= 0;		//When error is getting bigger than 1e4
	int ErrEstPoor		= 0;		//When poor error but not exploding
	int IntSteps		= 0;		//Integration steps
	int InternalSteps	= 0;		//Total internal steps (good and bad)
	realtype PhiNorm	= 0;
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
	int retVal		= 0;
	int ForceRej 	= 0;
	realtype TOL 	= 5e-7; //9e-2;
	//TOL 			= (absTol+relTol)/2.0;
	//Main integration loop
	while(t<tFinal)
    {
		realtype Err=5.0;
		realtype ErrEst=0;
		IntSteps ++;
		f(t, y, fy, userData); 							// f(t, y) = fy
		this->jacf(t, y, y,     pb->Mat, this->userData, y,y,y);
		JTimesV jtimesv(jtv, f, delta, t, y, fy, userData, tmpVec);
		while(Err > 1.0 )//was 1
		{//Iterate until error is low enough
			
			if(InternalSteps>1)//Output the step analysis to a new file
			{
				StepFile << t <<"\t"<< InternalSteps << "\t" <<hNew;
				StepFile << "\t" <<  h/hNew <<"\t";
			}
			InternalSteps ++;
			hOldOld= hOld;
			hOld= h;
			h=hNew;//h = k;
			N_VScale(h, fy, hfy); 				//Scale f(y);
			
			//Mayya's new method//
			// Y1= y_n + phi_1(6/8 hJ) h f(y_n)
			// R(z)= f(z)-f(y_n) - J*(z-y_n)
			// y(n+1)= y_n + phi_1(hJ) h f(y_n) + 2 phi_3(hj) h r(Y1)
			N_Vector stage1InputVecs[] 	= {zeroVec, hfy}; //Set the b vector
			N_Vector stage1OutputVecs[] = {r1,r2}; //Set output vectors
			N_Vector stage2OutputVecs[]	= {r3};
			realtype timePts[] 			= {6.0/8.0, 1.0};
			realtype timePts2[]			= {1.0};

			retVal = NewKrylov->ComputeKry(2, stage1InputVecs, timePts, 2, stage1OutputVecs, &jtimesv,
				h, krylovTol, basisSizes[0], &integratorStats->krylovStats[0]);
			if(retVal!=0)
			{
				ForceRej 	= 1;
				retVal		= 1;
				KiopsErrors ++;
			}

			//Set  r1= 6/8 phi_1(6/8 hJ), r2=phi_1(hJ)
			N_VScale(8.0/6.0, r1, r1);              //Set r1=phi_1 (6/8 hJ) h f_n //PBFlag
			N_VLinearSum(1.0, y, 1.0, r1, Y1);      //Set Y1= y_n + phi_1(6/8 hJ) h f(y_n)= y_n + r1
			f(t,Y1, fY1, userData);                 //f(t, Y1)= fY1
			jtimesv.ComputeJv(r1,Remainder);        //Set Remainder = J (Y1-y_n)= J (r1)
			N_VLinearSum(-1.0, fy, -1.0, Remainder, Remainder);//Remainder = J(r1) - f(y_n)
			N_VLinearSum(1.0, Remainder, 1.0, fY1, Remainder);//Remainder = R(Y1) = J(r1) - f(y_n) + f(Y1)
			N_VScale(h,Remainder,Remainder);        //set Remainder= h R(Y1)
			N_Vector stage2InputVecs[]= {zeroVec, zeroVec, zeroVec, Remainder}; //[0,0,0,hR(Y1)]
			//Run Kiops again.
			retVal = NewKrylov->ComputeKry(4, stage2InputVecs, timePts2, 1, stage2OutputVecs, &jtimesv,
					h, krylovTol, basisSizes[0], &integratorStats->krylovStats[0]);

			if(retVal!=0)
			{
				//cout << "We have returned an error in Phi3\n";
				ForceRej 	= 1;
				retVal		= 0;
				KiopsErrors ++;
			}
			if(ForceRej ==1)//If we forced a step rejection
			{
				Err			= 1000;
				hNew 		= h/2;
				ForceRej 	= 0;
			}
			else
			{
				//this->CheckNanBlowUp(r3, N_VGetLength(r3));
				PhiNorm = N_VL1Norm(r3);
				N_VScale(2.0, r3, r3);                       	//R3 is also the error estimate
				//get final err est for next step
				//Testing estimator
				ErrNorm 	= sqrt(N_VDotProd(r3,r3));
				//Mayya's estimator

				N_VAbs(y, Scratch1);							//Scratch1 sp to be high order
				N_VScale(relTol, Scratch1, Scratch1);			//relTol*|y|->Scratch1
				N_VAddConst(Scratch1, absTol, Scratch1);		//relTol*|y|+absTol ->Scratch1
				N_VDiv(r3, Scratch1, Scratch1);					//ErrEst/(relTol*|y| + absTol)->Scratch1
				ErrEst = N_VDotProd(Scratch1, Scratch1);		//dot(ErrEst)
				ErrEst = ErrEst/ N_VGetLength(r3);				//Normalize ErrEst
				ErrEst = EPICRSqrt(ErrEst);						//sqrt ErrEst
				Err = ErrEst;                               	//Finalize Error Estimate


				if(PhiNorm >5)
					NaNPhiCount ++;
				else if(PhiNorm<5 && Err>1)
					ErrEstPoor ++;
				//============================
				//Past this point errors arise
				//============================
				//This one helps Gri but none of the others
				//hNew	= max(.5, 1-pow( (absTol+relTol)/2.0, 1.0/4.0 ) )* h * pow(ErrEst, -1.0/ 2.0)  ; 

				
				//Soderlind's control (but with EPS Error Per Step)
				//hNew	= 0.9 * h *pow( TOL/(h*ErrEst), 1.0/4.0);   //Slow
				
				//hNew	= max(.5, 1-pow( (absTol+relTol)/2.0, 1.0/4.0 ) )* h *pow( /(ErrEst), 5.0/12.0) * pow( h/hOld, -1.0/4.0)*pow(hOld/hOldOld, -1.0/8.0);   //Slow

				//hNew	= 0.9 * h * pow(1.0/(ErrEst), 1.0/2.0);    		//Works but slow

				hNew = 0.9 * h * 1.0 * pow(ErrEst, -1.0/ 4.0);                  //Standard method
				ErrOldOld = ErrOld;
				ErrOld = pow(ErrEst, -1.0/ 2.0);
				if(hNew>100*h)//we increase the time
					hNew= 2*h;			//throttle 2x
				if(1000*hNew<h)
					hNew= 0.01*h;		//throttle
				if( hNew>hMax)
					hNew=hMax;
				if(hNew<1e-15)//if( hNew <= ZERO)
				{
					//Perform Data dump
					printf("There is a possible singularity in the solution\n");
					std :: cout << "time stamp: " << t << std :: endl;
					std :: cout << "hNew: " << hNew<< std :: endl;
					std :: cout << "ErrEst: " << Err << std :: endl;
					std :: cout << "Conf: " << ConfR << std :: endl;
					std :: cout << "y: \n";
					//N_VPrint_Serial(y);
					exit(EXIT_FAILURE);
				}
			}
		}//Exit Adaptive loop
		//Step accepted, Do the following update
		N_VLinearSum(1.0, y, 1.0, r2, y); 		// Second Order = y + r2 = y + phi( ) * hfy
		N_VLinearSum(1.0 ,y, 1.0, r3, y); 		// Add the third order component.
		integratorStats->Step();                //Recompute integrator statistics

		this->Clean(y, N_VGetLength(y));		//Clean the data
		t = t + h;                             	//Advance the time 

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
	}//end Loop

	std :: cout << std :: endl;
	pb->InternalSteps	= InternalSteps;
	pb->BadErrSteps		= ErrEstPoor;
	pb->BlowupSteps		= NaNPhiCount;
	pb->KiopsBlowups	= KiopsErrors;

	if(pb->Movie==1)
	{
		myfile.close();
		datafile.close();
	}
	
	StepFile << std :: endl;
	StepFile.close();
	return integratorStats;
}

