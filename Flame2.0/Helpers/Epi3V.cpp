#include "Epi3V.h"
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
	SUNMatrix				Mat;
	int 					Movie;
	std :: string			dumpJacFile;
	int						InternalSteps;
	int						BadErrSteps;
	int						BlowupSteps;
	int						KiopsBlowups;
	//functions
	ordinal_type			get_num_equations(void)
	{
		return this->num_equations;
	}
};

//Minifactorial function
int factorial2(int n) { return (n == 1 || n == 0) ? 1 : factorial2(n - 1) * n; }



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

Epi3VChem_KIOPS::Epi3VChem_KIOPS(CVRhsFn f, CVSpilsJacTimesVecFn jtv, CVLsJacFn jacf, void *userData, int maxKrylovIters,
					N_Vector tmpl, const int vecLength) : NEQ(vecLength)
{
    	this->f = f;
    	this->jtv = jtv;
		this->jacf = jacf;
    	this->userData = userData;

    	krylov = new Kiops(MaxPhiOrder1 + 2, maxKrylovIters, tmpl, NEQ);//Originally MaxPhiOrder1+1
		NewKrylov= new NewKiops(MaxPhiOrder1 + 2, maxKrylovIters, tmpl, NEQ);//Originally MaxPhiOrder1+1

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

//Main integration file
IntegratorStats *Epi3VChem_KIOPS::Integrate(const realtype hStart, const realtype hMax, const realtype absTol,
                const realtype relTol, const realtype t0, const realtype tFinal,
                const int numBands, int basisSizes[], N_Vector y)
{
	realtype * 	data	= N_VGetArrayPointer(y);		//Query the state in debugger with this
	realtype fac		= pow(0.25, .5);
	realtype PercentDone= 0;
	int PercentDots		= 0;
	int ProgressDots	= 0;
	//myPb2 * pb{static_cast<myPb2 *> (userData)};    //Recast
	myPb * pb{static_cast<myPb *> (userData)};    //Recast
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
			{//Iterate until error is low enough
				count ++;
				h=hNew;//h = k;
				// f is RHS function from the problem; y = u_n
				f(t, y, fy, userData); 							// f(t, y) = fy
				N_VScale(h, fy, hfy); 							//Scale f(y);
				this->jacf(t, y, y,     pb->Mat, this->userData, y,y,y);
				JTimesV jtimesv(jtv, f, delta, t, y, fy, userData, tmpVec);

				//Mayya's new method//
				// Y1= y_n + phi_1(6/8 hJ) h f(y_n)
				// R(z)= f(z)-f(y_n) - J*(z-y_n)
				// y(n+1)= y_n + phi_1(hJ) h f(y_n) + 2 phi_3(hj) h r(Y1)
				N_Vector stage1InputVecs[] 	= {zeroVec, hfy}; //Set the b vector
				N_Vector stage1OutputVecs[] = {r1,r2}; //Set output vectors
				N_Vector stage2OutputVecs[]	= {r3};
				realtype timePts[] 			= {6.0/8.0, 1.0};
				realtype timePts2[]			= {1.0};

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
				//get final err est for next step
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




































//Playing with a new cleaned up integration with a modified kiops
IntegratorStats *Epi3VChem_KIOPS::NewIntegrate(const realtype hStart, const realtype hMax, const realtype absTol,
			const realtype relTol, const realtype t0, const realtype tFinal,
			int basisSizes[], N_Vector y)
{
	realtype * 	data	= N_VGetArrayPointer(y);		//Query the state in debugger with this
	realtype fac		= pow(0.25, .5);
	realtype PercentDone= 0;
	int PercentDots		= 0;
	int ProgressDots	= 0;
	myPb * pb{static_cast<myPb *> (userData)};    		//Recast
	//Now need to dump the jac, rhs and y to file on a failure
	//realtype * jacPtr	= N_VGetArrayPointer(pb->Jac);
	//realtype * Rem 		= NV_DATA_S(Remainder);
	ofstream myfile;
	pb->MaxStepTaken		= 0;
	pb->MinStepTaken		= 1;

    //myfile.open(pb->dumpJacFile, std:: ios_base::app);
	ofstream datafile;
	//datafile.open("EPI3VData.txt", std::ios_base::app);

	//Debugging output files
	//Need RHS(y), hJac, y, and Remainder.

	// ofstream	hRHS;
	// ofstream    hJAC;
	// ofstream    Y;
	// ofstream 	REM;

	// hRHS.open("FailedIsoRHS.txt",  std::ios_base::app);
	// hJAC.open("FailedIsohJac.txt",  std::ios_base::app);
	// Y.open("FailedIsoY.txt", std::ios_base::app);
	// REM.open("FailedIsoRemainder.txt",  std::ios_base::app);

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
	//Main integration loop
	while(t<tFinal)
    {
		realtype Err=5.0;
		realtype ErrEst=0;
		IntSteps ++;
		//f(t, y, fy, userData); 							// f(t, y) = fy
		//this->jacf(t, y, y,     pb->Mat, this->userData, y,y,y);
		while(Err > 1 )//was 1
		{//Iterate until error is low enough
			InternalSteps ++;
			h=hNew;//h = k;
			// f is RHS function from the problem; y = u_n
			f(t, y, fy, userData); 							// f(t, y) = fy
			N_VScale(h, fy, hfy); 							//Scale f(y);
			this->jacf(t, y, y,     pb->Mat, this->userData, y,y,y);
			JTimesV jtimesv(jtv, f, delta, t, y, fy, userData, tmpVec);

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
			this->CheckNanBlowUp(r1, N_VGetLength(r1));
			this->CheckNanBlowUp(r2, N_VGetLength(r2));
			//std :: cout << "Second order cleared\n";
			//New comments out next line
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
				//Optional additional printing
				// for (int i = 0; i < N_VGetLength(pb->Jac); i++ )
				// 	hJAC << jacPtr[i] << endl;
				
				// hJAC << endl << endl;

				// realtype * hfyData= NV_DATA_S(hfy);
				// for (int i = 0; i < N_VGetLength(y); i++)
				// {
				// 	hRHS << hfyData[i] << endl;
				// 	Y   << data[i] << endl;
				// 	REM << Rem[i] << endl;
				// }
				// hRHS << endl << endl;
				// Y << endl << endl;

			}
			else
			{
				this->CheckNanBlowUp(r3, N_VGetLength(r3));
				PhiNorm = N_VL1Norm(r3);
				//N_VScale(32.0/9.0,r3, r3);						//New
				N_VScale(2.0, r3, r3);                       	//R3 is also the error estimate
				//get final err est for next step
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
				else if(PhiNorm<5 &&Err>1)
					ErrEstPoor ++;
					//cout << "Rejected Error Step\n";
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
				if(hNew<1e-15)//if( hNew <= ZERO)
				{
					//Perform Data dump
					printf("There is a possible singularity in the solution\n");
					std :: cout << "time stamp: " << t << std :: endl;
					std :: cout << "hNew: " << hNew<< std :: endl;
					std :: cout << "ErrEst: " << Err << std :: endl;
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
		//this->CheckNanBlowUp(y, N_VGetLength(y));		//Check for errors
		t = t + h;                              		//Advance the time 

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
		
		// if(pb->Movie==1)//If we want a movie
		// {
		// 	for(int i = 0 ; i < N_VGetLength(pb->Jac); i ++)
		// 	{
		// 		myfile << jacPtr[i];
		// 		myfile.flush();
		// 		myfile << "\n";
		// 	}
		// 	myfile << h << "\n";
		// 	//Print y data to file
		// 	for(int i = 0 ; i < N_VGetLength(y); i++)
		// 		datafile << data[i] << endl;

		// 	datafile << h << endl;
		// 	datafile << Err << endl;
		// 	datafile << IgnDelay << endl;
		// 	datafile << t << endl;
		// 	datafile << relTol << endl;
		// 	datafile << absTol << endl;
		// }
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

	// std :: cout << "Internal step:" << InternalSteps << std :: endl;
	// std :: cout << "Kiops Errors: " << KiopsErrors << std :: endl;
	// std :: cout << "Integration steps: " << IntSteps << std :: endl;
	// std :: cout << "Phi Blowups: " << NaNPhiCount << std :: endl;
	// std :: cout << "Poor Error Est: " << ErrEstPoor << std :: endl;
	pb->InternalSteps	= InternalSteps;
	pb->BadErrSteps		= ErrEstPoor;
	pb->BlowupSteps		= NaNPhiCount;
	pb->KiopsBlowups	= KiopsErrors;

	if(pb->Movie==1)
	{
		myfile.close();
		datafile.close();
	}
	
	// hRHS.close();
	// hJAC.close();
	// Y.close();
	// REM.close();

	return integratorStats;
}











//The "new" stuff
//Super experimental stuff
//Kiops class
NewKiops::NewKiops(int maxNumVectors, int m_max, N_Vector templateVector, int vecLength) 
    : MaxNumInputVectors(maxNumVectors),
      MaxPhiOrder(maxNumVectors-1),
      VecLength(vecLength),
      M_max(m_max),
      M_min(10),
      Orth_len(vecLength),
      MatrixSize(m_max + 1),
      PhiMatrixSize(MatrixSize + 1)
{
    V = N_VCloneVectorArray(MatrixSize, templateVector);
    if ( V == NULL ) {
        throw std::bad_alloc();
    }

    V_aug = new realtype[MatrixSize*MaxPhiOrder];
    if ( V_aug == NULL ) {
        throw std::bad_alloc();
    }

    scratchVector = N_VClone(templateVector);
    if ( scratchVector == NULL ) {
        throw std::bad_alloc();
    }

    B = N_VCloneVectorArray(MaxPhiOrder, templateVector);
    if ( B == NULL ) {
        throw std::bad_alloc();
    }

    H = new realtype[MatrixSize*MatrixSize];
    if ( H == NULL ) {
        throw std::bad_alloc();
    }
    for (int i = 0; i < MatrixSize*MatrixSize; i++) {
        H[i] = 0.0;
    }

    phiMatrix = new realtype[PhiMatrixSize*PhiMatrixSize];
    if ( phiMatrix == NULL ) {
        throw std::bad_alloc();
    }
    
    phiMatrixSkipped = new realtype[MatrixSize*MatrixSize];
    if ( phiMatrixSkipped == NULL ) {
        throw std::bad_alloc();
    }
    
    w_aug = new realtype[MaxPhiOrder];
    if ( w_aug == NULL ) {
        throw std::bad_alloc();
    }

    expm = new ExpPade(PhiMatrixSize);
}

NewKiops::~NewKiops() {
    N_VDestroyVectorArray(V, M_max);
    N_VDestroyVectorArray(B, MaxPhiOrder);
    N_VDestroy(scratchVector);
    delete[] H;
    delete[] phiMatrix;
    delete[] phiMatrixSkipped;
    delete[] w_aug;
    delete expm;
}


//==================
//Kiops compute phi
//==================
int NewKiops::ComputeKry(const int numVectors, N_Vector* inputVectors, const realtype timePoints[], 
	const int numTimePoints, N_Vector* outputVectors, JTimesV* jtimesv, realtype h, 
	realtype tol, int& m, KrylovStats* krylovStats) 
{
	int p = numVectors - 1;
	//Set optional dump filestreams
	ofstream Hfile;
	ofstream Phifile;
	//Open said filestreams
	// Hfile.open("FailedIsoHMat.txt", std::ios_base::app);
	// Phifile.open("FailedIsoPhiMat.txt", std::ios_base::app);
    // We only allow m to vary between mmin and mmax
    m = max(M_min, min(M_max, m));

    // Initial condition
    N_VScale(1.0, inputVectors[0], outputVectors[0]);
    
    realtype normU = 0.0;
    for (int i = 1; i < numVectors; i++) {
        realtype sum = N_VL1Norm(inputVectors[i]);
        normU = max(normU, sum);
    }

    realtype nu, mu;
    if (normU != 0.0) {
        double ex = ceil(log2(normU));
        nu = pow(2, -ex);
        mu = pow(2, ex);
    } else {
        nu = 1;
        mu = 1;
    }

    for (int i = 0; i < p; i++) {
        N_VScale(nu, inputVectors[p-i], B[i]);
    }

    realtype sign = timePoints[numTimePoints-1] > 0 ? 1 : -1;
    realtype t_now = 0.0;
    realtype t_out = abs(timePoints[numTimePoints-1]);
    realtype tau = t_out;
    bool happy = false;
    int totalMatrixExponentials = 0;

    // Setting the safety factors and tolerance requirements
    realtype gamma, gamma_mmax;
    if (t_out > 1.0) {
        gamma = 0.2;
        gamma_mmax = 0.1;
    } else {
        gamma = 0.9;
        gamma_mmax = 0.6;
    }

    realtype delta = 1.4;

    // Used in the adaptive selection
    int oldm = -1;
    double oldtau = -1;
    double omega = -1;
    bool orderold = true;
    bool kestold = true;

    int l = 0;
    int j = 0;
    double beta;
    int ireject = 0;
    int reject = 0;

    while (t_now < t_out) {
        if (j == 0) {
            for (int k = 1; k < p; k++) {
                int i = p - k;
                w_aug[k - 1] = pow(t_now, i) / factorial2(i) * mu;
            }
            w_aug[p - 1] = mu;

            double sum = N_VDotProd(outputVectors[l],outputVectors[l]);
            for (int i = 0; i < p; i++) {
                sum += w_aug[i] * w_aug[i];
            }
            beta = sqrt(sum);

            N_VScale(1.0/beta, outputVectors[l], V[0]);
            for (int i = 0; i < p; i++) {
                V_aug[i] = w_aug[i] / beta;
            }
        }

        // Incomplete orthogonalization process
        while (j < m) 
		{
            j++;

            // Augmented matrix - vector product
            jtimesv->ComputeJv(V[j-1], V[j]);
            N_VLinearCombination(p, V_aug + (j-1)*p, B, scratchVector);
            N_VLinearSum(h, V[j], 1.0, scratchVector, V[j]);

            for (int k = 0; k < p - 1; k++) {
                V_aug[j*p + k] = V_aug[(j-1)*p + k + 1];
            }
            V_aug[j*p + p - 1] = 0.0;

            // Modified Gram-Schmidt
            for (int i = max(0, j - Orth_len); i < j; i++) {
                int idx = i + (j - 1) * MatrixSize;
                H[idx] = N_VDotProd(V[i], V[j]);
                for (int k = 0; k < p; k++) {
                    H[idx] += V_aug[i*p + k] * V_aug[j*p + k];
                }

                N_VLinearSum(1.0, V[j], -H[idx], V[i], V[j]);
                for (int k = 0; k < p; k++) {
                    V_aug[j*p + k] -= H[idx] * V_aug[i*p + k];
                }
            }

            double nrm = N_VDotProd(V[j],V[j]);
            for (int k = 0; k < p; k++) {
                nrm += V_aug[j*p + k] * V_aug[j*p + k];
            }
            nrm = sqrt(nrm);
			if(isnan(nrm))
			{
				//cout << "norm is NaN\n";
			}
            // Happy breakdown
            if (nrm < tol) {
                happy = true;
                break;
            }

            H[j + (j - 1) * MatrixSize] = nrm;
            N_VScale(1.0 / nrm, V[j], V[j]);
            for (int k = 0; k < p; k++) {
                V_aug[j*p + k] /= nrm;
            }
        }


        // setup PhiMatrix
        for (int i = 0; i < j; i++) {
           for (int k = 0; k < j; k++) {
                phiMatrix[k + i * (j + 1)] = tau * sign * H[k + i * MatrixSize];
            }
            phiMatrix[j + i * (j + 1)] = 0.0;
        }
        
        phiMatrix[j * (j + 1)] = tau * sign;
        for (int k = 1; k < j + 1; k++) {
            phiMatrix[k + j * (j + 1)] = 0.0;
        }
        
        
        // Compute the exponential of the augmented matrix
        expm->Compute(j + 1, phiMatrix);
        totalMatrixExponentials++;
        
        double m_new, tau_new;
        if (happy) {
            // Happy breakdown; wrap up
            omega = 0;
            happy = false;
            m_new = m;
            tau_new = min(t_out - (t_now + tau), tau);
        } else 
		{
            // Local truncation error estimation
            double err = abs(beta * H[j + (j - 1) * MatrixSize] * phiMatrix[(j - 1) + j * (j + 1)]);
			//======================
			//Additional error state check
			//======================
			if(isnan(err))//New, put a better condition in here
			{
				// cout << "Dumping the phi matrix to PhiMat.txt\n";
				// cout << "Dumping the H matrix to HMat.txt\n"; 
				// std :: cout << "Error State has occurred: Err NaN\n";
				// std :: cout << "Beta: " << beta << "\n";
				// std :: cout << "Phi Mat value: " <<phiMatrix[j-1 + j * (j+1)] << "\n";
				// std :: cout << "H Mat value: " << H[j + (j-1) * MatrixSize] << "\n";
				// cout << expm->Compute(j + 1, phiMatrix) << endl;
				// cout << t_now << endl;
				// for( int k = 0 ; k < MatrixSize*MatrixSize; k++)
				// {
				// 	Phifile << phiMatrix[k] << endl;
				// 	Hfile   << tau*H[k] << endl;
				// }
				// Hfile   << endl << endl;
				// Phifile << endl << endl;
				// cout << "Error: LTE NaN" << endl;
				//exit(EXIT_FAILURE);
				return 1;
			}
			//===============
			//End error state check
			//===============
            // Error for this step
			double oldomega = omega;
			omega = t_out * err / (tau * tol);
            // Estimate order
			double order;
			if (m == oldm && tau != oldtau && ireject >= 1) {
				order = max(1.0, log(omega / oldomega) / log(tau / oldtau));
				orderold = false;
			} else if (orderold || ireject == 0) {
                orderold = true;
                order = j / 4;
            } else {
                orderold = true;
            }

            // Estimate k
            double kest;
            if (m != oldm && tau == oldtau && ireject >= 1) {
                kest = max(1.1, pow(omega / oldomega, 1 / (oldm - m)));
                kestold = false;
            } else if (kestold || ireject == 0) {
                kestold = true;
                kest = 2;
            } else {
                kestold = true;
            }

            double remaining_time = omega > delta ? t_out - t_now : t_out - (t_now + tau);
            // Krylov adaptivity
            double same_tau = min(remaining_time, tau);
            double tau_opt = tau * pow(gamma / omega, 1.0 / order);

			if(isnan(tau_opt)||tau_opt<1e-10)		//Added by me to try to stop endless cycling conditions.
			{
				tau_opt= tau/100;	//Same as above.
			}
            tau_opt = min(remaining_time, max(tau / 5.0, min(5.0 * tau, tau_opt)));
	    	int m_opt = ceil(j + log(omega / gamma) / log(kest));
            m_opt = max(M_min, min(M_max, max((int)floor(3.0 / 4.0 * m),
                            min(m_opt, (int)ceil(4.0 / 3.0 * m)))));

            if (j == M_max) {
                if (omega > delta) {
                    m_new = j;
                    tau_new = tau * pow(gamma_mmax / omega, 1 / order);
                    tau_new = min(t_out - t_now, max(tau / 5, tau_new));
                } else {
                    tau_new = tau_opt;
                    m_new = m;
                }
            } else {
                m_new = m_opt;
                tau_new = same_tau;
            }


        }

        // Check error against target
        if (omega <= delta) {
            // Yep, got the required tolerance; update
            krylovStats->NewProjection(j, ireject);
            reject += ireject;

            // Udate for t in the interval (t_now, t_now + tau)
            int blownTs = 0;
            double nextT = t_now + tau;
            for (int k = l; k < numTimePoints; k++) {
                if (abs(timePoints[k]) < abs(nextT)) {
                    blownTs++;
                }
            }

            if (blownTs != 0) {
                // Copy current w to w we continue with.
                N_VScale(1.0, outputVectors[l], outputVectors[l+ blownTs]);

                for (int k = 0; k < blownTs; k++) {
                    double tauPhantom = timePoints[l + k] - t_now;

                    // setup PhiMatrixSkipped
                    for (int i = 0; i < j; i++) {
                        for (int k = 0; k < j; k++) {
                            phiMatrixSkipped[k + i*j] = tauPhantom * sign * H[k + i * MatrixSize];
                        }
                    }
                    
                    // Compute the exponential of the augmented matrix
                    expm->Compute(j, phiMatrixSkipped);
                    totalMatrixExponentials++;

                    N_VLinearCombination(j, phiMatrixSkipped, V, outputVectors[l+k]);
                    N_VScale(beta, outputVectors[l+k], outputVectors[l+k]);
                }

                // Advance l.
                l = l + blownTs;
            }

            // Using the standard scheme
            N_VLinearCombination(j, phiMatrix, V, outputVectors[l]);
            N_VScale(beta, outputVectors[l], outputVectors[l]);
            
            // Update t
            t_now += tau;

            j = 0;
            ireject = 0;

        } else {
            // Nope, try again
            ireject = ireject + 1;
        }

        oldtau = tau;
        tau = tau_new;
		oldm = m;
		m = m_new;
		if (ireject > 2*VecLength)
		{
			//std::cout<<"==========KIOPS has stalled, rejects========="<<std::endl;
			// std::cout<<"==============Data=================="<<std::endl;
			// std::cout<<"t_now="<<t_now<<"\t\t"<< "tau="<<tau<<endl;
            //             std::cout<<"Pade norm="<<expm->Compute(j+1,phiMatrix)<<endl;
			// std::cout<<"rejects "<< ireject<<endl;
			// std::cout << "Matrix Size: " << MatrixSize << endl;
			// std::cout << "Omega: " << omega << endl;
			return 1;
		}
		if (tau_new==0&&t_now!=t_out)
		{
			
			//std::cout << "\n============KIOPS stalled, tau==========" << std::endl;
			// std::cout << "t_now=" << t_now << "\t\t" << "tau=" << tau << endl;
			// std::cout << "Pade norm=" << expm->Compute(j+1,phiMatrix) << endl;
			// std::cout << "Omega: " << omega << endl;
			return 1;
		}
    }
    krylovStats->numMatrixExponentials = totalMatrixExponentials;

	return 0;
}




//Playing with a new cleaned up integration with a modified kiops and full jtv
IntegratorStats *Epi3VChem_KIOPS::NewIntegrateNoTChem(const realtype hStart, const realtype hMax, const realtype absTol,
			const realtype relTol, const realtype t0, const realtype tFinal,
			int basisSizes[], N_Vector y)
{
	realtype * 	data	= N_VGetArrayPointer(y);		//Query the state in debugger with this
	realtype fac		= pow(0.25, .5);
	realtype PercentDone= 0;
	int PercentDots		= 0;
	int ProgressDots	= 0;
	//int IgnDelay		= 0;
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
	//Main integration loop
	auto Start=std::chrono::high_resolution_clock::now();
	while(t<tFinal)
	{
		realtype Err=5.0;
		realtype ErrEst=0;
		while(Err > 1 )//was 1
		{//Iterate until error is low enough
			h=hNew;//h = k;
			// f is RHS function from the problem; y = u_n
			f(t, y, fy, userData); 							// f(t, y) = fy
			N_VScale(h, fy, hfy); 							//Scale f(y);
			//this->jacf(t, y, y,     pb->Mat, this->userData, y,y,y);
			JTimesV jtimesv(jtv, f, delta, t, y, fy, userData, tmpVec);

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
				ForceRej = 1;

			//Set  r1= 6/8 phi_1(6/8 hJ), r2=phi_1(hJ)
			//Trying Val's new method
			//N_VScale(8.0/6.0, r1, r1);              //Old Set r1=phi_1 (6/8 hJ) h f_n //PBFlag
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
				ForceRej = 1;

			if(ForceRej ==1)//If we forced a step rejection
			{
				Err			= 1000;
				hNew 		= h/2;
				ForceRej 	= 0;
			}
			else
			{
				N_VScale(32.0/9.0, r3, r3);
				//N_VScale(2.0, r3, r3);                       	//Old R3 is also the error estimate
				//get final err est for next step
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
				if(hNew<1e-15)//if( hNew <= ZERO)
				{
					//Perform Data dump
					printf("There is a possible singularity in the solution\n");
					std :: cout << "time stamp: " << t << std :: endl;
					std :: cout << "hNew: " << hNew<< std :: endl;
					std :: cout << "ErrEst: " << Err << std :: endl;
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
		t = t + h;                              //Advance the time 
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
	auto Stop=std::chrono::high_resolution_clock::now();
	auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
	cout << endl;
	cout << Pass.count()/1e9 << endl;
	return integratorStats;
}