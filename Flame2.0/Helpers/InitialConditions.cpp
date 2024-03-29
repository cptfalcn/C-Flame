//#include "Derivatives.h"
#include "InitialConditions.h"
#define BAR "===================="
using namespace std;
//========================================
//   ____  _  _  ____  ____        __   //
//  (_  _)| \| |(_  _)(_  _)      / _)  //
//    ||  |  \ |  ||    ||   _   / /    //
//   _||_ | |  | _||_   ||  (_)  \ \_   //
//  (____)|_|\_|(____)  []        \__)  //
//========================================
void IntConKappa(realtype * DATA, realtype EPS)
{
	DATA[0]=1.0-.5*EPS;
	DATA[1]=.5*EPS;
	DATA[2]=1.0;
	cout << BAR << "\tKapila Experiment\t" << BAR <<endl;
}

void IntConHydro(realtype * data, int SampleNum)
{
	data[0]=900;
 	switch(SampleNum)
        {
        case 1:
                data[1]=2.7431213550e-01;
                data[2]=1-data[1];
                cout<< BAR <<"Hydrogen Experiment Sample 1" << BAR << endl;
                break;
        case 2:
                data[1]=5.926665910131902193e-02;
                data[2]=9.407333408986811030e-01;
                cout << BAR <<"Hydrogen Experiment Sample 2" << BAR << endl;
                break;
        case 3:
                data[1]=1.438365693957185942e-01;
                data[2]=8.561634306042813503e-01;
                cout << BAR <<"Hydrogen Experiment Sample 3" << BAR << endl;
                break;
        case 4: //1.000000000000000000e+03 1.013250000000000000e+05 5.926665910131902193e-02 9.40733340898681103>                //Run to 1e-5 time;
                data[0] = 1000.0;
                data[1] = 5.926665910131902193e-02;
                data[2] = 9.407333408986811030e-01;
                cout << BAR <<"Hydrogen Experiment Sample 2" << BAR << endl;
                break;
        default:
                cout<< BAR << "\tUnknown Sample\t\t" << BAR << endl;
                exit(EXIT_FAILURE);
        }
}

void IntConGri30(realtype * data, int SampleNum)
{
        data[0] = 1000.0;//use 1000 for base
        switch(SampleNum)
        {
        case 1:
                data[0] = 1000.0;//use 1000 for base
                data[14] =  5.483858845025559731e-02;//CH4
                data[4]  =  2.187578062376045740e-01;//O2
                data[48] =  7.137587863547695255e-01;//N2
                data[49] =  1.264481895737025810e-02;//AR
                cout << BAR <<"GRI Experiment Sample 1\t" << BAR <<endl;
                break;
	case 2://Ignition phase
		data [0]  =	1016.00647630225341799814;
		data [1]  =	0.00003165371310189990;
		data [2]  =	0.00000000049134365240;
		data [3]  = 	0.00000000169386903351;
		data [4]  =	0.21644178630274113484;
		data [5]  = 	0.00000000683229870383;
		data [6]  =	0.00152738571981635217;
		data [7]  = 	0.00000575687135275362;
		data [8]  =	0.00001006287733677639;
		data [9]  =	0.0;
		data [10] =	0.00000000000000000005;
		data [11] = 	0.00000000000036293814;
		data [12] = 	0.00000000000005331292;
		data [13] =	0.00000285917677468789;
		data [14] = 	0.05340734010105312179;
		data [15] = 	0.00072507393455497953;
		data [16] =	0.00000820940893000156;
		data [17] = 	0.00000000029929522812;
		data [18] =          0.00095827516681306409;
		data [19] =          0.00000000000491228838;
		data [20] =          0.00000001122933969694;
		data [21] =          0.00002716037638680501;
		data [22] =          0.00000000000000000285;
		data [23] =          0.00000002109699641356;
		data [24] =          0.00000000000145641804;
		data [25] =          0.00004468806554316932;
		data [26] =          0.00000000207838516599;
		data [27] =          0.00040552040884691296;
		data [28] =          0.00000000000005688428;
		data [29] =          0.00000050354078975932;
		data [30] =          0.00000000000086813166;
		data [31] =          0.0;
		data [32] =          0.00000000000000000010;
		data [33] =          0.00000000000000000001;
		data [34] =          0.00000000000000000016;
		data [35] =          0.00000000000000602579;
		data [36] =          0.00000000000000062553;
		data [37] =          0.00000000000000421165;
		data [38] =          0.00000000039778419778;
		data [39] =          0.00000000000000000025;
		data [40] =          0.0;
		data [41] =          0.00000000000000000228;
		data [42] =          0.0;
		data [43] =          0.0;
		data [44] =          0.0;
		data [45] =          0.0;
		data [46] =          0.00000000000000000146;
		data [47] =          0.00000000000000000003;
		data [48] =          0.71375878610153931092;
		data [49] =          0.01264481895737025810;
		data [50] =          0.00000000001282017043;
		data [51] =          0.00000007455627788664;
		data [52] =          0.00000000007115097556;
		data [53] =          0.00000000050972778972;
		cout << BAR << "Gri experiment, Sample 1: Ignition" << BAR <<endl;
		break;
	case 3: //Cooling phase.
		data[0]  = 		2541.73398420370449457550;
		data[1]  =              0.00079713818858765122;
		data[2]	 =          	0.00011615649331851274;
		data[3]  =	        0.00119009194142694042;
		data[4]  =	        0.01524429305687882708;
		data[5]  =          	0.00665025350382286631;
		data[6]  =          	0.11147589800627141776;
		data[7]  = 	        0.00000310465234751613;
		data[8]  =          	0.00000018516455058009;
		data[9] =	        0.00000000000000483520;
		data[10] =		0.00000000000000058037;
		data[11] = 	        0.00000000000000069035;
		data[12] =      	0.00000000000000005258;
		data[13] =           	0.00000000000000133873;
		data[14] =          	0.00000000000000027699;
		data[15] =          	0.02750254224752582410;
		data[16] =          	0.10722384830012078594;
		data[17] =          	0.00000001201934283241;
		data[18] =          	0.00000000012801367522;
		data[19] =          	0.00000000000000133435;
		data[20] =          	0.00000000000000002604;
		data[21] =          	0.00000000000000005443;
		data[22] =          	0.0;
		data[23] =          	0.00000000000000000017;
		data[24] =          	0.0;
		data[25] =       	0.0;
		data[26] =         	0.0;
		data[27] =          	0.0;
		data[28] =          	0.00000000000000001073;
		data[29] =          	0.00000000000000000533;
		data[30] =          	0.00000000000000000001;
		data[31] =          	0.00000018420177294650;
		data[32] =          	0.00000002451731891133;
		data[33] =          	0.00000000576347090342;
		data[34] =           	0.00000000594523165504;
		data[35] =           	0.00000000720817952472;
		data[36] =          	0.00635989207120909918;
		data[37] =          	0.00000206721163002734;
		data[38] =          	0.00000048695089732093;
		data[39] =          	0.00000022656890009984;
		data[40] =          	0.00000000000465967205;
		data[41] =          	0.00000000030656999855;
		data[42] =          	0.00000000000000033243;
		data[43] =          	0.00000000000000000030;
		data[44] =          	0.00000000000000569858;
		data[45] =         	0.00000000001778692080;
		data[46] =          	0.00000000309257649779;
		data[47] =          	0.00000000038330845161;
		data[48] =          	0.71078875309729439014;
		data[49] =          	0.01264481895737025810;
		data[50] =          	0.0;
		data[51] =          	0.0;
		data[52] =          	0.0;
		data[53] =          	0.0;
		cout << BAR << "Gri Sample 1: Post Ign" << BAR << endl;
		break;
	case 4: //Aditya's case, runs at 4 ATM
		data[0] 	=  800;//Temp
        data[14]	=  0.0392661; //CH4
        data[4] 	=  0.2237715;//O2
        data[48]	=  0.736962400000000; //N2
		cout << BAR << "Gri Sample: Hi Pressure" << BAR << endl;
		break;
	case 5: //High temp setup
		data[0]  =  1300;
		data[14] =  5.483858845025559731e-02;//CH4
        data[4]  =  2.187578062376045740e-01;//O2
        data[48] =  7.137587863547695255e-01;//N2
        data[49] =  1.264481895737025810e-02;//AR
        cout << BAR <<"GRI Experiment Sample 2\t" << BAR <<endl;
		break;
	case 6: //Low temp long burn
		cout << BAR << "GRI Experiment: Long burn" << BAR << endl;
		data[0] =   810.0;//use 1000 for base
        data[14] =  5.483858845025559731e-02;//CH4
        data[4]  =  2.187578062376045740e-01;//O2
        data[48] =  7.137587863547695255e-01;//N2
        data[49] =  1.264481895737025810e-02;//AR
        cout << BAR <<"GRI Experiment Sample 1\t" << BAR <<endl;
        break;
	case 60:
		cout << BAR << "Gri Sample: Hi-Argon" << BAR << endl;
		data[0] 	= 1300;
		data[14] 	= 5e-2;//CH4
		data[4]		= 2e-1;//O2
		data[48] 	= 7e-1;//N2
		data[49] 	= 5e-2;//AR
		break;
	case 61:
		cout << BAR << "Gri Sample: Super-Argon" << BAR << endl;
		data[0] 	= 1300;
		data[14] 	= 5e-2;//CH4
		data[4]		= 2e-1;//O2
		data[48] 	= 5e-2;//N2
		data[49] 	= 7e-1;//AR
		break;
	case 62:
		cout << BAR << "Gri Sample: Low-Argon" << BAR << endl;
		data[0] 	= 1400;
		data[14] 	= 3e-1;//CH4
		data[4]		= 5e-1;//O2
		data[48] 	= 1e-1;//N2
		data[49] 	= 1e-1;//AR
		break;		
	case 69:
		cout << BAR << "Gri Sample: 1e-4 Hi-P" << BAR << endl;
		ConditionsFromFile(data, "GriIgnition.txt");
		break;
	case 70:
		cout << BAR << "Gri Sample: Hi-Fuel" << BAR << endl;
		data[0] 	= 1400;
		data[14] 	= 8e-2;//CH4
		data[4]		= 1.2e-1;//O2
		data[48] 	= 8e-1;//N2
		break;
	case 71:
		cout << BAR << "Gri Sample: Super-Fuel" << BAR << endl;
		data[0] 	= 1300;
		data[14] 	= 1e-1;//CH4
		data[4]		= 1e-1;//O2
		data[48] 	= 8e-1;//N2
		break;
	case 80:
		cout << BAR << "Gri Sample: Hydrogen Burn" << BAR << endl;
		data[0] 	= 1300;
		data[2] 	= 1e-1;//H
		data[4]		= 1e-1;//O2
		data[48] 	= 8e-1;//N2
		break;
	case 81:
		cout << BAR << "Gri Sample: Hydrogen Burn low phi" << BAR << endl;
		data[0] 	= 1300;
		data[2] 	= 1e-2;//H
		data[4]		= 1.8e-1;//O2
		data[48] 	= 8e-1;//N2
		break;
	case 10:
		data[0] = 1000.0;//use 1000 for base
		data[14] =  5.483858845025559731e-02;//CH4
		data[4]  =  2.187578062376045740e-01;//O2
		data[48] =  7.137587863547695255e-01;//N2
		data[49] =  1.264481895737025810e-02;//AR
		cout << BAR <<"GRI Experiment Sample 1--Hi Pressure\t" << BAR <<endl;
		break;

    }
}

void IntConIsoOct(realtype * data, int SampleNum)
{
	switch(SampleNum)
	{
	case 1:
		data[0] 	=  800.0;//use 1000 for base
		data[742] 	=  1.654852956820578869e-01;//IC8H18 
		data[4]  	=  1.931486027241477166e-01;//O2
		data[7] 	=  6.302015669181122437e-01;//N2
		data[10] 	=  1.116453467568208867e-02;//AR
		cout << BAR <<"Iso-Oct Experiment Sample 10 atm\t" << BAR << endl;
		break;
	case 2:
		data[0]		= 1250.0;//use 1300 for hot gas
		data[742]	= 1.654852956820578869e-01;//IC8H18
		data[4]		= 1.931486027241477166e-01;//O2
		data[7] 	= 6.302015669181122437e-01;//N2
		data[10] 	= 1.116453467568208867e-02;//AR
		cout << BAR << "Iso-Oct Experiment Sample 10 atm Hi temp\t" << BAR << endl;
		break;
	case 3:
		data[0] 	=  800.0;//use 1000 for base
		data[742] 	=  1.654852956820578869e-01;//IC8H18 
		data[4]  	=  1.931486027241477166e-01;//O2
		data[7] 	=  6.302015669181122437e-01;//N2
		data[10] 	=  1.116453467568208867e-02;//AR
		cout << BAR <<"Iso-Oct Experiment Sample 16 atm\t" << BAR << endl;
		break;
	case 4:
		data[0]		= 1250.0;//use 1300 for hot gas
		data[742]	= 1.654852956820578869e-01;//IC8H18
		data[4]		= 1.931486027241477166e-01;//O2
		data[7] 	= 6.302015669181122437e-01;//N2
		data[10] 	= 1.116453467568208867e-02;//AR
		cout << BAR << "Iso-Oct Experiment Sample 16 atm Hi temp\t" << BAR << endl;
		break;
	}

}

void IntConNDodec(realtype * data, int SampleNum)
{
	switch(SampleNum)
	{
		case 1:
		data[0] 	= 1200;
		data[4]		= 2.169982215727128128e-01;	 //O2
		data[7] 	= 7.080176471629753276e-01;  //N2
		data[67]	= 6.244102126951704690e-02;  //NC12H26
		break;
	}


}

void IntConNButane(realtype * data, int SampleNum)
{
	switch(SampleNum)
	{
		case 1:
		data[0] 	= 1200;
		data[11]	= 2.173895119224689421e-01;	 //O2
		data[155] 	= 7.092943418325631244e-01;  //N2
		data[154] 	= 1.256572768198739760e-02;  //AR
		data[58]	= 6.075041856298054460e-02;  //C4H10
		break;
	}	
}

void IntConNHeptane(realtype * data, int SampleNum)
{//10 atm
	switch(SampleNum)
	{
		case 1:
		data[0] 	= 800;
		data[4]		= 2.171425680562193705e-01;	 //O2
		data[7] 	= 7.084886181086718082e-01;  //N2
		data[10] 	= 1.255145362916586065e-02;  //AR
		data[396]	= 6.181736020594310460e-02;  //C4H10
		break;
	}
}


void SetIntCons(int Experiment, int Sample, realtype* data)
{
	switch(Experiment)
	{
		case 0:
			IntConKappa(data, 1e-2);
			break;
		case 1:
			IntConHydro(data, Sample);
			break;
		case 2:
			IntConGri30(data, Sample);
			break;
		case 3:
			IntConIsoOct(data, Sample);
			break;
		case 4:
			IntConNDodec(data, Sample);
			break;
		case 5:
			IntConNButane(data, Sample);
			break;
		case 6:
			IntConNHeptane(data, Sample);
			break;
    }
}


void SetSuperInitCons(realtype * data, realtype * StateData, int numEqs, int numPts)
{
	for( int i = 0 ; i < numEqs ; i ++ )
	{
		for (int j = 0 ; j < numPts; j ++)
			StateData[j + i * numPts] = data[i];
	}
}


void SinWave(realtype * data, int numPts, int myLength, realtype h)
{//Bugged, solve later
	realtype x = 0;
        int myIndex;
	realtype temp1 = 0;
        for (int i = 0; i < myLength; i++)
        {
                myIndex = i % numPts;
                //cout << myIndex << "\t\t";
                x =  myIndex*h + h/2;
                //cout << x << endl;
                temp1 = cos(2 * M_PI * x)+1;
		//cout << temp1 << endl;
                //realtype temp2 = pow(1-x, 3.0/2.0);
                data[i] = temp1;
        }
}

void SinWave2(realtype * data, int numPts, int myLength, realtype h)
{
        realtype x = 0;
        int myIndex;
        realtype temp1 = 0;
        for (int i = 0; i < myLength; i++)
        {
                myIndex = i % numPts;
                //cout << myIndex << "\t\t";
                x =  myIndex*h + h/2;
                //cout << x << endl;
                temp1 = sin(2 * M_PI * x)+1;
                //cout << temp1 << endl;
                //realtype temp2 = pow(1-x, 3.0/2.0);
                data[i] = temp1;
        }
}

void GuassWave(realtype * data, int numPts, int myLength, realtype h)
{
        realtype x = 0;
        int myIndex;
        realtype temp1 = 0;
	realtype OffSet = h*( round (0.5 * numPts) -0.5);
        for (int i = 0; i < myLength; i++)
        {
                myIndex = i % numPts;
                //cout << myIndex << "\t\t";
                x =  myIndex*h + h/2;
                //cout << x << endl;
                temp1 = 1e3 * exp(-1e8 * pow(x - OffSet, 2 ) );
                //cout << temp1 << endl;
                //realtype temp2 = pow(1-x, 3.0/2.0);
                data[i] = temp1;
        }
}



void Constant(realtype * data, int vecLength)
{
	for(int i = 0; i < vecLength ; i ++)
		data[i] = 1;
}

void TestingInitCons(int SampleNum, int numPts, int numEqs, int vecLength, realtype Delx,
			realtype * data, realtype * StateData)
{
	if(SampleNum == 10)
	{
		SinWave(StateData, numPts, vecLength, Delx);
		for(int i = 0 ; i < numEqs; i++)
			data[i] = StateData[0];
	}
	if(SampleNum == 11)
		Constant(data, vecLength);

	if(SampleNum == 12)
	{
		SinWave2(StateData, numPts, vecLength, Delx);
		for( int i = 0 ; i < numEqs; i ++)
		{
			data[i] = StateData[0];
		}
	}

	if(SampleNum ==13)
	{
		GuassWave(StateData, numPts, vecLength, Delx);
		for( int i = 0 ; i < numEqs; i++)
		{
			data[i] = StateData[0];
		}
	}
}

void ConditionsFromFile(realtype * Ref, std :: string file)
{
        std :: string line;
        cout << "Getting Reference from: " << file << endl;
        //realtype * REF = N_VGetArrayPointer(Ref);
        //realtype Test = 0;
        std :: ifstream inputFile;
        inputFile.open(file);
        cout << "Is the file valid? " << (inputFile.good()? "Yes" : "No" ) << endl;
        int i = 0;
        while(!inputFile.eof())
        {
                getline(inputFile,line);
                Ref[i] = stod(line);          	//This is the problem line, shifts the solution by one.
                i++;
                inputFile.peek();               //Set the flag to see if the next entry is good.
        }
        inputFile.close();


}
