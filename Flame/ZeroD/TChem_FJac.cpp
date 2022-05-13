/* =====================================================================================
TChem version 2.0
Copyright (2020) NTESS
https://github.com/sandialabs/TChem

Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

This file is part of TChem. TChem is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Cosmin Safta at <csafta@sandia.gov>, or
           Kyungjoo Kim at <kyukim@sandia.gov>, or
           Oscar Diaz-Ibarra at <odiazib@sandia.gov>

Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */
#include "TChem_CommandLineParser.hpp"

#include "TChem_Util.hpp"
#include "TChem_KineticModelData.hpp"
#include "TChem_Impl_IgnitionZeroD_Problem.hpp" // here is where Ignition Zero D problem is implemented
#include <iostream>

void IntConGri30(real_type * data, int SampleNum);


int main(int argc, char* argv[])
{
  ///
  /// 1. Command line input parser
  ///

  /// default input filename
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

  ///
  /// 2. Kokkos initialization
  ///

  /// Kokkos initialization requires local scoping to properly deallocate global variables.
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

    ///
    /// 4. Construction of a kinetic model
    ///

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
    //x_std[0] = 1200; // temperature
    //x_std[1] = 0.1;  // mass fraction species 1
    //x_std[2] = 0.9;  // mass fraction species 2
   

    /// output rhs vector and J matrix
    std::vector<real_type> rhs_std(number_of_equations);
    std::vector<real_type> J_std(number_of_equations * number_of_equations);

    /// get points for std vectors
    real_type * ptr_x = x_std.data();
    real_type * ptr_rhs = rhs_std.data();
    real_type * ptr_J = J_std.data();
    IntConGri30(ptr_x, 1);

    ///
    /// 7. Compute right hand side vector and Jacobian
    ///
    {
      problem_type problem;

      /// initialize problem
      problem._p = pressure; // pressure
      problem._work = work;  // problem workspace array
      problem._kmcd = kmcd;  // kinetic model

      /// create view wrapping the pointers
      real_type_1d_view_type x(ptr_x,   number_of_equations);
      real_type_1d_view_type f(ptr_rhs, number_of_equations);
      real_type_2d_view_type J(ptr_J,   number_of_equations, number_of_equations);

      /// create a fake team member
      const auto member = Tines::HostSerialTeamMember();

      /// compute rhs and Jacobian
      problem.computeFunction(member, x, f);
      problem.computeJacobian(member, x, J);

      /// change the layout from row major to column major
      for (ordinal_type j=0;j<number_of_equations;++j)
	      for (ordinal_type i=0;i<j;++i)
	        std::swap(J(i,j), J(j,i));
    }

    ///
    /// 8. Check the rhs and Jacobian
    ///

    printf("RHS std vector \n" );
    for (ordinal_type i=0;i<number_of_equations;++i) {
      printf("%e\n", rhs_std[i]);
	std::cout<<"\n";
    }

    printf("Jacobian std vector \n" );
    for (ordinal_type i=0;i<number_of_equations;++i) {
      for (ordinal_type j=0;j<number_of_equations;++j) {
        printf("%e ", J_std[i+j*number_of_equations] );
      }
      printf("\n" );
    }

    /// all Kokkos varialbes are reference counted objects. they are deallocated
    /// within this local scope.
  }

  /// Kokkos finalize checks any memory leak that are not properly deallocated.
  Kokkos::finalize();

  return 0;
}

void IntConGri30(real_type * data, int SampleNum)
{
        
        #define BAR "===================="
        data[0] = 1000.0;//use 1000 for base
        switch(SampleNum)
        {
        case 1:
                data[0] = 1000.0;//use 1000 for base
                data[14] =  5.483858845025559731e-02;//CH4
                data[4]  =  2.187578062376045740e-01;//O2
                data[48] =  7.137587863547695255e-01;//N2
                data[49] =  1.264481895737025810e-02;//AR
                std::cout << BAR <<"GRI Experiment Sample 1\t" << BAR <<std::endl;
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
		std::cout << BAR << "Gri experiment, Sample 1: Ignition" << BAR <<std::endl;
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
		std::cout << BAR << "Gri Sample 1: Post Ign" << BAR << std::endl;
		break;
	case 4: //Aditya's case, runs at three ATM
		data[0] 	=  800;//Temp
                data[14]	=  0.0392661; //CH4
                data[4] 	=  0.2237715;//O2
                data[48]	=  0.736962400000000; //N2
		std::cout << BAR << "Gri Sample: Hi Pressure" << BAR << std::endl;
		break;
	case 5: //High temp setup
		data[0]  =  1300;
		data[14] =  5.483858845025559731e-02;//CH4
                data[4]  =  2.187578062376045740e-01;//O2
                data[48] =  7.137587863547695255e-01;//N2
                data[49] =  1.264481895737025810e-02;//AR
                std::cout << BAR <<"GRI Experiment Sample 2\t" << BAR <<std::endl;
		break;
	case 6: //Low temp long burn
		std::cout << BAR << "GRI Experiment: Long burn" << BAR << std::endl;
		data[0] =   810.0;//use 1000 for base
                data[14] =  5.483858845025559731e-02;//CH4
                data[4]  =  2.187578062376045740e-01;//O2
                data[48] =  7.137587863547695255e-01;//N2
                data[49] =  1.264481895737025810e-02;//AR
                std::cout << BAR <<"GRI Experiment Sample 1\t" << BAR <<std::endl;
                break;
	// case 69:
	// 	std::cout << BAR << "Gri Sample: 1e-4 Hi-P" << BAR << std::endl;
	// 	ConditionsFromFile(data, "GriIgnition.txt");
	// 	break;
        }
}
