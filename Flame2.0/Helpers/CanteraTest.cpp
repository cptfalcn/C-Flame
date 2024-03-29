#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include <iostream>

using namespace Cantera;

// The actual code is put into a function that can be called from the main program.
void simple_demo()
{
    // Create a new Solution object
    auto sol = newSolution("h2o2.yaml", "ohmech", "None");
    auto gas = sol->thermo();

    // Set the thermodynamic state by specifying T (500 K) P (2 atm) and the mole
    // fractions. Note that the mole fractions do not need to sum to 1.0 - they will
    // be normalized internally. Also, the values for any unspecified species will be
    // set to zero.
    gas->setState_TPX(500.0, 2.0*OneAtm, "H2O:1.0, H2:8.0, AR:1.0");

    // Print a summary report of the state of the gas.
    std::cout << gas->report() << std::endl;
}

// The main program just calls function simple_demo within a 'try' block, and catches
// CanteraError exceptions that might be thrown.
int main()
{
    try {
        simple_demo();
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
    }
}

// #include "cantera/base/Array.h"
// #include "cantera/base/plots.h"

// // Save the temperature, density, pressure, and mole fractions at one
// // time
// template<class G, class A>
// void saveSoln(int i, double time, const G& gas, A& soln)
// {
//     soln(0,i) = time;
//     soln(1,i) = gas.temperature();
//     soln(2,i) = gas.density();
//     soln(3,i) = gas.pressure();
//     gas.getMoleFractions(&soln(4,i));
// }

// template<class G, class A>
// void saveSoln(double time, const G& gas, A& soln)
// {
//     soln.resize(soln.nRows(), soln.nColumns() + 1);
//     int back = soln.nColumns() - 1;
//     soln(0,back) = time;
//     soln(1,back) = gas.temperature();
//     soln(2,back) = gas.density();
//     soln(3,back) = gas.pressure();
//     int nsp = gas.nSpecies();
//     for (int k = 0; k < nsp; k++) {
//         soln(4+k,back) = gas.moleFraction(k);
//     }
// }

// template<class G, class V>
// void makeDataLabels(const G& gas, V& names)
// {
//     int nsp = gas.nSpecies();
//     names.resize(nsp + 4);
//     names[0] = "time (s)";
//     names[1] = "Temperature (K)";
//     names[2] = "Density (kg/m3)";
//     names[3] = "Pressure (Pa)";
//     int k;
//     for (k = 0; k < nsp; k++) {
//         names[4+k] = gas.speciesName(k);
//     }
// }

// template<class G, class A>
// void plotSoln(const std::string& fname, const std::string& fmt,
//               const std::string& title, const G& gas, const A& soln)
// {
//     std::vector<std::string> names;
//     makeDataLabels(gas, names);
//     writePlotFile(fname, fmt, title, names, soln);
// }


// /*!
//  * @file kinetics1.cpp
//  *
//  * Zero-dimensional kinetics
//  *
//  * This example simulates autoignition of hydrogen in a constant pressure
//  * reactor and saves the time history to files that can be used for plotting.
//  *
//  * Keywords: combustion, reactor network, ignition delay, saving output
//  */

// // This file is part of Cantera. See License.txt in the top-level directory or
// // at https://cantera.org/license.txt for license and copyright information.

// #include "cantera/zerodim.h"
// #include "cantera/thermo/IdealGasPhase.h"
// #include "cantera/numerics/Integrator.h"
// //#include "example_utils.h"

// using namespace Cantera;
// using std::cout;
// using std::endl;

// int kinetics1(int np, void* p)
// {
//     cout << "Constant-pressure ignition of a "
//          << "hydrogen/oxygen/nitrogen"
//          " mixture \nbeginning at T = 1001 K and P = 1 atm." << endl;

//     // create an ideal gas mixture that corresponds to OH submech from GRI-Mech 3.0
//     auto sol = newSolution("h2o2.yaml", "ohmech", "None");
//     auto gas = sol->thermo();

//     // set the state
//     gas->setState_TPX(1001.0, OneAtm, "H2:2.0, O2:1.0, N2:4.0");
//     int nsp = gas->nSpecies();

//     // create a reactor
//     IdealGasConstPressureReactor r;

//     // 'insert' the gas into the reactor and environment.  Note
//     // that it is ok to insert the same gas object into multiple
//     // reactors or reservoirs. All this means is that this object
//     // will be used to evaluate thermodynamic or kinetic
//     // quantities needed.
//     r.insert(sol);

//     double dt = 1.e-5; // interval at which output is written
//     int nsteps = 100; // number of intervals

//     // create a 2D array to hold the output variables,
//     // and store the values for the initial state
//     Array2D soln(nsp+4, 1);
//     saveSoln(0, 0.0, *(sol->thermo()), soln);

//     // create a container object to run the simulation
//     // and add the reactor to it
//     ReactorNet sim;
//     sim.addReactor(r);

//     // main loop
//     clock_t t0 = clock(); // save start time
//     for (int i = 1; i <= nsteps; i++) {
//         double tm = i*dt;
//         sim.advance(tm);
//         cout << "time = " << tm << " s" << endl;
//         saveSoln(tm, *(sol->thermo()), soln);
//     }
//     clock_t t1 = clock(); // save end time


//     // make a Tecplot data file and an Excel spreadsheet
//     std::string plotTitle = "kinetics example 1: constant-pressure ignition";
//     plotSoln("kin1.dat", "TEC", plotTitle, *(sol->thermo()), soln);
//     plotSoln("kin1.csv", "XL", plotTitle, *(sol->thermo()), soln);


//     // print final temperature and timing data
//     double tmm = 1.0*(t1 - t0)/CLOCKS_PER_SEC;
//     cout << " Tfinal = " << r.temperature() << endl;
//     cout << " time = " << tmm << endl;
//     cout << " number of residual function evaluations = "
//          << sim.integrator().nEvals() << endl;
//     cout << " time per evaluation = " << tmm/sim.integrator().nEvals()
//          << endl << endl;
//     cout << "Output files:" << endl
//          << "  kin1.csv    (Excel CSV file)" << endl
//          << "  kin1.dat    (Tecplot data file)" << endl;

//     return 0;
// }


// int main()
// {
//     try {
//         int retn = kinetics1(0, 0);
//         appdelete();
//         return retn;
//     } catch (CanteraError& err) {
//         // handle exceptions thrown by Cantera
//         std::cout << err.what() << std::endl;
//         cout << " terminating... " << endl;
//         appdelete();
//         return -1;
//     }
// }
