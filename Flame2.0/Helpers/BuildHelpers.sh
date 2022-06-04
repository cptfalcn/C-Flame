#compile helper functions
echo ========================================Building helpers=============================================
#/usr/bin/c++ -I ~/Code/C++/epic-cpp/Integrators/AdaptiveKrylov/ -I ~/Code/C++/epic-cpp/Integrators/EpiRK -c ~/Code/C++/epic-cpp/Integrators/EpiRK/Epi3_KIOPS.cpp -o trash.o

#/usr/bin/c++ -I ~/Code/C++/epic-cpp/Integrators/AdaptiveKrylov/ -I ~/Code/C++/epic-cpp/Integrators/EpiRK -c ~/Code/C++/epic-cpp/Integrators/EpiRK/Epi3SC_KIOPS.cpp -o Epi3SC.o

#/usr/bin/c++ -I ~/Code/C++/epic-cpp/Integrators/AdaptiveKrylov/ -I ~/Code/C++/epic-cpp/Integrators/EpiRK -c ~/Code/C++/epic-cpp/Integrators/EpiRK/Epi3Ros.cpp -o Epi3Ros.o

#/usr/bin/c++ -DKOKKOS_DEPENDENCE -I/home/jstewart23/Code/C++/TChem/build/build/tchem -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core/impl -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core/problems -isystem /home/jstewart23/Code/C++/TChem/build/install/kokkos/include -isystem /home/jstewart23/Code/C++/TChem/build/install/tines/include/tines -isystem /home/jstewart23/Code/C++/TChem/build/install/openblas/include -I ~/Code/C++/Flame/OneD/ -I ~/Code/C++/epic-cpp/Integrators/AdaptiveKrylov/ -I ~/Code/C++/epic-cpp/Integrators/EpiRK -c ~/Code/C++/epic-cpp/Integrators/EpiRK/EpiP2_KIOPS.cpp -o EpiP2.o

/usr/bin/c++  -g -pg -O0 \
-I /home/jstewart23/Code/C++/Sundials_6_2/include \
-c /home/jstewart23/Code/C++/Flame/OneD/InitialConditions.cpp \
-o InitConditions.o

#/usr/bin/c++  -g -pg -O0 -c /home/jstewart23/Code/C++/Flame/OneD/Print.cpp -I/home/jstewart23/Code/C++/epic-cpp/Integrators/AdaptiveKrylov -I/home/jstewart23/Code/C++/epic-cpp/Integrators/EpiRK -o Print.o

#/usr/bin/c++  -g -pg -O0 -DKOKKOS_DEPENDENCE -I/home/jstewart23/Code/C++/TChem/build/build/tchem -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core/impl -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core/problems -isystem /home/jstewart23/Code/C++/TChem/build/install/kokkos/include -isystem /home/jstewart23/Code/C++/TChem/build/install/tines/include/tines -isystem /home/jstewart23/Code/C++/TChem/build/install/openblas/include -DNDEBUG -fopenmp -I/home/jstewart23/Code/C++/epic-cpp/Integrators/AdaptiveKrylov -I/home/jstewart23/Code/C++/epic-cpp/Integrators/EpiRK -c /home/jstewart23/Code/C++/Flame/OneD/Chemistry.cpp -o Chemistry.o

#/usr/bin/c++  -g -pg -O0 -DKOKKOS_DEPENDENCE -I/home/jstewart23/Code/C++/TChem/build/build/tchem -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core/impl -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core/problems -isystem /home/jstewart23/Code/C++/TChem/build/install/kokkos/include -isystem /home/jstewart23/Code/C++/TChem/build/install/tines/include/tines -isystem /home/jstewart23/Code/C++/TChem/build/install/openblas/include -DNDEBUG -fopenmp -I/home/jstewart23/Code/C++/epic-cpp/Integrators/AdaptiveKrylov -I/home/jstewart23/Code/C++/epic-cpp/Integrators/EpiRK -c /home/jstewart23/Code/C++/Flame/OneD/CreateIntegrators.cpp -o Integrators.o

echo =========================================Helpers built===============================================
