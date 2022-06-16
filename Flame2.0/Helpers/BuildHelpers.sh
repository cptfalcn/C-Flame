#compile helper functions
echo ========================================Building helpers=============================================
export TP=~/Code/C++/TCHEM3.0

#EPIP2
#/usr/bin/c++ -DKOKKOS_DEPENDENCE -I/home/jstewart23/Code/C++/TChem/build/build/tchem -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core/impl -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core/problems -isystem /home/jstewart23/Code/C++/TChem/build/install/kokkos/include -isystem /home/jstewart23/Code/C++/TChem/build/install/tines/include/tines -isystem /home/jstewart23/Code/C++/TChem/build/install/openblas/include -I ~/Code/C++/Flame/OneD/ -I ~/Code/C++/epic-cpp/Integrators/AdaptiveKrylov/ -I ~/Code/C++/epic-cpp/Integrators/EpiRK -c ~/Code/C++/epic-cpp/Integrators/EpiRK/EpiP2_KIOPS.cpp -o EpiP2.o

#Chem Initial Conditions
/usr/bin/c++  -g -pg -O3 \
-I /home/jstewart23/Code/C++/Sundials_6_2/include \
-c /home/jstewart23/Code/C++/Flame2.0/Helpers/InitialConditions.cpp \
-o InitConditions.o

#General Print functions
/usr/bin/c++  -g -pg -O3 \
-I /home/jstewart23/Code/C++/Sundials_6_2/include \
-c /home/jstewart23/Code/C++/Flame2.0/Helpers/Print.cpp \
-o Print.o

#OneD Problem Object
#/usr/bin/c++  -g -pg -O0 -DKOKKOS_DEPENDENCE -I/home/jstewart23/Code/C++/TChem/build/build/tchem -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core/impl -I/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core/problems -isystem /home/jstewart23/Code/C++/TChem/build/install/kokkos/include -isystem /home/jstewart23/Code/C++/TChem/build/install/tines/include/tines -isystem /home/jstewart23/Code/C++/TChem/build/install/openblas/include -DNDEBUG -fopenmp -I/home/jstewart23/Code/C++/epic-cpp/Integrators/AdaptiveKrylov -I/home/jstewart23/Code/C++/epic-cpp/Integrators/EpiRK -c /home/jstewart23/Code/C++/Flame/OneD/Chemistry.cpp -o Chemistry.o


#EPI3V
/usr/bin/c++  -g -pg -O3 -w -DKOKKOS_DEPENDENCE \
-I $TP/install/tchem/include/tchem \
-I $TP/install/kokkos/include \
-I $TP/install/tines/include/tines \
-I $TP/install/openblas/include \
-I /home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/include \
-DNDEBUG -fopenmp \
-I/home/jstewart23/Code/C++/epic-cpp/Integrators/AdaptiveKrylov \
-I/home/jstewart23/Code/C++/epic-cpp/Integrators/EpiRK \
-c  /home/jstewart23/Code/C++/Flame2.0/Helpers/Epi3V.cpp -o Epi3V.o

echo =========================================Helpers built===============================================
