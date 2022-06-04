echo ====================================Compiling TestZeroDRHSJac.cpp================================
export TP=~/Code/C++/TCHEM3.0
#================
#Build ZeroDVar.o
#================
#Tons of depreciation warnings, suppressed
/usr/bin/c++ -g -fopenmp -O0  \
-I $TP/install/tchem/include/tchem \
-I $TP/install/kokkos/include \
-I $TP/install/tines/include/tines \
-I $TP/install/openblas/include \
-I /home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/include \
-c /home/jstewart23/Code/C++/Flame2.0/ZeroD/TestZeroDSundials.cpp -o ZeroDTest.o
# /usr/bin/c++ -g -O3 \
# -I /home/jstewart23/Code/C++/TChem/build/build/tchem \
# -I /home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core \
# -I /home/jstewart23/Code/C++/Flame/OneD \
# -I /home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core/impl \
# -I /home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core/problems \
# -I /home/jstewart23/Code/C++/TChem/build/install/kokkos/include \
# -I /home/jstewart23/Code/C++/TChem/build/install/tines/include/tines \
# -I /home/jstewart23/Code/C++/TChem/build/install/openblas/include \
# -I /home/jstewart23/Code/C++/epic-cpp/Integrators/AdaptiveKrylov \
# -I /home/jstewart23/Code/C++/epic-cpp/Integrators/EpiRK \
# -c /home/jstewart23/Code/C++/Flame/ZeroD/ZeroDIgnVariable.cpp -o ZeroDVars.o

#link in the new file: /home/jstewart23/Code/C++/Sundials3_0/cvode-6.1.1/installDir
#Void
echo ============================================Linking==================================================
#=======================================
#Link all the helpers into OneDIgnVar.o
#=======================================
#==============
#Modern version
#==============
/usr/bin/c++ ZeroDTest.o \
-g -O0 -w -fopenmp -DNDEBUG \
$TP/install/tchem/lib/libtchem.a \
$TP/install/tines/lib/libtines.a \
$TP/install/kokkos/lib/libkokkoscontainers.a \
$TP/install/kokkos/lib/libkokkoscore.a \
/usr/lib/x86_64-linux-gnu/libdl.so \
$TP/install/yaml/lib/libyaml-cpp.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_generic.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_nvecserial.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_cvode.a \
/home/jstewart23/Code/C++/TChem/build/install/openblas/lib/libopenblas.so \
-o TestZeroDSundials.x

# /usr/bin/c++ Epi3SC.o  ZeroDVars.o trash.o InitConditions.o \
# Integrators.o Chemistry.o Print.o -g -O3 -DNDEBUG \
# ~/Code/C++/TChem/build/build/tchem/core/libtchem.a \
# /home/jstewart23/Code/C++/TChem/build/install/tines/lib/libtines.a \
# /home/jstewart23/Code/C++/TChem/build/install/kokkos/lib/libkokkoscontainers.a \
# /home/jstewart23/Code/C++/TChem/build/install/kokkos/lib/libkokkoscore.a \
# /usr/local/lib/libsundials_cvode.so \
# /usr/local/lib/libsundials_nvecserial.so \
# /home/jstewart23/Code/C++/epic-cpp/build/Integrators/libepic1.0.0.a \
# /usr/lib/x86_64-linux-gnu/libdl.so \
# /home/jstewart23/Code/C++/TChem/build/install/openblas/lib/libopenblas.so \
# -o ZeroDIgnVariable.x

#/usr/bin/c++  EpiP2.o Epi3SC.o Integrators.o ZeroDVar2.o trash.o Print.o Derivatives.o InitConditions.o Chemistry.o -pg -O0 -DNDEBUG -lgfortran -rdynamic -DKOKKOS_DEPENDENCE -fopenmp -Wl,-rpath,$TCHEM_PATH/openblas/lib: ~/Code/C++/TChem/build/build/tchem/core/libtchem.a $TCHEM_PATH/tines/lib/libtines.a $TCHEM_PATH/gtest/lib/libgtest_main.a $TCHEM_PATH/kokkos/lib/libkokkoscontainers.a $TCHEM_PATH/kokkos/lib/libkokkoscore.a /usr/lib/x86_64-linux-gnu/libdl.so $TCHEM_PATH/openblas/lib/libopenblas.so $TCHEM_PATH/gtest/lib/libgtest.a -lpthread -rdynamic -pthread -Wl,-rpath,/usr/local/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib -lopenblas /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so $EPIC_PATH/Integrators/libepic1.0.0.a -lopenblas $LIB_PATH/libsundials_cvode.so $LIB_PATH/libsundials_cvode.so.6 $LIB_PATH/libsundials_nvecparallel.so /usr/local/lib/libsundials_nvecserial.so /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so -o ZeroDIgnVar2.x
echo =======================================Removing intermediary objects=================================

#Clean up
rm ZeroDTest.o
#rm ZeroDVar2.o
# echo ===========================================Built ZeroDIgnCons.x=======================================

# inputs=gri3.0/

# #Use the default TChem names
# chemfile=$inputs"chem.inp"
# thermfile=$inputs"therm.dat"

# #Select an end time for the experiment, 10 for option 0 (Kapila), or something less than 1 otherwise (1e-4)
# EndT=1e-4

# #Select the time integration method
# Method="EPI3V"

# #Select Krylov tolerances 
# Tol=1e-6

# #Set the output file name
# #OutputFile="Scrapt.txt"

# #Set Experiment
# Experiment=2

# #Set Sample
# Sam=1

# #Set Profiling. Add Profiling=1 if you want to see the profile.
# Prof=1

# OutputFile="Scrap.txt"
# #./ZeroDIgnVariable.x --chemfile=$chemfile --thermfile=$thermfile --StepSize=1e-8 --FinalTime=1e-5 --MyFile=$OutputFile  --KrylovTol=$Tol --UseJac=1 --SampleNum=$Sam --Method="EPI2" --Experiment=2 --Profiling=$Prof
# TF=1.3e0 
# Out="Scrap.txt"
# ./ZeroDIgnVariable.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method="CVODEKry" --Profiling=1 --Experiment=2 --maxSS=5e-2 --absTol=1e-9 --relTol=1e-8 --Movie=0
# ./ZeroDIgnVariable.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=7.5e-2 --Movie=0
# ./ZeroDIgnVariable.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=7.5e-2 --Movie=0
# #./ZeroDIgnVariable.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=2.18e-2 --Movie=0
# #./ZeroDIgnVariable.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-4 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method=$Method --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=3.13e-3 --Movie=0
./TestZeroDSundials.x 1e-1 "CVODE" "Kapila.txt" 
