echo ====================================Compiling ZeroDIgnVar.cpp================================

export TP=~/Code/C++/TCHEM3.0Serial
#This build can be either
#export TP=~/Code/C++/TCHEM3.0
echo "building from $TP"
#================
#Build ZeroDCons.o
#================
#Tons of depreciation warnings, suppressed
/usr/bin/c++ -g -fopenmp -O3 -w  \
-I /home/jstewart23/Code/C++/epic-cpp/Integrators/AdaptiveKrylov \
-I /home/jstewart23/Code/C++/epic-cpp/Integrators/EpiRK/ \
-I $TP/install/tchem/include/tchem \
-I $TP/install/kokkos/include \
-I $TP/install/tines/include/tines \
-I $TP/install/openblas/include \
-I /home/jstewart23/Code/C++/Sundials_6_2/include/nvector/ \
-I /home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/include \
-I /home/jstewart23/Code/C++/Flame2.0/Helpers \
-c /home/jstewart23/Code/C++/Flame2.0/ZeroD/ZeroDIgnVar.cpp -o ZeroDIgnVar.o

echo ============================================Linking==================================================
#=======================================
#Link all the helpers into OneDIgnVar.o
#=======================================
#==============
#Modern version
#==============
/usr/bin/c++ ZeroDIgnVar.o \
-g -O3 -w -fopenmp -DNDEBUG \
~/Code/C++/Flame2.0/Helpers/Epi3VCntrl.o \
~/Code/C++/Flame2.0/Helpers/KiopsRetCode.o \
~/Code/C++/Flame2.0/Helpers/InitConditions.o \
~/Code/C++/Flame2.0/Helpers/Print.o \
~/Code/C++/Flame2.0/Helpers/Epi3V.o \
$TP/install/tchem/lib/libtchem.a \
$TP/install/tines/lib/libtines.a \
$TP/install/kokkos/lib/libkokkoscontainers.a \
$TP/install/kokkos/lib/libkokkoscore.a \
/home/jstewart23/Code/C++/epic-cpp/build/Integrators/libepic1.0.0.a \
/usr/lib/x86_64-linux-gnu/libdl.so \
$TP/install/yaml/lib/libyaml-cpp.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_generic.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_nvecserial.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_nvecopenmp.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_cvode.a \
$TP/install/openblas/lib/libopenblas.a \
-o ZeroDIgnVar.x
echo =======================================Removing intermediary objects=================================

#Clean up
rm ZeroDIgnVar.o

# inputs=gri3.0/

# #Use the default TChem names
# chemfile=$inputs"chem.inp"
# thermfile=$inputs"therm.dat"

# #Select an end time for the experiment, 10 for option 0 (Kapila), or something less than 1 otherwise (1e-4)
TF=1.3

# #Select the time integration method
# Method="EPI3V"

# #Select Krylov tolerances 
# Tol=1e-6

# #Set the output file name
# #OutputFile="Scrapt.txt"

# #Set Experiment
Experiment=2

# #Set Sample
Sam=1

# #Set Profiling. Add Profiling=1 if you want to see the profile.
# Prof=1

Out="ScrapVar.txt"
# #./ZeroDIgnVariable.x --chemfile=$chemfile --thermfile=$thermfile --StepSize=1e-8 --FinalTime=1e-5 --MyFile=$OutputFile  --KrylovTol=$Tol --UseJac=1 --SampleNum=$Sam --Method="EPI2" --Experiment=2 --Profiling=$Prof
# TF=1.3e0 
# Out="Scrap.txt"
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method="CVODEKry" --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-9 --relTol=1e-8 --Movie=0 --Input=gri3.0/sample.dat

./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method="EPI3V" --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=1e-1 --Movie=0 --Input=gri3.0/sample.dat
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method="EPI3VCntrl" --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=1e-1 --Movie=0 --Input=gri3.0/sample.dat
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method="EPI3V" --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=5e-2 --Movie=0 --Input=gri3.0/sample.dat
echo ============================================Testing 10 atm as input =======================================
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method="EPI3V" --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=5e-2 --Movie=0 --Input=gri3.0/sample.dat --PressMult=10
echo ============================================Testing JacFile output =======================================
./ZeroDIgnVar.x --chemfile=gri3.0/"chem.inp" --thermfile=gri3.0/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=$Sam --Method="EPI3V" --Profiling=1 --Experiment=2 --maxSS=5e-1 --absTol=1e-10 --relTol=5e-2 --Movie=1 --Input=gri3.0/sample.dat --PressMult=10 --JacFile="TestJacDump.txt"

#Both these stall and will not run.  But TChem TrBDF2 will run this sample okay.
TF=5e-2
#./ZeroDIgnVar.x --chemfile=iso-oct/"chem.inp" --thermfile=iso-oct/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-5 --UseJac=1 --SampleNum=1 --Method="CVODEKry" --Profiling=1 --Experiment=3 --maxSS=1e-2 --absTol=1e-10 --relTol=1e-10 --Movie=0 --Input=iso-oct/"sample1_10atm.dat"
#./ZeroDIgnVar.x --chemfile=iso-oct/"chem.inp" --thermfile=iso-oct/"therm.dat" --StepSize=1e-5 --FinalTime=$TF --MyFile=$Out --KrylovTol=1e-7 --UseJac=1 --SampleNum=1 --Method="EPI3V" --Profiling=1 --Experiment=3 --maxSS=1e-2 --absTol=1e-10 --relTol=1e-1 --Movie=0 --Input=iso-oct/"sample1_10atm.dat"