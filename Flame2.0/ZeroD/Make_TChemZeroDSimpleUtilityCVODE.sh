echo ====================================Compiling TChemZeroDUtility.cpp================================
export TP=~/Code/C++/TCHEM3.0CVODE
export Base=~/Code/C++/TCHEM3.0
echo "Building with $TP"
#================
#Build ZeroDVar.o
#================
#Tons of depreciation warnings, suppressed
/usr/bin/c++ -g -fopenmp -O3  -w \
-I $TP/install/tchem/include/tchem \
-I $TP/install/kokkos/include \
-I $TP/install/tines/include/tines \
-I $Base/install/openblas/include \
-I /home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/include \
-c /home/jstewart23/Code/C++/Flame2.0/ZeroD/TChemZeroDSimpleUtilityCVODE.cpp -o TChemZeroDSimpleUtility.o

echo ============================================Linking==================================================
#=======================================
#Link all the helpers into OneDIgnVar.o
#=======================================
#==============
#Modern version
#==============
/usr/bin/c++ TChemZeroDSimpleUtility.o \
-g -O3 -w -fopenmp \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_generic.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_nvecserial.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_cvode.a \
$TP/install/tchem/lib/libtchem.a \
$TP/install/tines/lib/libtines.a \
$TP/install/kokkos/lib/libkokkoscontainers.a \
$TP/install/kokkos/lib/libkokkoscore.a \
$Base/install/yaml/lib/libyaml-cpp.a \
$Base/install/openblas/lib/libopenblas.a \
/usr/lib/x86_64-linux-gnu/libdl.so \
/usr/lib/x86_64-linux-gnu/libgfortran.so.5.0.0 \
-o TChemZeroDSimpleUtilityCVODE.x
#$TP/install/openblas/lib/libopenblas.a \

echo =======================================Removing intermediary objects=================================

#Clean up
rm TChemZeroDSimpleUtility.o
#rm ZeroDVar2.o
# echo ===========================================Built ZeroDIgnCons.x=======================================
#./Run_TChemZeroDUtility.sh
# Runs a sample
# inputs=iso-oct
# #Default settings, build into a new utilty
# this="$exec --chemfile=$inputs/chem.inp \
#             --thermfile=$inputs/therm.dat \
#             --samplefile=$inputs/sample.dat \
#             --outputfile=IgnSolution.dat \
#             --atol-newton=1e-18 \
#             --rtol-newton=1e-8\
#             --max-newton-iterations=20 \
#             --tol-time=1e-6 \
#             --dtmax=1e-3 \
#             --dtmin=1e-20 \
#             --tend=2 \
#             --time-iterations-per-interval=10 \
#             --max-time-iterations=260 \
#             --ignition-delay-time-file=IgnitionDelayTime.dat \
#             --ignition-delay-time-w-threshold-temperature-file=IgnitionDelayTimeTthreshold.dat
#             --threshold-temperature=1500"
# #echo $this
# eval $this