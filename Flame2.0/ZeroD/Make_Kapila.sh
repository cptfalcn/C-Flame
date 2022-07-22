echo ====================================Compiling Kapila.cpp================================
export TP=~/Code/C++/TCHEM3.0Serial
#================
#Build ZeroDVar.o
#================
#Tons of depreciation warnings, suppressed
/usr/bin/c++ -w -O3  \
-I /home/jstewart23/Code/C++/epic-cpp/Integrators/AdaptiveKrylov \
-I /home/jstewart23/Code/C++/epic-cpp/Integrators/EpiRK/ \
-I $TP/install/tchem/include/tchem \
-I $TP/install/kokkos/include \
-I $TP/install/tines/include/tines \
-I $TP/install/openblas/include \
-I /home/jstewart23/Code/C++/Sundials_6_2/include/nvector/ \
-I /home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/include \
-I /home/jstewart23/Code/C++/Flame2.0/Helpers \
-c /home/jstewart23/Code/C++/Flame2.0/ZeroD/Kapila.cpp -o Kapila.o

#link in the new file: /home/jstewart23/Code/C++/Sundials3_0/cvode-6.1.1/installDir
#Void
echo ============================================Linking==================================================
#=======================================
#Link all the helpers into OneDIgnVar.o
#=======================================
#==============
#Modern version
#==============
/usr/bin/c++ Kapila.o \
-O3 -w \
~/Code/C++/Flame2.0/Helpers/Epi3V.o \
$TP/install/tchem/lib/libtchem.a \
$TP/install/tines/lib/libtines.a \
$TP/install/kokkos/lib/libkokkoscontainers.a \
$TP/install/kokkos/lib/libkokkoscore.a \
/usr/lib/x86_64-linux-gnu/libdl.so \
$TP/install/yaml/lib/libyaml-cpp.a \
/home/jstewart23/Code/C++/epic-cpp/build/Integrators/libepic1.0.0.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_generic.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_nvecserial.a \
/home/jstewart23/Code/C++/Sundials_6_2/INSTDIR/lib/libsundials_cvode.a \
/home/jstewart23/Code/C++/TChem/build/install/openblas/lib/libopenblas.so \
-o Kapila.x

echo =======================================Removing intermediary objects=================================

#Clean up
rm Kapila.o
./Kapila.x 1e-1 "CVODE" "Kapila.txt" 1e-10 1e-5 
./Kapila.x 1e-1 "EPI3V" "Kapila.txt" 1e-10 1e-1