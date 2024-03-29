#cmake -S . -B ./BuildProject/
cmake_minimum_required(VERSION 3.1)
project (OneDConsIgn)

#include ("/usr/include/cantera/Cantera.mak")
set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CMAKE_CXX_STANDARD 19)

find_package(Threads REQUIRED)

include_directories("/home/jstewart23/Code/C++/TChem/build/build/tchem\
" "/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core\
" "/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core/impl\
" "/home/jstewart23/Code/C++/TChem/build/repos/tchem/src/core/problems\
" "/home/jstewart23/Code/C++/TChem/build/install/kokkos/include\
" "/home/jstewart23/Code/C++/TChem/build/install/tines/include/tines\
" "/home/jstewart23/Code/C++/TChem/build/install/openblas/include\
" "/home/jstewart23/Code/C++/epic-cpp/Integrators/AdaptiveKrylov\
" "/home/jstewart23/Code/C++/epic-cpp/Integrators/EpiRK\
" "home/jstewart23/Code/C++/Flame2.0/OneDOldCantera")
include_directories(-lfmt)
include_directories( "/opt/cantera/include")

#link_directories("/usr/lib")

add_compile_options(-g -pg -O2 -fopenmp)
add_link_options(-pg -fopenmp -O3 -Wall -DNDEBUG)

add_executable(OneDIgnCons.x /home/jstewart23/Code/C++/Flame2.0/OneDOldCantera/OneDIgnCons.cpp)

target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/Flame2.0/OneDOldCantera/Epi3SC.o")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/Flame2.0/OneDOldCantera/Integrators.o")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/Flame2.0/OneDOldCantera/trash.o")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/Flame2.0/OneDOldCantera/Print.o")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/Flame2.0/OneDOld/InitConditions.o")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/Flame2.0/OneDOldCantera/Chemistry.o")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/TChem/build/build/tchem/core/libtchem.a")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/TChem/build/install/tines/lib/libtines.a")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/TChem/build/install/kokkos/lib/libkokkoscontainers.a")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/TChem/build/install/kokkos/lib/libkokkoscore.a")
target_link_libraries(OneDIgnCons.x "/usr/local/lib/libsundials_cvode.so")
target_link_libraries(OneDIgnCons.x "/usr/local/lib/libsundials_nvecserial.so")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/epic-cpp/build/Integrators/libepic1.0.0.a")
target_link_libraries(OneDIgnCons.x "/usr/lib/x86_64-linux-gnu/libdl.so")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/TChem/build/install/openblas/lib/libopenblas.so")

#Added to link Cantera
target_link_libraries(OneDIgnCons.x cantera sundials_cvodes sundials_ida sundials_nvecserial sundials_sunlinsollapackdense sundials_sunlinsollapackband fmt yaml-cpp lapack blas Threads::Threads)


#Sample runnable CMAKE
#find_package(Threads REQUIRED)
#include_directories("/usr/include" "/usr/include/eigen3")
#link_directories("/usr/lib")
#add_executable(demo CanteraTest3.cpp)
#target_link_libraries(demo cantera sundials_cvodes sundials_ida sundials_nvecserial sundials_sunlinsollapackdense sundials_sunlinsollapackband fmt yaml-cpp lapack blas Threads::Threads)
