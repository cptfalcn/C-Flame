#cmake -S . -B ./BuildProject/
cmake_minimum_required(VERSION 3.1)
project (OneDConsIgn)

#include ("/usr/include/cantera/Cantera.mak")
set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CMAKE_CXX_STANDARD 19)

#find_package(Threads REQUIRED)

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

#link_directories("/usr/lib")

add_compile_options(-g -pg -O2 -fopenmp)
add_link_options(-pg -fopenmp -O3 -Wall -DNDEBUG)

add_executable(OneDIgnCons.x /home/jstewart23/Code/C++/Flame2.0/OneDOldCantera/OneDIgnCons.cpp)

target_link_libraries(OneDIgnCons.x cantera)
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/Flame2.0/OneDOld/Epi3SC.o")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/Flame2.0/OneDOld/Integrators.o")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/Flame2.0/OneDOld/trash.o")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/Flame2.0/OneDOld/Print.o")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/Flame2.0/OneDOld/InitConditions.o")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/Flame2.0/OneDOld/Chemistry.o")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/TChem/build/build/tchem/core/libtchem.a")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/TChem/build/install/tines/lib/libtines.a")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/TChem/build/install/kokkos/lib/libkokkoscontainers.a")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/TChem/build/install/kokkos/lib/libkokkoscore.a")
target_link_libraries(OneDIgnCons.x "/usr/local/lib/libsundials_cvode.so")
target_link_libraries(OneDIgnCons.x "/usr/local/lib/libsundials_nvecserial.so")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/epic-cpp/build/Integrators/libepic1.0.0.a")
target_link_libraries(OneDIgnCons.x "/usr/lib/x86_64-linux-gnu/libdl.so")
target_link_libraries(OneDIgnCons.x "/home/jstewart23/Code/C++/TChem/build/install/openblas/lib/libopenblas.so")
# fmt yaml-cpp lapack blas Threads::Threads)
