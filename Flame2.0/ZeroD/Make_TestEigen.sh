/usr/bin/c++ -g -O3 -w  \
-I /home/jstewart23/Code/C++/Flame2.0/Helpers/eigen/Eigen \
-c /home/jstewart23/Code/C++/Flame2.0/ZeroD/TestEigen.cpp -o TestEigen.o

/usr/bin/c++ -g -O3 -w  \
TestEigen.o -o TestEigen.x
echo "Compiled"

./TestEigen.x