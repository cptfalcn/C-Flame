#include <iostream>
#include <Eigen>
#include <chrono>

using Eigen::MatrixXd;
 
int main()
{

  MatrixXd m = MatrixXd::Random(54,54);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  //std::cout << m << std::endl;
  auto Start=std::chrono::high_resolution_clock::now();//Time integrator
  //std :: cout << m.eigenvalues() << std :: endl;
  m.eigenvalues();
	auto Stop=std::chrono::high_resolution_clock::now();
  auto Pass = std::chrono::duration_cast<std::chrono::nanoseconds>(Stop-Start);
  std :: cout << "Eigenvalue time: " << Pass.count()/1e9 << "seconds\n";

}


