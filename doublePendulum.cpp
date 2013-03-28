#include "ODESolution.h"
#include <sstream>
using namespace std;

int main()
{
  int nSteps = 1000000;
  double h = 0.0001;
  double G = 1.;
  double R = 1.;
  double theta;
  double omega;
  cout << "G: ";
  cin >> G;
  cout << "R: ";
  cin >> R;
  cout << "h: ";
  cin >> h;
  cout << "n: ";
  cin >> nSteps;
  cout << "theta: ";
  cin >> theta;
  cout << "omega: ";
  cin >> omega;

  std::ostringstream dataOut;
  dataOut << "./data/double/" << "RungeKutta4_" << h << "_" << nSteps << "_" << G << "_" << R << "_theta_" << theta << "_omega_" << omega  << ".txt";
  std::string dataOutDir = dataOut.str();

  std::vector<double> init(4);
  init[0] = theta;
  init[1] = 0.;
  init[2] = omega;
  init[3] = 0.;

  RungeKuttaFour a(dataOutDir,init,h,nSteps,G,R);
  a.solve();
  
  
  return 0;
}
