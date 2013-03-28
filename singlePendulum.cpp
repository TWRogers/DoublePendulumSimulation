#include "ODESolution.h"
#include <sstream>
using namespace std;

int main()
{
  int nSteps = 1000000;
  double h = 0.0001;
  double gamma = 0.;
  int method = 0;
  std::string methodName;
  cout << "Method (0) Eul (1) Leap (2) RK4: ";
  cin >> method;
  cout << "gamma: ";
  cin >> gamma;
  cout << "h: ";
  cin >> h;
  cout << "n: ";
  cin >> nSteps;
  if(method == 0)
  {
    methodName = "Euler";
  } 
  else if(method == 1)
  {
    methodName = "Leapfrog";
  }
  else
  {
    methodName = "RungeKutta4";
  }
  std::ostringstream dataOut;
  dataOut << "./data/single/" << methodName << "_" << h << "_" << nSteps << "_" << gamma << ".txt";
  std::string dataOutDir = dataOut.str();

  std::vector<double> init(2);
  init[0] = 0.1;
  init[1] = 0.;

  
  if(method == 0)
  {
    Euler a(dataOutDir,init,h,nSteps,gamma);
    a.solve();
  } 
  else if(method == 1)
  {
    Leapfrog a(dataOutDir,init,h,nSteps,gamma);
    a.solve();
  }
  else
  {
    RungeKuttaFour a(dataOutDir,init,h,nSteps,gamma);
    a.solve();
  }
  
  return 0;
}
