#include "ODESolution.h"
#include <fstream>
using namespace std;

int main()
{
  ofstream chaosOut;
  chaosOut.open("chaos.txt");
  std::vector<double> init(4);
  init[0] = 0.5;
  init[1] = 0.;
  init[2] = 0.1;
  init[3] = 0.1;

  for(int i = 0; i < 100; i++)
  {
      init[0] += 0.05;
      init[1] = 0;
    for(int j = 0; j< 100; j++)
    {
      init[1] += -0.05;

      RungeKuttaFour a("",init,0.01,1200,1,1);
      a.solve();  
      chaosOut << init[0] << "\t" << init[1] << "\t" << a.getTheta(1199) << endl;

    }
  }
  return 0;
 
}
