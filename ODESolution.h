#include <cmath>
#include <fstream>
#include "matrix.h"
#include <cassert>
class ODEMatrixSolution
{ 
public:
ODEMatrixSolution(std::string dataOutDir, std::vector<double> init,const double h, const int nSteps):
                  m_h(h),m_nSteps(nSteps),m_dataOutDir(dataOutDir)
{
  m_vectorHolder = new std::vector<double> [nSteps];                      //creates the dynamic array to store vectors in
  m_vectorHolder[0] = init;
  m_thetaInit = init[0];
}
~ODEMatrixSolution()
{
 if(m_vectorHolder != NULL) 
 {
   delete [] m_vectorHolder;                                              //deletes the vector array to prevent segmentation faults
   m_vectorHolder= NULL;
 }
}
void deleteArray()                                                        //useful for freeing up system memory if solving multiple objects
{
  delete [] m_vectorHolder;
  m_vectorHolder =NULL;  
}
/*
 * void solve() : the general procedure for all ODE methods
 */
void solve()
{
  std::ofstream saveFile;                      
  saveFile.open(m_dataOutDir.c_str());
  for(int i = 0; i < m_nSteps-1;i++)
  {
    findNext(i);                                                            //uses virtual function, overloaded in derived class, to store next step vector in the array
    saveFile << i*m_h << "\t";                                              //prints t' to the first collumn in data file
    saveFile << totalMechEnergy(i) << "\t";                                 //prints E to the second collumn in data file
    print(m_vectorHolder[i],saveFile);                                      //prints the vector components to he next two collumns of data file 
    if(!m_isDoublePendulum && i != 0) saveFile << "\t" << fabs((analyticSolution(i*m_h)-m_vectorHolder[i][0])/(analyticSolution((i-1)*m_h)-m_vectorHolder[i-1][0])) << std::endl;    //calculates analytic solution for single pendulum
    else saveFile << std::endl;
    if (((i+1) % 10) == 0) {                                                //simplistic counter to give solution progress
      double frac = static_cast<double> (i)/m_nSteps*100.0;
      std::cout << " Progress: " << frac << " %      \xd";
    }
  }
  saveFile.close();
}
double getTheta(int i)
{
  std::vector<double> temp = m_vectorHolder[i];
  return temp[0];
} 
protected:
  std::vector<double> * m_vectorHolder; //protected dynamic array for storing the vector steps in
  const double m_h; 
  virtual void findNext(int i) = 0;     // virtual functions to be overloaded  in derived classes
  virtual double totalMechEnergy(int i) = 0;
  virtual double analyticSolution(double t) = 0;
  const int m_nSteps;
  std::string m_dataOutDir;
  double m_thetaInit;
  bool m_isDoublePendulum;              //check for whether method is double or single

/*
 * std::vector<double> add(std::vector<double>, std::vector<double>) : adds vector a and b return a third vector c
 */

std::vector<double> add(std::vector<double> a, std::vector<double> b)
{
  std::vector<double> c(a.size());
  for(int i=0; i < a.size(); i++)
  {
    c[i] = a[i]+b[i];
  } 
  return c;
}

/*
 * std::vector<double> scale(std::vector<double>, double) : scales vector a by b
 */

std::vector<double> scale(std::vector<double> a, double b)
{
  std::vector<double> c(a.size());
  for(int i=0; i < a.size(); i++)
  {
    c[i] = a[i]*b;
  } 
  return c;
}

/*
 * void print(std::vector<double>, std::ofstream) : prints the vector components u and v to data file with tabs for easy reading by gnuplot
 */
void print(std::vector<double> _vector, std::ofstream &file)
{
	for(int i = 0; i < _vector.size(); i++)
	{
  file << _vector[i] << "\t";
  }
}

};

class Pendulum : public ODEMatrixSolution
{
public:

Pendulum(std::string dataOutDir, std::vector<double> init,const double h, const int nSteps,double gamma):
         ODEMatrixSolution(dataOutDir,init,h,nSteps),m_gamma(gamma),m_m(1),m_l(1),m_g(9.81),m_G(0),m_R(0),m_paraMatrix(2,2) //constructor for single pendulum
{
  if(init.size() !=2)
  {
    std::cout << "You entered a vector of the wrong size" << std::endl;
    assert(0);
  }
  m_paraMatrix.set_element( 0, 0, 0.                                  );                      //initialises the parameter matrix
  m_paraMatrix.set_element( 0, 1, 1                                   );
  m_paraMatrix.set_element( 1, 0, -1                                  );
  m_paraMatrix.set_element( 1, 1, -(m_gamma/(m_m*std::sqrt(m_l*m_g))) );
  m_isDoublePendulum = 0;
}

Pendulum(std::string dataOutDir, std::vector<double> init,const double h, const int nSteps,double G,double R):
         ODEMatrixSolution(dataOutDir,init,h,nSteps),m_gamma(0),m_m(1),m_l(1),m_g(9.81),m_G(G),m_R(R),m_paraMatrix(4,4) //constructor for double pendulum
{
  if(init.size() !=4)                                                                         //sanity check
  {
    std::cout << "You entered a vector of the wrong size" << std::endl;
    assert(0);
  }
  m_paraMatrix.set_element( 0, 2, 1             );                                             //initialises the parameter matrix
  m_paraMatrix.set_element( 1, 3, 1             );
  m_paraMatrix.set_element( 2, 0, -(R+1)        );
  m_paraMatrix.set_element( 2, 1, R             );
  m_paraMatrix.set_element( 2, 2, -G            );
  m_paraMatrix.set_element( 3, 0, R+1           );
  m_paraMatrix.set_element( 3, 1, -(R+1)        );
  m_paraMatrix.set_element( 3, 2, G*(1.-(1./R)) );
  m_paraMatrix.set_element( 3, 3, -(G/R)        );
  m_isDoublePendulum = 1;
}


protected:
Matrix m_paraMatrix;                                                                            //parameter matrix for the system of ODE's
private:
double m_gamma;
double m_m;
double m_l;
double m_g;
double m_G;
double m_R;

/*
 * double totalMechEnergy(int i) : checks whether system is single or double pendulum and then calculates the total mechanical energy accordingly
 */
double totalMechEnergy(int i)
{
std::vector<double> current = m_vectorHolder[i];
if(m_isDoublePendulum)
{
double theta = current[0];
double phi = current[1];
double omega = current[2];
double nu = current[3];
double kinetic = 0.5*m_l*m_g*m_m*m_R*(((1+(1./m_R))*(omega*omega))+(nu*nu)+(2.*omega*nu*cos(theta-phi)));
double potential = m_m*m_g*m_l*m_R*((1+(1./m_R))*(2-cos(theta))-cos(phi));
return kinetic + potential;
}
else return (0.5*m_m*m_l*m_g*current[1]*current[1])+(m_g*m_m*m_l*(1-cos(current[0]))); // only works for the case that theta(t=0) = 0.1;
}

/*
 * double analyticSolution(double time) : checks whether parameter imply crital, under, or over damping and then calculates the analytical value at time
 * accordingly
 */
double analyticSolution(double time) // need to consider 3 cases
{
double alpha = m_gamma*(1./m_m)*sqrt(1./(m_l*m_g));
double test = (alpha*alpha) - 4;
if( test > 0) //over damping
{
  double m1 = 0.5*((-1)*alpha - sqrt(test));
  double m2 = 0.5*((-1)*alpha + sqrt(test));
  double theta = (m_thetaInit/(m1-m2))*(m1*exp(m2*time)-m2*exp(m1*time));
  return theta;
}
else if( test == 0) //crtical damping
{
  double m1 = (-0.5)*alpha;
  double theta = m_thetaInit*(1+(0.5*alpha*time))*exp((-0.5)*alpha*time);
  return theta;
}
else if( test < 0) //under damping
{
  double theta = (((2)*m_thetaInit)/sqrt((-1)*test))*exp((-0.5)*alpha*time)*sin(acos((-0.5)*alpha)+(0.5*sqrt((-1)*test)*time));
  return theta;
}
return 0;
}
};
/*
 * RungeKuttaFour class : inherits from singlePedulum class
 */
class RungeKuttaFour : public Pendulum
{
public:
RungeKuttaFour(std::string dataOutDir,std::vector<double> init,const double h, const int nSteps,double gamma):
               Pendulum(dataOutDir,init,h,nSteps,gamma)                                //constructor for single pendulum
{}
RungeKuttaFour(std::string dataOutDir,std::vector<double> init,const double h, const int nSteps,double G,double R):
               Pendulum(dataOutDir,init,h,nSteps,G,R)                                  //constructor for double pendulum
{}
private:
void findNext(int i)
{
	std::vector<double> k1 = (m_paraMatrix*m_h)*m_vectorHolder[i];                        
  std::vector<double> k2 = (m_paraMatrix*m_h)*add(m_vectorHolder[i],scale(k1,0.5));
  std::vector<double> k3 = (m_paraMatrix*m_h)*add(m_vectorHolder[i],scale(k2,0.5));
  std::vector<double> k4 = (m_paraMatrix*m_h)*add(m_vectorHolder[i],k3);

  m_vectorHolder[i+1] = add(m_vectorHolder[i],scale(add(add(k1,scale(k2,2.)),add(scale(k3,2.),k4)), 1./6.));
}
};
/*
 * Euler class : inherits from singlePedulum class
 */
class Euler : public Pendulum
{
public:
Euler(std::string dataOutDir,std::vector<double> init,const double h, const int nSteps,double gamma):
      Pendulum(dataOutDir,init,h,nSteps,gamma),m_updateMatrix(m_paraMatrix*m_h)
{
  Matrix unit(m_paraMatrix.m_rows,m_paraMatrix.m_cols);
  unit.set_unit();
  m_updateMatrix += unit;
}

protected:
Matrix m_updateMatrix;
private:
void findNext(int i)
{
  m_vectorHolder[i+1] = m_updateMatrix*m_vectorHolder[i];
}
};
/*
 * Leapfrog class : inherits from Euler class to make use of m_updateMatrix
 * findNext is overloaded so that Leapfrog calculates according to the Leapfrog algorithm after the initial euler step
 */
class Leapfrog : public Euler
{
public:
Leapfrog(std::string dataOutDir,std::vector<double> init,const double h, const int nSteps,double gamma):Euler(dataOutDir,init,h,nSteps,gamma)
{
  m_vectorHolder[1] = m_updateMatrix*m_vectorHolder[0];  //initial euler step
  m_paraMatrix.scale(m_h*2.);                            //scales the parameter matrix as required by the leapfrog algorithm
}
private:
void findNext(int i)
{
  if(i != 0) m_vectorHolder[i+1] = add(m_vectorHolder[i-1],m_paraMatrix*m_vectorHolder[i]); //Cannot calculate n = 1 value (done by euler), so skips it
}
};
