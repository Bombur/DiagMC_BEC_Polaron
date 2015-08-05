#include <exception>
#include <stdexcept>
#include "DiagMC.h"
  
//tests
class oor_Gp: public std::exception
{
  virtual const char* what() const throw()
  {
    return "G(p) is out of range!";
  }
};

class oor_tau: public std::exception
{
  virtual const char* what() const throw()
  {
    return "tau is out of range!";
  }
};

class data_empty: public std::exception
{
  virtual const char* what() const throw()
  {
    return "G(p, tau_i) is empty!";
  }
};

class greenerr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Error in the Data of the Green's function!";
  }
};

class weight_check: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Weight Check failed!";
  }
};

class staterr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Statistics do not match!";
  }
};

class ins_rem: public std::exception
{
  virtual const char* what() const throw()
  {
    return "More Removes than Inserts!";
  }
};


void DiagMC::test() {
  try{
	//Diagram
	diag.test();
	
	//tau
	if (tau > taumax || tau < 0) {throw oor_tau();}
	
	//Weight Check and Green's Function
	if (Data((int)(drnd()*taubin), 0) == 0) {throw data_empty();}
	int CG0p = 0;
	for (int i=0; i< taubin; i++) {
	  if (Data(i, 0)< (Data(i, 1) + Data(i, 2) + Data(i, 3))) {throw  greenerr();}
	  CG0p += Data(i, 1);	  
	}
	for (int i=0; i< taubin; i++) {
	  double G0deltatau = (exp(-(E*i*taumax/taubin))/E) *(1 - exp(-E*(taumax/taubin)));
	  //std::cout << (G0deltatau/G0p * CG0p)<<'\t' << double(Data(i,1)) << '\t' << fabs((G0deltatau/G0p * CG0p) - double(Data(i,1))) << std::endl;
	  if (((G0deltatau/G0p) * CG0p) > 100 && fabs((G0deltatau/G0p * CG0p) - double(Data(i,1))) > 0.5*(G0deltatau/G0p) * CG0p )  {throw weight_check();}
	  if (((G0deltatau/G0p) * CG0p) > 1000 && fabs((G0deltatau/G0p * CG0p) - double(Data(i,1))) > 0.16*(G0deltatau/G0p) * CG0p )  {throw weight_check();}
	  if (((G0deltatau/G0p) * CG0p) > 10000 && fabs((G0deltatau/G0p * CG0p) - double(Data(i,1))) > 0.05*(G0deltatau/G0p) * CG0p )  {throw weight_check();}
	  if (((G0deltatau/G0p) * CG0p) > 100000 && fabs((G0deltatau/G0p * CG0p) - double(Data(i,1))) > 0.016*(G0deltatau/G0p) * CG0p )  {throw weight_check();}
	}
	//Statistics
	for (int i = 0; i< 4; i++) {
	  if (stats(i,1) != stats(i,2)+stats(i, 3)) {throw staterr();}
	}
	if (stats(1,3) <stats(2,3)) {throw ins_rem();}
  }
  catch (std::exception& e){
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}
  

	

