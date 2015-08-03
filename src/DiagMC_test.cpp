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


void DiagMC::test(const int & which) {
  try{
	//Diagram
	diag.test();
		
	//Weight Check and Green's Function
	if (Data((int)drnd()*taubin, 1) == 0) {throw data_empty();}
	int CG0p = 0;
	for (int i=0; i< taubin; i++) {
	  if (Data(i, 1)< (Data(i, 2) + Data(i, 3) + Data(i, 4))) {throw  greenerr();}
	  CG0p += Data(i, 2);	  
	}
	for (int i=0; i< taubin; i++) {
	  double G0deltatau = (exp(-(E*i*taumax/taubin))/E) *(1 - exp(-E*(taumax/taubin)));
	  if ((G0deltatau/G0p * CG0p) < (Data(i,2)-1) || (G0deltatau/G0p * CG0p) > (Data(i,2)+1)) {throw weight_check();}
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
  

	

