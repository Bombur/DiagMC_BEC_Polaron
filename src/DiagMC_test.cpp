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


void DiagMC::test(const int & which) {
  try{
	if (G0p < (mean(Data*(G0p*taubin/stats(0,0))) - 0.0000001) || G0p > (mean(Data*(G0p*taubin/stats(0,0))) + 0.0000001)) {throw oor_Gp();} 
	if (Data((int)drnd()*taubin) == 0) {throw data_empty();} 
  }
  catch (std::exception& e){
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}
  
  
  
/*
void Ising_markov_method::initialize() {
  //E=0;
  //sigma = (irnd() < 0.5 ? VectorXi::Ones(N) : VectorXi::Constant(N, -1)) ;
  //m=sigma.dot(VectorXi::Ones(N));
  //Data.setZero((Write_its/Meas_its), 2);
  stats.setZero(0,6);
}

*/



/*
void Ising_markov_method::Etest() {
  //weight checking
  try{
	if (meanbuffer(0)<0) {throw Elesszero;}
	if (stats(4) != exp(- beta*meanbuffer(0))) {throw;}
	
	std:cout<< "Energy Test ok!" <<std::endl;
  }
  catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
	
  
}

void Ising_markov_method::magtest() {
  //weight checking
  try{
	if (abs(meanbuffer(0))>N) {throw ;}
	
	
	std:cout<< "Magnetization Test ok!" <<std::endl;
}
  catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
	
  
}

*/

  
  
  
  
  
	

