#include "DiagMC.h"
  

/*
void Ising_markov_method::initialize() {
  //E=0;
  //sigma = (irnd() < 0.5 ? VectorXi::Ones(N) : VectorXi::Constant(N, -1)) ;
  //m=sigma.dot(VectorXi::Ones(N));
  //Data.setZero((Write_its/Meas_its), 2);
  stats.setZero(0,6);
}

*/

void DiagMC::change_tau() {
  double ntau = - ln(drnd())/E;			//new tau
  stats(0,0) +=1;						//attempted
  stats(0,1) +=1;						//possible
  
  if (drnd() < taumax) { 			
	stats(0,3) +=1;						//accepted
	tau = ntau;
  else {
	stats(0,2) +=1;						//rejected
  }
}

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

  
  
  
  
  
	

