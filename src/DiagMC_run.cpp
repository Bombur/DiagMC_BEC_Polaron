// Froehlich Polaron
// Hans Peter Guertner
#include <iostream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "DiagMC.h"


namespace pt = boost::property_tree;
using namespace std::chrono;

class oor_Probs: public std::exception {
  virtual const char* what() const throw()
  {
    return "No possible action chosen";
  }
};

void DiagMC::change_tau() {
  double ntau = - log(drnd())/E;		//new tau
  //std::cout<< ntau<< '\n';
  stats(0,0) +=1;						//attempted
    
  if (ntau < taumax) {
	if (diag.set_tau(ntau) == 0) {
	  stats(0,1) +=1;						//possible
	  stats(0,3) +=1;						//accepted
	  tau = ntau;
	}
  }
}

int DiagMC::insert() {
  stats(1,0) +=1;		//attempted
  if(diag.propose_insert()!=0) {return -1;};
  stats(1,1) +=1;		//possible
  if (drnd() < ((Prem*diag.high_weight()*diag.P_hilo())/(Pins*diag.low_weight()*diag.P_lohi()))) {
    stats(1,3) +=1;		//accepted
    diag.insert();
  }
  else {
    stats(1,2) +=1;		//rejected
  }  
  return 0;    
}

int DiagMC::remove() {
  stats(2,0) += 1;		//attempted
  if(diag.propose_remove()!=0) {return -1;};
  stats(2,1) +=1;		//possible
  //std::cout<< Pins << '\t' << diag.low_weight() << '\t' << diag.P_lohi() << '\t' << Prem <<'\t' << diag.high_weight() << '\t' << diag.P_hilo()<<std::endl;
  if (drnd() < ((Pins*diag.low_weight()*diag.P_lohi())/(Prem*diag.high_weight()*diag.P_hilo()))) {
    stats(2,3) +=1;		//accepted
    diag.remove();
  }
  else {
    stats(2,2) +=1;		//rejected
  }
    
  return 0;
}

int DiagMC::swap() {
  stats(3,0) += 1; //attempted
  stats(4,0) += 1;
  stats(5,0) += 1;
  stats(6,0) += 1;
  int which = diag.propose_swap();
  if (which == -1) {return -1;} //not possible
  
  // both vertices are open or closed
  if (which == 0) {
	stats(4,1) +=1; //possible
	if (drnd() < (diag.G0el(diag.get_prp())/diag.G0el(diag.get_p()))) {
	  stats(4,3) +=1; //accepted
	  diag.swap();
	}
	else {
	  stats(4,2)+= 1; 	// rejected
	}
  }
  
  // first vertex is start point, second one is end point 
  if (which == -2) {
	stats(5,1) +=1; //possible
	if (drnd() < (diag.G0el(diag.get_prp())/diag.G0el(diag.get_p())) * diag.Dph(double(which))) {
	  stats(5,3) +=1; //accepted
	  diag.swap();
	}
	else {
	  stats(5,2)+= 1; 	// rejected
	}
  }
  
   // first vertex is end point, second one is start point 
  if (which == 2) {
	stats(6,1) +=1; //possible
	if (drnd() < (diag.G0el(diag.get_prp())/diag.G0el(diag.get_p())) * diag.Dph(double(which))) {
	  stats(6,3) +=1; //accepted
	  diag.swap();
	}
	else {
	  stats(6,2)+= 1; 	// rejected
	}
  } 
  
  for (int i =1 ; i<4 ; i++){
	stats(3,i) = stats(4,i) + stats(5,i) + stats(6,i);
  }
  
  return 0;
}


int main() {
  try{
	
	pt::ptree config;
	pt::read_json("DiagMC_FP.json", config);
  
	DiagMC fp(config);
  
	steady_clock::time_point time_begin = steady_clock::now();  //start time
	double nseconds;
	do {	
	  for (int i=0; i<fp.Write_its; i++) {
		for (int j=0; j<fp.Test_its; j++) {
		  for (int k=0; k<fp.Meas_its; k++) {
			double action = fp.drnd();
			if (action < fp.Prem) {
			  fp.remove();
			}
			else if ((action-fp.Prem)<fp.Pins) {
			  fp.insert();
			}
			else if ((action-fp.Prem-fp.Pins) < fp.Pct) {
			  fp.change_tau();
			}
			else if ((action- fp.Prem- fp.Pins - fp.Pct) < fp.Psw) {
			  fp.swap();
			}
			else {
			  throw oor_Probs();
			}
		  }
		  fp.measure(j);
		}
		fp.status();
		fp.test();
	  }
	  fp.status();
	
	  fp.write();
	
	  steady_clock::time_point time_end = steady_clock::now();
	  steady_clock::duration time_span = time_end-time_begin;
	  nseconds = double(time_span.count()) * steady_clock::period::num / steady_clock::period::den;
	} while (nseconds < fp.RunTime);
	
	fp.Stattofile();
	
	return 0;
  } catch (std::exception& e){
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}



	
