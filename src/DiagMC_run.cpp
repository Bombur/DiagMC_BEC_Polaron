// Froehlich Polaron
// Hans Peter Guertner
#include <iostream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "DiagMC.h"


namespace pt = boost::property_tree;
using namespace std::chrono;


void DiagMC::change_tau() {
  double ntau = - log(drnd())/E;		//new tau
  //std::cout<< ntau<< '\n';
  stats(0,0) +=1;						//attempted
    
  if (ntau < taumax) {
	stats(0,1) +=1;						//possible
	stats(0,3) +=1;						//accepted
	tau = ntau;
  }
  else {
	stats(0,2) +=1;						//rejected
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
    stats(2,2) +=1;
  }
    
  return 0;
}


int main() {
  pt::ptree config;
  pt::read_json("DiagMC_FP.json", config);
  
  DiagMC fp(config);
  // check
  
  steady_clock::time_point time_begin = steady_clock::now();  //start time
  double nseconds;
  do {	
	for (int i=0; i<fp.Write_its; i++) {
	  for (int j=0; j<fp.Test_its; j++) {
		for (int k=0; k<fp.Meas_its; k++) {
		  if (fp.drnd() < fp.Prem) {
			fp.remove();
		  }
		  else {
			fp.insert();			  
		  }
		  
		  //fp.change_tau();
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
  
  return 0;
}



	
