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

void DiagMC::insert() {
  stats(1,0) +=1;		//attempted
  stats(1,1) +=1;		//possible
  diag.propose_insert();
  
  if (drnd() < ((Prem*diag.high_weigth()*diag.P_hilo())/(Pins*diag.low_weight()*diag.P_lohi()))) {
    stats(1,3) +=1;
    diag.insert();
  }
  else {
    stats(1,2) +=1;
  }
    
}

int DiagMC::remove() {
  stats(2,0) +=1;		//attempted
  if(diag.propose_remove()!=0) {return -1;};
  stats(2,1) +=1;		//possible
  
  if (drnd() < ((Pins*diag.low_weight()*diag.P_lohi())/(Prem*diag.high_weigth()*diag.P_hilo()))) {
    stats(2,3) +=1;
    diag.remove();
  }
  else {
    stats(1,2) +=1;
  }
    
  return 0;
}


int main() {
  pt::ptree config;
  pt::read_json("DiagMC_FP.json", config);
  
  DiagMC fp(config);
  
  steady_clock::time_point time_begin = steady_clock::now();  //start time
  double nseconds;
  do {	
	//initialize();	
	for (int i=0; i<fp.Write_its; i++) {
	  for (int j=0; j<fp.Test_its; j++) {
		for (int k=0; k<fp.Meas_its; k++) {
		  if (fp.drnd() < Prem) {
			fp.remove();
		  }
		  else {
			fp.insert();			  
		  }
		  
		  //fp.change_tau();
		}
		fp.measure(j);
	  }
	  fp.test(i);
	}
	fp.status();
	
	fp.write();
	
	/*
	stats(0,4)= stats(0,3)/stats(0,1);
	stats(0,5)= stats(0,3)/stats(0,0);
		
	//Energy mean
	mean(Data.col(0));
	bin_ana(Data.col(0));
	printmean("Energy");
	
	//Magnetization mean
	mean(Data.col(1));
	bin_ana(Data.col(1));
	printmean("Magnetization");
	
	printstats();
	
	//Writing to file
	Datatofile();
	Stattofile();
	  
	cycle+=1;
	*/

	
	//Time Check
	steady_clock::time_point time_end = steady_clock::now();
	steady_clock::duration time_span = time_end-time_begin;
	nseconds = double(time_span.count()) * steady_clock::period::num / steady_clock::period::den;
  } while (nseconds < fp.RunTime);
  
  return 0;
}



	
