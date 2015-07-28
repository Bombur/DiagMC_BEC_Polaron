// Advanced Computational Physics
// Exercise 3.1
#include <iostream>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "DiagMC.h"


namespace pt = boost::property_tree;
using namespace std::chrono;


int main() {
  pt::ptree config;
  pt::read_json("Froehlich.json", config);
  
  DiagMC fp(config);
  
  steady_clock::time_point time_begin = steady_clock::now();  //start time
  int cycle=0;
  double nseconds;
  do {	
	//initialize();	
	for (int i=0; i<fp.Write_its); i++) {
	  for (int j=0; j<fp.Test_its; j++) {
		for (int k=0; k<fp.Meas_its; k++) {
		  if (arcs==0) {
			fp.change_tau();
		  }
		}
		fp.measure(j);
	  }
	  fp.test();
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
  } while (nseconds < RunTime);
  
  Datafile.close();
  Statsfile.close();

 
  return 0;
}



	
