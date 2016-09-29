//config
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

//io
#include <iostream>
#include <boost/concept_check.hpp>

//string
#include <string>

//find bottleneck
#include <array>
#include <algorithm>

//exceptions
#include <stdexcept>

#ifndef __ADAPTION_H_INCLUDED__
#define __ADAPTION_H_INCLUDED__

namespace pt = boost::property_tree;

enum bottleneck {NORM = 0, END, TIME};

class Adaption{
  private:
	int ordstsz; //current order step size
	const double nnorms; //Number in NormDiag
	const double nends; //Number in EndDiag
	const double tls; //Time needed for last step
	
	//Expected Value
	const double normmin; 
	const double endmin;
	const double texp;
	
	//Relative Deviations
	double countdev; //end or norm dev
	double normdev; 
	double enddev;
	double tdev;
  public:
	Adaption(const int & oss, const double & time, const double & normal, const double & ende, const pt::ptree & config);
	bottleneck whichbtlnk(); //which deviation is the biggest
	int ordstszadapt();
	
	void printall();

};


#endif