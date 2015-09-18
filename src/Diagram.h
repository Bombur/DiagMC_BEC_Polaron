//config
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
 
//io
#include <iostream>
#include <fstream>
#include <string>
 
//random
#include <random>
#include <functional>

//exceptions
#include <exception>
#include <stdexcept>

//math and container
#include "/project/theorie/h/H.Guertner/lib/Eigen/Eigen/Dense"
#define _USE_MATH_DEFINES
#include <cmath>
#include "dvector.h"

#ifndef __DIAGRAM_H_INCLUDED__
#define __DIAGRAM_H_INCLUDED__

namespace pt = boost::property_tree;
using namespace Eigen;

class Diagram {
  private:
	int order;
	std::function<double()> drnd;
	std::vector< std::vector<double> > times;
	std::vector< std::vector<double> > phprop;
	std::vector< std::vector<double> > elprop;
	int sw_pos; 					// possible vertices for swap
	  
  public:
	Diagram();
	Diagram(const double & p, const double & tau, const std::function<double()> & rnd);
	void set(const double & p, const double & tau, const std::function<double()> & rnd);
	 
	int get_order() {return order;}
	double get_tau() {return times[2*order+1][0];}
	double get_tinit(const int & arc) {return times[arc][0];}
	double get_tfin(const int & arc) {return times[arc+1][0];}
	int get_link(const int & arc) {return (int)(times[arc][1]+0.5);}
	double get_sw_pos() {return static_cast<double>(sw_pos);}
	std::vector<double> get_q(const int & arc) {return phprop[arc];}
	std::vector<double> get_p(const int & arc) {return elprop[arc];}
	
	 
	//proposing pr_
	int pr_arc;
	double pr_tauin, pr_taufin;
	std::vector<double> pr_tau1;
	std::vector<double> pr_tau2;
	std::vector<double> pr_q;
	std::vector<double> pr_p;
	

	//proposing
	void random_arc();
	int propose_insert();
	int propose_remove();
	int propose_swap();
	int propose_ct(const double & taumax, const double & ctcor);
	int propose_dq(const double & qcor);
	

	//changes
	void insert();
	void remove();
	int set_tau(double tau);
	void swap(const double &);
	void ct();
	void dq();
	
	//tests
	void test();
	void printall(); 
	
	
};



   
#endif