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
#include <assert.h>
#include "DiagramException.h"

//math and container
#include "/project/theorie/h/H.Guertner/lib/Eigen/Eigen/Dense"
#define _USE_MATH_DEFINES
#include <cmath>
#include "dvector.h"
#include <array>
#include "mystructs.h"
#include <algorithm>

#ifndef __DIAGRAM_H_INCLUDED__
#define __DIAGRAM_H_INCLUDED__

namespace pt = boost::property_tree;
using namespace Eigen;


 
class Diagram {
  private:
	std::function<double()> drnd;
	int order;
	std::vector< vertex > times;
	std::vector< std::array<double,3> > phprop; // to control the electron prop vector
	std::vector< std::array<double,3> > elprop;
	  
  public:
	Diagram();
	void set(const double & p, const double & tau, const std::function<double()> & rnd);
	 
	int get_order() {return order;}
	double get_tau() {return times[2*order+1].t;}
	double get_tinit(const int & arc) {return times[arc].t;}
	double get_tfin(const int & arc) {return times[arc+1].t;}
	int get_link(const int & arc) {return times[arc].link;}
	std::array<double,3> get_q(const int & arc) {return phprop[arc];}
	std::array<double,3> get_p(const int & arc) {return elprop[arc];}
	
	  
	//proposing pr_ 
	int pr_arc;
	double pr_tauin, pr_taufin;
	vertex pr_tau1;
	vertex pr_tau2;
	std::array<double,3> pr_q;
	std::array<double,3> pr_p;
	 

	//proposing
	void random_arc();
	int propose_insert(const double & dqins, const double & sigfac);
	int propose_remove(const double & dqins);
	int propose_swap();
	int propose_ct(const double & taumax, const double & ctcor);
	int propose_dq(const double & qcor);
	int propose_insatend(const double & taumax, const double & dtins, const double & dqins);
	int propose_rematend(const double & dtins, const double & dqins);
	

	//changes
	void insert();
	void remove();
	int set_tau(double tau);
	void swap();
	void ct();
	void dq();
	void insatend();
	void rematend();
	
	//tests
	void test(const double & qc);
	void printall(); 
	bool is_reducible();
	int arch_num(const int &);
	
	
};



   
#endif
