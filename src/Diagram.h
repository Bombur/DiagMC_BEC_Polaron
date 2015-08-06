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

#ifndef __DIAGRAM_H_INCLUDED__
#define __DIAGRAM_H_INCLUDED__

namespace pt = boost::property_tree;
using namespace Eigen;

class Diagram {
  private:
	double mu;
	double wp;
	double alpha;
	int order;
	std::function<double()> drnd;
	std::vector< std::vector<double> > times;
	std::vector< std::vector<double> > phprop;
	std::vector< std::vector<double> > elprop;
	
	//proposing pr_
	int pr_arc;
	double pr_tauin, pr_taufin;
	std::vector<double> pr_tau1;
	std::vector<double> pr_tau2;
	std::vector<double> pr_q;
	std::vector<double> pr_p;
	
  
  public:
	Diagram();
	Diagram(const double & p, const double & tau, const double & my, const double & omega, const double & alph, const std::function<double()> & rnd);
	void set(const double & p, const double & tau, const double & my, const double & omega, const double & alph, const std::function<double()> & rnd);
	
	int get_order() {return order;}
	double get_tinit(const int & arc) {return times[arc][0];}
	double get_tfin(const int & arc) {return times[arc+1][0];}
	std::vector<double> get_q(const int & arc) {return phprop[arc];}
	std::vector<double> get_p(const int & arc) {return elprop[arc];}
	
	//just for swap
	std::vector<double> get_p() {return elprop[pr_arc];}
	std::vector<double> get_prp() {return pr_p;} 
	
	//proposing
	void random_arc();
	int propose_insert();
	int propose_remove();
	int propose_swap();
	int propose_ct();
	
	//weights
	double high_weight();
	double low_weight();
	double P_lohi();
	double P_hilo();
	double G0el(const std::vector<double> &);
	double Dph(const double & factor);
	
	//changes
	void insert();
	void remove();
	int set_tau(double tau);
	void swap();
	double ct();
	
	//tests
	void test();
	void printall();
	
	
};



   
#endif