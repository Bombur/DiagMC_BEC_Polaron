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
#include "myexceptions.h"

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
	
	//proposing data for insert
	int arc;
	std::vector<double> tau1;
	std::vector<double> tau2;
	std::vector<double> q;
	
  
  public:
	Diagram(double p, double tau, double my, double omega, double alph, int ord, std::function<double()> rnd);
	
	int get_order() {return order;};
	double get_tinit(int arc) {return times[arc][0]};
	double get_tfin(int arc) {return times[arc+1][0]};
	std::vector<double> get_q(int arc) {return phprop[arc];}
	std::vector<double> get_p(int arc) {return elprop[arc];} 
	
	void propose();
	
	double high_weigth();
	double low_weight();
	double P_lohi();
	double P_hilo();
	
	void insert();
	void remove();
	
};



   
#endif