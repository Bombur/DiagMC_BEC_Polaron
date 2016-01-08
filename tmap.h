//io
#include <iostream>

//string
#include <string>

//exceptions
#include <stdexcept>
#include <assert.h>

//math and container
#include "/project/theorie/h/H.Guertner/lib/Eigen/Eigen/Dense"
#include <vector>
#include <map>
#define _USE_MATH_DEFINES
#include <cmath>

//functions
#include <functional>

//Iterator
#include <iterator>

#ifndef __TMAP_H_INCLUDED__
#define __TMAP_H_INCLUDED__

using namespace Eigen;

enum funcs {LIN = 0, EXP, LOG};

class tmap {
  private:
	std::map<double, int> mymap;
	
  public:
	const int taubin; //How many bins in total
	
	//create map
	tmap(const std::vector<std::function<double(int)> > & fvec, const std::vector<int> bins, const std::vector<double> taus);
	std::vector<double> discretize(const std::function<double(int)> & f, const int & binmin , const int & binmax);
	
	
	int bin(const double & tau);  //to find right bin for tau --> DiagMC Measure
	ArrayXd print(); //for output DiagMC io
	
	//Tests
	void print_all();
	
};



#endif
