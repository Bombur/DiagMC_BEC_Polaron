//config
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
 
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
	const double taumax;
	
	//create map
	tmap(const std::vector<std::function<double(int)> > & fvec, const std::vector<int> bins, const std::vector<double> taus);
	std::vector<double> discretize(const std::function<double(int)> & f, const int & binmin , const int & binmax);
	
	
	int bin(const double & tau);  //to find right bin for tau --> DiagMC Measure
	ArrayXd print(); //for output DiagMC io
	ArrayXd norm_table();
	
	//Tests
	void print_all();
	
};


//To read in a vector from a property tree
namespace pt = boost::property_tree;

template <typename T>
std::vector<T> as_vector(const pt::ptree & ptree, const pt::ptree::key_type & key)
{
    std::vector<T> r;
    for (const auto & item : ptree.get_child(key))
        r.push_back(item.second.get_value<T>());
    return r;
}


//Functions

double mylin(const double & x, const double & m, const double & x0, const double & y0);

//Log moved to (0,0)
double mylog(const double & x, const double & a, const double & x0, const double & y0); 

//Exp moved to (0,0)
double myexp(const double & x, const double & a, const double & x0, const double & y0);

//To Read in Functions from Property Tree
std::vector<std::function<double(int)> > create_fvec(const pt::ptree & mypt, const pt::ptree::key_type & key);




#endif
