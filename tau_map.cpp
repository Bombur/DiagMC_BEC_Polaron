//exceptions
#include <exception>
#include <stdexcept>

//io
#include <iostream>
#include <string>

//math and container
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

//functions
#include <functional>



class fminus1: public std::exception
{
  virtual const char* what() const throw()
  {
    return "func and its inverse do not fit!";
  }
};

template<typename T> 
std::vector<T> discretize(std::function<T(T) > func, std::function<T(T)> funcinv,  const T & xmin, const T & xmax, const int & bins) {
try{  
  //the function needs to be monoton increasing or decresing  
  std::vector<T> out (bins, 0);
  const T ymin = func(xmin);
  const T ymax = func(xmax);
  const T dy = (ymax -ymin)/ static_cast<T>(bins);
  if ((static_cast<double>(fabs(funcinv(ymin) - xmin)) < 0.0001*static_cast<double>(xmin) )|| (fabs(funcinv(ymax) - xmax)< 0.0001*static_cast<double>(xmin))) {throw fminus1();} 
  for (unsigned int i = 0; i != out.size(); i++) {
	out[i] = funcinv(ymin+ (static_cast<T>(i)*dy));	
  }
  return out;
  
  
  } catch (std::exception& e) {
 	std::cerr << e.what() << std::endl;
 	exit(EXIT_FAILURE);
  }
   
}


double mylog(const double & a, const double & b, const double & x) {
  return log(a * x) / b;
}

double myexp(const double & a, const double & b, const double & x) {
  return a * exp(b * x);
}


int main() {
  std::function<double(double)> fl = std::bind(mylog, 2, 3, std::placeholders::_1);
  std::function<double(double)> fe = std::bind(myexp, 2, 3, std::placeholders::_1);
  
  std::vector<double> tmp1 = discretize<double>(fl, fe, 1., 2., 5);
  
  for (unsigned int i = 0; i != tmp1.size(); i++) {
	std::cout<< tmp1[i] << std::endl;
  }
  
  return 0;
}