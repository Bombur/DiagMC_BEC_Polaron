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

std::vector<double> discretize(std::function<double(double) > func, std::function<double(double)> funcinv,  const double & xmin, const double & xmax, const int & bins) {
try{  
  //the function needs to be monoton increasing or decresing  
  std::vector<double> out (bins, 0);
  const double ymin = func(xmin);
  const double ymax = func(xmax);
  const double dy = (ymax -ymin)/ static_cast<double>(bins);
  if ((fabs(funcinv(ymin) - xmin) < 0.0001*xmin )|| (fabs(funcinv(ymax) - xmax)< 0.0001*xmin)) {throw fminus1();} 
  for (unsigned int i = 0; i != out.size(); i++) {
	out[i] = funcinv(ymin+ (static_cast<double>(i)*dy));	
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


class ar_size: public std::exception
{
  virtual const char* what() const throw()
  {
    return "tau map array sizes do not fit!";
  }



std::vector<double> create_map(const std::array<std::string> & funcs, const std::array<int> & bins, const std::array<double> & xlims, const std::Array<double> & para, const std::array<double> & parb) {
try{ 
  unsigned int sizetest=funcs.size();
  if ((bins.size()!= sizetest) || (xlims.size()!= sizetest+1) || (para.size()!= sizetest) || (parb.size() != sizetest)) {throw ar_size();}
  
  int total_bins=0;
  for (unsigned int i = 0; i< sizetest; i++) {
	total_bins += bins[i];
  }
  std::vector<double> map;
  map.reserve(total_bins);

  for (unsigned int = 0; i< sizetest; i++){
	if (funcs[i] == "lin") {
		std::function<double(double)> lin = std::bind(mylin, para[i], parb[i], std::placeholders::_1);
  std::function<double(double)> lininv = std::bind(mylin, 1/para[i], (-parb[i])/para[i], std::placeholders::_1);
  
		std::vector<double> linvec=discretize(lin, lininv, xlims[i], xlims[i+1], bins[i]);
		map.insert(map.end(), linvec.begin(), linvec.end());
	}

	if (funcs[i] == "log") {
		std::function<double(double)> logbla = std::bind(mylog, para[i], parb[i], std::placeholders::_1);
  std::function<double(double)> loginv = std::bind(myexp, para[i], parb[i], std::placeholders::_1);
  
		std::vector<double> logvec=discretize(logbla, loginv, xlims[i], xlims[i+1], bins[i]);
		map.insert(map.end(), logvec.begin(), logvec.end());
	}

	if (funcs[i] == "exp") {
		std::function<double(double)> expbla = std::bind(myexp, para[i], parb[i], std::placeholders::_1);
  std::function<double(double)> expinv = std::bind(mylog, para[i], parb[i], std::placeholders::_1);
  
		std::vector<double> logvec=discretize(expbla, expinv, xlims[i], xlims[i+1], bins[i]);
		map.insert(map.end(), logvec.begin(), logvec.end());
	}
  }
  
  return map;
} catch (std::exception& e) {
 	std::cerr << e.what() << std::endl;
 	exit(EXIT_FAILURE);
  }
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
