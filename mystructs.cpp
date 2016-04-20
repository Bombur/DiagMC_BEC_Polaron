#include "mystructs.h"


double expfun(const double & x) {
  if (x < -400.) {
	std::cerr<< "Exp() Underflow!" << '\t' << x << std::endl;
	return 0.;
  } else if (x > 400.){
	std::cerr<< "Exp() Overflow!"<< '\t' << x<< std::endl;
	return std::exp(100);
  }
  else{
	return std::exp(x);
  }
}