#include "tmap.h"


double mylin(const double & x, const double & m = 5./200., const double & x0 = 0., const double & y0 = 0.) {
  return (m * (x - x0)) + y0;
}

//Log moved to (0,0)
double mylog(const double & x, const double & a, const double & x0 = 0. , const double & y0 = 0.) {
  return a * log(x + 1. - x0) + y0;
}

//Exp moved to (0,0)
double myexp(const double & x, const double & a, const double & x0 = 0. , const double & y0 = 0.) {
  return exp(a*(x - x0))- 1. + y0;
}



int main() {
using namespace std::placeholders;

std::vector<std::function<double(int)>> fvec;
fvec.push_back(std::bind(mylin, _1, 5./200., 0., 0.));
fvec.push_back(std::bind(myexp, _1, 5./200., 100., 2.5));
fvec.push_back(std::bind(mylog, _1 , 0.29 , 200., 13.6825));
fvec.push_back(std::bind(mylin, _1, 0.007, 250., 14.8227));
std::vector<int> bins = {0, 100, 200, 250, 300};
std::vector<double> taus = {0, 2.5, 13.6825, 14.8227, 15.1727};

tmap testmap(fvec, bins, taus);

testmap.print_all();

std::cout <<testmap.bin(14.9)<< "\n\n" << testmap.print()<< std::endl;

return 0;
}
