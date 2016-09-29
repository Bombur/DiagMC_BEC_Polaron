//math and container
#include <vector>
#include <array>
#include <cmath>

#include <iostream>

#ifndef __MYSTRUCTS_H_INCLUDED__
#define __MYSTRUCTS_H_INCLUDED__



struct arch {
  std::array<double,3> mom;
  int beg;
  int link;
};

struct vertex {
  double t;
  int link;
};

double expfun(const double & x);

#endif