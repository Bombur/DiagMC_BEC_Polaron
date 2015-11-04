//math and container
#include <vector>
#include <array>

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

#endif