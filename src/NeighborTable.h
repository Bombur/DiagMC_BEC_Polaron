#include<iostream>
#include "/project/theorie/h/H.Guertner/lib/Eigen/Eigen/Dense"
#ifndef __NeighborTable_H_INCLUDED__
#define __NeighborTable_H_INCLUDED__

using namespace Eigen;
using namespace std;
 
class NeighborTable {
  protected:
	int dim;
	MatrixX4i Nbr;
  public:
	NeighborTable();
	NeighborTable(int);
	MatrixX4i get_Table() const;
	void print_Table();
};

#endif