#include "NeighborTable.h"

NeighborTable::NeighborTable() {
  dim = 3;
  Nbr = (MatrixX4i(9,4) << 2, 4, 0, 0,
						  3, 5, 1, 0,
						  0, 6, 2, 0,
						  5, 7, 0, 1,
						  6, 8, 4, 2,
						  0, 9, 5, 3,
						  8, 0, 0, 4,
						  9, 0, 7, 5,
						  0, 0, 8, 6).finished();
}

NeighborTable::NeighborTable(int d) {
	dim = d;
	Nbr = MatrixX4i(dim*dim, 4);
	for (int i=1; i< 1+d*d; i++) {
	  Nbr(i-1, 0)= ((i % dim) == 0 ? 0 : i+1);
	  Nbr(i-1, 1)= (i > (dim*(dim-1)) ? 0 :  i+d);
	  Nbr(i-1, 2)= (((i-1) % dim) == 0 ? 0 : i-1);
	  Nbr(i-1, 3)= (i < (dim+1)) ? 0 : i-d;
	}
}

MatrixX4i NeighborTable::get_Table() const {
  return Nbr;  
}

void NeighborTable::print_Table() {
  cout << Nbr << endl;
}