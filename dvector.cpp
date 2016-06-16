#include "dvector.h"

void eigen_push_back(ArrayXXd & array, const VectorXd & values) {
	assert(values.size() == array.cols());
	int row = array.rows();
	array.conservativeResize(row+1, NoChange);
	array.row(row) = values;
}

