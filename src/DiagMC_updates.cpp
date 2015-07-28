#include "DiagMC.h"

void DiagMC::measure(const int & whichstep) {
  Data(whichstep, 0) = tau;
  if (arcs==0) {
	Data(whichstep, 1) = exp(-E * tau);
  }
}


/*


*/