#include "DiagMC.h"

void DiagMC::measure(const int & whichstep) {
  Data((int)(tau/taumax*taubin)) += 1;
}

