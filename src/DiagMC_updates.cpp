#include "DiagMC.h"

void DiagMC::measure(const int & whichstep) {
  Data((int)(tau/taumax*taubin), 1) +=1;
  if (get_order() < 2) {Data((int)(tau/taumax*taubin), get_order()+2) += 1;}
   
}

