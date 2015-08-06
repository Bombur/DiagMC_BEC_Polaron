#include "DiagMC.h"

void DiagMC::measure(const int & whichstep) {
  Data((int)(tau/taumax*taubin), 0) +=1;
  if (get_order() < 3) {Data((int)(tau/taumax*taubin), get_order()+1) += 1;}
   
} 

