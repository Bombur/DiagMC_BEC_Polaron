#include "DiagMC.h"
 
double DiagMC::G0el(const std::vector< double > & p, const double & tfin, const double & tinit) {
  return exp(-((vsq(p)/2) - mu)*(tfin-tinit));
}

double DiagMC::Dph(const double & tfin, const double & tinit) {
  return exp((-omegap)*(tfin-tinit));
}