#include "DiagMC.h"
 
double DiagMC::G0el(const std::array< double, 3 > & p, const double & tfin, const double & tinit) {
  return exp(-((vsq(p)/2.) - mu)*(tfin-tinit));
}

#ifdef FP
double DiagMC::Dph(const std::array< double, 3> & q, const double & tfin, const double & tinit) {
  return exp((-wp)*(tfin-tinit));
}

double DiagMC::Vq2(const std::array< double, 3> & q) {
  double pref =  alpha*2.*sqrt(2.)*M_PI;
  return pref/vsq(q);  
}
#else
double DiagMC::Dph(const std::array< double, 3> & q, const double & tfin, const double & tinit) {
  double omega = sqrt(vsq(q)/2 * (1+(vsq(q)/2)));
  return exp((-omega)*(tfin-tinit));
}

double DiagMC::Vq2(const std::array< double, 3> & q) {
  double pref = M_PI * alpha / 2 *(1+(1/relm));
  return pref * sqrt(vsq(q)/(vsq(q)+2));  
}
#endif