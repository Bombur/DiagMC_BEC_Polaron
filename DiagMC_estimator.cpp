#include "DiagMC.h"

//Propagators
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
#endif

#ifdef BEC
double DiagMC::Dph(const std::array< double, 3> & q, const double & tfin, const double & tinit) {
  double omega = sqrt(vsq(q) *(1.+(vsq(q)/2.)));
  return exp((-omega)*(tfin-tinit));
}

double DiagMC::Vq2(const std::array< double, 3> & q) {
  return alpha *M_PI * (1.+ 1./relm)*(1. + 1./relm)* sqrt(vsq(q)/(vsq(q)+2.));   
}
#endif

//Estimators
void DiagMC::meas_Epol(const double & cor){
  if (diag.get_order() > 0) { // && !(diag.is_reducible())){ //&& (vsq(diag.get_p(3)-diag.get_p(1))< 1e-8)) {
	double Eptau = diag.get_tinit(2*diag.get_order()) - diag.get_tinit(1);
	Epol.row(taumap.bin(Eptau)) += ((ws*Eptau - mu*Eptau).exp() *cor);
  }
}

ArrayXXd DiagMC::get_Eptest() {
  ArrayXXd output = Epol;
  for (int i=0 ; i< Epol.cols() ; i++) {
	Epol.col(i) *= taumap.norm_table();
	//std::cout<<taumap.norm_table() <<std::endl;
  }
  //std::cout <<G0p <<'\t' << Data.col(1).sum() << std::endl;
  return output /G0p/ Data.col(1).sum();
}

ArrayXd DiagMC::get_Ep() {
  ArrayXXd output = Epol;
  for (int i=0 ; i< Epol.cols() ; i++) {
	Epol.col(i) *= taumap.norm_table();
  }
#ifndef SECUMUL
  return output.colwise().sum() / G0p / Data.col(1).sum();
#else
  return output.colwise().sum() / G0p *pref_calc();
#endif
}
