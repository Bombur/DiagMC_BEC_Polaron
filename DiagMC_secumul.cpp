#include "DiagMC.h"
#ifdef SECUMUL

int DiagMC::check_ordstsz(const double & Timeperorderstep) {
  if (Timeperorderstep > RunTime && ordstsz > 1) {
	return ordstsz - 1;
  }
  if (Timeperorderstep < RunTime) {
	return  ordstsz + 1;
  }
  return ordstsz;
} 


void DiagMC::ord_step() {
  minord = maxord; 
  maxord += ordstsz;
  nnorms.push_back(0);
  nends.push_back(0);
  ordstep += 1;
  
  SEib= ArrayXd::Zero(taubin);
  Norms = ArrayXd::Zero(taubin);
  Ends = ArrayXd::Zero(taubin);
}


double DiagMC::normcalc() {
  nnorms[ordstep] = Norms.sum();
  return nnorms[ordstep]; 
}

double DiagMC::endcalc() {
  nends[ordstep] = Ends.sum();
  return nends[ordstep];
}

std::vector<double> DiagMC::get_minmax() {
  std::vector<double> tmp(2);
  tmp[0] = static_cast<double>(minord);
  tmp[1] = static_cast<double>(maxord);
  return tmp;
}

double DiagMC::pref_calc(){
  double pref = 1.;
  for (int i = 0; i < ordstep; i++) {
	pref *= nends[i];
	pref /= nnorms[i+1];
  }
  return pref * static_cast<double>(taubin)/taumax/nnorms[0];
}
  

ArrayXd DiagMC::get_SEib(){
  return SEib * pref_calc();
}

ArrayXd DiagMC::get_NormDiag(){
  return Norms * pref_calc();
}

ArrayXd DiagMC::get_EndDiag(){
  return Ends * pref_calc();
}
































#endif