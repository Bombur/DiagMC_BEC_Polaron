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
  
  SEib= ArrayXi::Zero(taubin);
  Norms = ArrayXi::Zero(taubin);
  Ends = ArrayXi::Zero(taubin);
}


int DiagMC::normcalc() {
  if(ordstep == 0) { // case where everything is still saved to data 
	nnorms.at(0) = Data.col(1).sum();
  } else {
	nnorms[ordstep] = Norms.sum();
  }  
  return nnorms[ordstep]; 
}

int DiagMC::endcalc() {
  nends[ordstep] = Ends.sum();
  return nends[ordstep];
}

std::vector<int> DiagMC::get_minmax() {
  std::vector<int> tmp(2);
  tmp[0] = minord;
  tmp[1] = maxord;
  return tmp;
}

double DiagMC::pref_calc(){
  double pref = 1.;
  for (int i = 0; i < ordstep; i++) {
	pref *= static_cast<double>(nends[i]);
	pref /= static_cast<double>(nnorms[i+1]);
  }
  return pref * static_cast<double>(taubin)/taumax/static_cast<double>(nnorms[0]);
}
  

ArrayXd DiagMC::get_SEib(){
  return SEib.cast<double>() * pref_calc();
}

ArrayXd DiagMC::get_NormDiag(){
  return Norms.cast<double>() * pref_calc();
}

ArrayXd DiagMC::get_EndDiag(){
  return Ends.cast<double>() * pref_calc();
}
































#endif