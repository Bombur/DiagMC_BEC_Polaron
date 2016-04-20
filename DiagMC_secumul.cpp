#include "DiagMC.h"
#ifdef SECUMUL

void DiagMC::ord_step() {
  minord = maxord; 
  maxord += ordstsz;
  nnorms.push_back(0);
  nends.push_back(0);
  ordstep += 1;
  
  SEib= ArrayXd::Zero(taumap.taubin);
  Norms = ArrayXd::Zero(taumap.taubin);
  Ends = ArrayXd::Zero(taumap.taubin);
   
  Epol = ArrayXXd::Zero(taumap.taubin, ws.size());
  SE = ArrayXXd::Zero(taumap.taubin, 3);
  G0SEiw = ArrayXXcd::Zero(wbin, 2);
  
  updatestat = ArrayXXd::Zero(updatestat.rows(),updatestat.cols());
  orderstat = ArrayXd::Zero(orderstat.size());
  qstat = ArrayXXd::Zero(qstat.rows(),2);
  tstat = ArrayXXd::Zero(taumap.taubin,2);

  diag.capacity_check();
}

double DiagMC::normcalc() {
  nnorms[ordstep] = Norms.sum();
  return nnorms[ordstep]; 
}

double DiagMC::endcalc() {
  nends[ordstep] = Ends.sum();
  return nends[ordstep];
}

void DiagMC::set_av_nei(const double & av_normi, const double & av_endi, const int & ordit){
  assert(ordit == ordstep);
  nnorms[ordstep] = av_normi;
  nends[ordstep] = av_endi;
}

std::array<double, 2> DiagMC::get_minmax() {
  std::array<double, 2> tmp;
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
  pref /= nnorms[0];
  pref *= fwzero; 
  pref *= fwone;
  return pref;
}
  

ArrayXd DiagMC::get_SEib(){
  return SEib * (taumap.norm_table()*pref_calc()) ;
}

ArrayXd DiagMC::get_NormDiag(){
  return Norms * (taumap.norm_table()*pref_calc());
}

ArrayXd DiagMC::get_EndDiag(){
  return Ends * (taumap.norm_table()*pref_calc());
}
































#endif
