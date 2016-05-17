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
  SE = ArrayXXd::Zero(taumap.taubin, SE.cols());
  G0SEiw = ArrayXXcd::Zero(wbin, G0SEiw.cols());
  
  updatestat = ArrayXXd::Zero(updatestat.rows(),updatestat.cols());
  orderstat = ArrayXd::Zero(orderstat.size());
  qstat = ArrayXXd::Zero(qstat.rows(),qstat.cols());
  tstat = ArrayXXd::Zero(taumap.taubin,tstat.cols());

  diag.capacity_check();
  
  if (desi_ordrat>0) {
	fw = 1.;
	fw_counts.fill(0.);
	fw_last[0] = 0;
	fw_last[1] = 0;
  }
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


//Each time write() is called this function checks the ratio between the counts in minimum and maximum order and changes the fake weight according to the desired ratio 
int DiagMC::fw_adapt(){
  assert(((orderstat(minord) - fw_last[0]) != 0) || ((orderstat(maxord) - fw_last[1]) != 0));
  fw_counts[0] = static_cast<double>(orderstat[minord]);
  fw_counts[1] += static_cast<double>(orderstat[maxord]-fw_last[1]); ///pow(fw, static_cast<double>(maxord-minord));
  if ((orderstat[minord] - fw_last[0])  == 0) {fw /= fw_max;}
  else if ((orderstat[maxord] - fw_last[1]) == 0) {fw *= fw_max;}
  else {
	double max_min_rat = fw_counts[1]/fw_counts[0]; //static_cast<double>(orderstat[maxord] - fw_last[1])/static_cast<double>(orderstat[minord] - fw_last[0]);
	double fw_new =  pow(1./max_min_rat*desi_ordrat, 1./static_cast<double>(maxord-minord));
	if ((fw_new/fw) > fw_max) {fw *= fw_max;} 
	else if ((fw/fw_new) > fw_max) {fw /= fw_max;} 
	else {fw *= fw_new;}
	std::cout << "Fake Weight Control!"  << '\n';
	std::cout << fw << '\t' << orderstat[maxord]/orderstat[minord]<< '\t' << orderstat[minord] << '\t' <<orderstat[maxord] <<std::endl;
	std::cout << fw_counts << '\t' << fw_last <<std::endl;
  }  
  
  fw_last[0] = orderstat[minord];
  fw_last[1] = orderstat[maxord];
  return 1;
}





























#endif
