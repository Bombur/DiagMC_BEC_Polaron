#include "DiagMC.h"
#ifdef SECUMUL

void DiagMC::ord_step() {
  ordstep += 1;
  if (ordstep < ord_tab.size()){
	if ((TotMaxOrd > (maxord+ord_tab[ordstep]))||(TotMaxOrd <0)) {ordstsz = ord_tab[ordstep];}
  }	
  minord = maxord; 
  maxord += ordstsz;
  nnorms.push_back(0);
  nends.push_back(0);
  
  SEib= ArrayXd::Zero(taumap.taubin);
  Norms = ArrayXd::Zero(taumap.taubin);
  Ends = ArrayXd::Zero(taumap.taubin);
  
  Epol.assign(ws.size(), 0.);
  SE = ArrayXXd::Zero(taumap.taubin, SE.cols());
  SEacc = create_empty_SE_acc("step"+std::to_string(ordstep)); 
  G0SEiw = ArrayXXcd::Zero(wbin, G0SEiw.cols());
  
  counts = create_empty_count_acc(ordstep);
  Epbin = create_empty_Ep_acc("step" + std::to_string(ordstep));
  last_g0_count = 0.;
  last_measured_Epol = Epol;
  Ep_intv = create_empty_Ep_acc("step" + std::to_string(ordstep));
  ordesti = create_empty_ordesti(ordstep);
  
  updatestat = ArrayXXd::Zero(updatestat.rows(),updatestat.cols());
  if (maxord > orderstat.size()) {orderstat.resize(2*orderstat.size());}
  orderstat = ArrayXd::Zero(orderstat.size());
  qstat = ArrayXXd::Zero(qstat.rows(),qstat.cols());
  tstat = ArrayXXd::Zero(taumap.taubin,tstat.cols());

  diag.capacity_check();
  
  if (which_fw_ad==1) {
    fw_last.assign(2, 0);
    fw_counts.assign(2, 0.);
  } else if (which_fw_ad ==2){
    fw_last.assign(ordstsz+1, 0);
	fw_counts.assign(ordstsz+1, 0.);
  }
  fw_vec.assign(ordstsz, 1.);  
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
  fw_counts[0] = orderstat(minord);
  fw_counts[1] += (orderstat(maxord)-fw_last[1])*fw_for_meas(maxord-minord, true);

  std::cout << "-------------------Fake Weight Adapt!"  << '\n';
  
  if ((orderstat(minord) - fw_last[0])  < 10) {fw /= fw_max;}
  else if ((orderstat(maxord) - fw_last[1]) < 10) {fw *= fw_max;}
  else {
	double max_min_rat = fw_counts[1]/fw_counts[0]; 
	double fw_new =  pow(1./max_min_rat, 1./static_cast<double>(maxord-minord));
	std::cout << "Calculated Min_Max_Ratio and fw_new: "<< 1./max_min_rat << '\t' << fw_new << std::endl;
	if ((fw_new/fw) > fw_max) {fw *= fw_max;} 
	else if ((fw/fw_new) > fw_max) {fw /= fw_max;} 
	else {fw = fw_new;}
  }  
  std::cout << fw << '\t' << orderstat(minord)-fw_last[0] << '\t' << orderstat(maxord)-fw_last[1] << '\t' << fw_counts <<std::endl;
  
  fw_last[0] = orderstat(minord);
  fw_last[1] = orderstat(maxord);
  return 1;
}


int DiagMC::fw_vec_adapt(){
	std::cout << "------------------Fake Weight Vec Adapt!"  << '\n';
  std::cout << "Old fw_vec:\n" << fw_vec << std::endl;
  std::cout << "Desired Ratio:\n" << pow(desi_rat, 1./(maxord-minord)) << std::endl;
  std::cout << "Calculated Min_Max_Ratio:\n(";
  fw_counts[0] = orderstat(minord);
  for (int orddiff = 1; orddiff <=(maxord-minord); orddiff++){
    fw_counts[orddiff] += (orderstat(minord+orddiff)-fw_last[orddiff])*fw_for_meas(orddiff, true);
	
	if ((fw_counts[orddiff-1]> 0.) && (fw_counts[orddiff]>0.)){	
	  double fw_new =  fw_counts[orddiff-1]/fw_counts[orddiff]*pow(desi_rat, 1./(maxord-minord));
	  std::cout << fw_new << ", ";
	  if ((fw_new/fw_vec[orddiff-1]) > fw_max) {fw_vec[orddiff-1] *= fw_max;} 
	  else if ((fw_vec[orddiff-1]/fw_new) > fw_max) {fw_vec[orddiff-1] /= fw_max;} 
	  else {fw_vec[orddiff-1] = fw_new;}
	} else{
	  std::cout <<  "NaN, ";
	}
  }
  std::cout << ")" << std::endl;
  std::cout << "New fw_vec:\n" << fw_vec << std::endl;
  
  std::cout << "Orderstat:\n(";
  for (int orddiff = 0; orddiff <=(maxord-minord); orddiff++){
	  std::cout << orderstat(minord+orddiff)-fw_last[orddiff]<< ", ";
		fw_last[orddiff] = orderstat(minord+orddiff);
  } 
  std::cout << ")" << std::endl;
	std::cout << "------------------"  << std::endl;
  return 1;
}



#endif
