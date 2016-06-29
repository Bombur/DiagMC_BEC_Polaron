#include "DiagMC.h"
  
int DiagMC::measure() {
#ifndef SELFENERGY
bool reduc= diag.is_reducible();
#endif 

int cur_order = diag.get_order();

//Fake weight reweighting
double cor = fw_for_meas(cur_order);  

// We are in Zero Order
if (cur_order == 0) {
	int bintau = taumap.bin(diag.get_tau());
#ifndef SELFENERGY
	Data(bintau, 0) += cor; //0. Order contributes to All Orders for Green Sampling
#endif
	Data(bintau, 1) += cor; //Zero Order container
#ifndef SECUMUL
	counts["norm0"] << cor;
} else if (ordstep ==0) {counts["norm0"] << 0.;}
#else
}
#endif

// We are in First Order
if (cur_order ==1) {
	int bintau = taumap.bin(diag.get_tau());
#ifdef SELFENERGY   // The First order can be measured as g0seg0 or g0se in Self ENergy Sampling
	if (!fog0seg0) {bintau = taumap.bin(diag.get_tinit(2));}
#else
	Data(bintau, 0) += cor; //1. Order contributes to All Orders for Green Sampling
#endif
	Data(bintau, 2) += cor; //First Order container
}
  
  
//We are in 2nd Order the contributions to >1 in Self Energy Sampling will be treated later
if (cur_order == 2) {
#ifdef SELFENERGY
	int bintau = taumap.bin(diag.get_tinit(4));
		//Rainbow
	if (vsq(diag.get_p(1)-diag.get_p(3)) < 0.00000000000001){
		Data(bintau, 3) += cor;
	} else {
		//Crossed
		Data(bintau, 4) += cor;
	}
#else // Greens Sampling
	int bintau = taumap.bin(diag.get_tau());
		//reducible 
	if (reduc) {
		Data(bintau, 3) += cor;
	} else if (!(vsq(diag.get_p(1)-diag.get_p(3)) < 0.00000000000001)) {
		  //crossed
		Data(bintau, 4) += cor;
	}
#endif
}
  

//We are in >1
#ifdef SELFENERGY
#ifdef SECUMUL //if it's cumulative Data.col(0)(all_orders) is replaced
  //Inbetween Diagram
if ((cur_order > 1) && (cur_order > minord) && (cur_order <= maxord)) {
	SEib(taumap.bin(diag.get_tinit(2*cur_order))) += cor;
}
	
  	//End Diagram
if (cur_order == maxord) {
	Ends(taumap.bin(diag.get_tinit(2*cur_order))) += cor;
	counts["end"+std::to_string(ordstep)] << cor;
} else {counts["end"+std::to_string(ordstep)] << 0.;}
  
  	//Norm Diagram
if (cur_order == minord) {
	if (cur_order == 0) {
		Norms(taumap.bin(diag.get_tau())) += cor;
	} else {
		Norms(taumap.bin(diag.get_tinit(2*cur_order))) += cor;
	}
	counts["norm"+std::to_string(ordstep)] << cor;
} else {counts["norm"+std::to_string(ordstep)] << 0.;}
#else 
  	//all G0SE
if (cur_order > 1 ){
	Data(taumap.bin(diag.get_tinit(2*cur_order)), 0) += cor;
}
#endif	
#else // Green Sampling
if (cur_order > 1 ){
	Data(taumap.bin(diag.get_tau()),0) += cor;
}
#endif

//Estimator measurement. We need usually Self Energy sampling for it
#ifndef SECUMUL
#ifndef SELFENERGY
if (cur_order == 1 || (cur_order > 1  && !reduc)) {
#else
if (cur_order > 0) {
#endif
#else
if ((cur_order > minord) && (cur_order <= maxord)) {
#endif
	fill_Epcont(cor);
	meas_Epol();
	meas_SE(cor);
	if (G0SEiw_meas) {meas_G0SEiw(cor);}
	//binning
	if (Ep_bin_each_step){bin_Epol();}
	meas_ordesti();
}
return 1;
}

void DiagMC::meas_histo(){
  //q Histogram an t Histogram
  int cur_order = diag.get_order(); 
  for (auto vert=1; vert < 2*cur_order; vert++){
	if (diag.get_link(vert) > vert) {
	  
	  //tau
	  double ttmp= diag.get_tinit(diag.get_link(vert)) - diag.get_tinit(vert);   //measure length of Phonon in tau
	  int tbin = taumap.bin(ttmp);
	  
	  if (diag.get_link(vert) == vert+1) {tstat(tbin,1)+=1;}
	  else if ((diag.get_link(vert) == vert + 2) && (diag.get_link(vert+1) == vert+3)){
		tstat(tbin,2) +=1;
	  } else if ((diag.get_link(vert) == vert + 3) && (diag.get_link(vert+1) == vert+2)){
		tstat(tbin,3) +=1;
	  }
	  tstat(tbin, 0)+=1;
		
	  //q
	  double qtmp = sqrt(vsq(diag.get_p(vert-1)-diag.get_p(vert)));
	  if (qtmp < qc) {
		int qbin = static_cast<int>((qtmp/qc*static_cast<double>(qstat.rows()))-0.5);
		qstat(qbin,0)+=1;
		if (diag.get_order()==1) {
		  qstat(qbin,1)+=1;
		}
	  }
	}
  }
}


void DiagMC::meas_ordstat() {
	//Orderstats
	assert(diag.get_order() < orderstat.size());
	orderstat(diag.get_order()) +=1;
}

void DiagMC::set_os_to_zero() {
	if ((orderstat(minord)> 1e-8) || (orderstat(minord+1)>1e-8)){orderstat.setZero(orderstat.size());}
}

double DiagMC::fw_for_meas(const int & ord, const bool & fw_ad) {
	double cor = 1.;		
	for (int i = minord; i < ord; i++){cor /= fw;}
	if (!fw_ad){
		if (ord < fwtab.size()) {for (int i = minord; i < ord; i++){cor /= fwtab[i];}}
		else if (minord < fwtab.size()) {for (int i = minord; i < fwtab.size(); i++){cor /= fwtab[i];}}
	}
#ifdef SECUMUL
	if (which_fw_ad == 2){
		assert((ord-minord) <= fw_vec.size());
   		for (int i = 0; i < (ord- minord); i++){cor /= fw_vec[i];}
	}
#endif
	return cor;
}
