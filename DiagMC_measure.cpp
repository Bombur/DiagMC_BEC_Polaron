#include "DiagMC.h"
  
int DiagMC::measure() {
  bool reduc= diag.is_reducible();
  int cur_order = diag.get_order();
  
#ifdef SELFENERGY
  if (reduc && cur_order>1) {	// In SELFENERGY sampling this should not be possible
	return 0;				// This will produce an error in the Test Function 
  }
#endif  

//Fake weight reweighting
  double cor = 1./pow(fw, cur_order);  
  
  
// We are in Zero Order
  if (cur_order == 0) {
#ifndef SELFENERGY
	Data(taumap.bin(diag.get_tau()), 0) += cor/fwzero/fwone; //0. Order contributes to All Orders for Green Sampling
#endif
	Data(taumap.bin(diag.get_tau()), 1) += cor/fwzero/fwone; //Zero Order container
  }
  

// We are in First Order
  if (cur_order ==1) {
	int bintau = taumap.bin(diag.get_tau());
	if(bintau > taumap.taubin || bintau <0){
	  std::cout << diag.get_tau() << '\t' << bintau <<std::endl;
	  assert(0);
	}
#ifdef SELFENERGY   // The First order can be measured as g0seg0 or g0se in Self ENergy Sampling
	if (!fog0seg0) {bintau = taumap.bin(diag.get_tinit(2));}
#else
	Data(bintau, 0) += cor/fwone; //1. Order contributes to All Orders for Green Sampling
#endif
	double tmp = Data.col(2).sum();
	Data(bintau, 2) += cor/fwone; //First Order container
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
  }
  
  //Norm Diagram
  if (cur_order == minord) {
	if (cur_order == 0) {
	  Norms(taumap.bin(diag.get_tau())) += cor;
	} else {
	  Norms(taumap.bin(diag.get_tinit(2*cur_order))) += cor;
	}
  }
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
  if (cur_order == 1 || (cur_order > 1  && !reduc)) {
	meas_Epol(cor);
	meas_SE(cor);
	meas_G0SEiw(cor);
  }
#else
  if ((cur_order > minord) && (cur_order <= maxord)) {
	meas_Epol(cor);
	meas_SE(cor);
	meas_G0SEiw(cor);
  }
#endif

//Tests
//for (int vert=1; vert<(2*diag.get_order()+1); vert++){
  //if (diag.get_link(vert+1) == (vert-1)) {testhisto(taumap.bin(diag.get_tfin(vert)- diag.get_tinit(vert)),2)+=1;}
//}
//Orderstats
  if (cur_order < orderstat.size()){
	orderstat(cur_order) +=1;
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
