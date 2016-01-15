#include "DiagMC.h"
  
int DiagMC::measure(const int & ordstp) {
//Fake weight reweighting
  double cor = 1./pow(fw, diag.get_order());  
  
#ifdef SELFENERGY
  if (diag.get_order() > 1) {
	if (diag.is_reducible()) {return 0;}
  }
  //zero order (fake check)
  if (diag.get_order() == 0) {
	Data(taumap.bin(diag.get_tau()), 1) += cor;
  }
  
//To measure G0Se from full Diagram sampling we have to reweight the measurement
//In case of else we have to reweight in IO
  double corbefore = cor;
#ifdef MEASREWEIGHT
  cor /= G0el(diag.get_p(0), diag.get_tau(), diag.get_tinit(2*diag.get_order()));
#endif
  
 //first order G0Se
  if (diag.get_order() == 1) {
//Measure either G0SE or G0SEG0
#ifdef FOG0SE
	Data(taumap.bin(diag.get_tau()), 2) += cor;
#else
	Data(taumap.bin(diag.get_tau()), 2) += corbefore;
#endif
  }
  
  //second order G0SE
  if (diag.get_order() == 2)  {
	if (vsq(diag.get_p(1)-diag.get_p(3)) < 0.00000000000001){
	  //Diagram 2a
	  Data(taumap.bin(diag.get_tinit(4)), 3) += cor;
	} else {
	  //Diagram 2b
	  Data(taumap.bin(diag.get_tinit(4)), 4) += cor;
	}
  }
  
  

#ifndef SECUMUL
  //all G0SE
  if (diag.get_order() > 1 ){
	Data(taumap.bin(diag.get_tinit(2*diag.get_order())), 0) += cor;
  }
  
#else //if it's cumulative we do not need to measure all diagrams at once
    //Norm Diagram
	double tau;
	if (diag.get_order() > 1) {
	  tau = diag.get_tinit(2*diag.get_order());
	  //Inbetween Diagram
	  if (diag.get_order() > minord && diag.get_order() <= maxord) {
		SEib(taumap.bin(diag.get_tinit(2*diag.get_order()))) += cor;
	  }
	} else {
	  tau = diag.get_tau();
	}
	
	if (diag.get_order() == minord) {
	  Norms(taumap.bin(tau)) += cor;
	}
  
	//End Diagram
	if (diag.get_order() == maxord) {
	  Ends(taumap.bin(tau)) += cor;
	}
#endif

  
#else
  Data(taumap.bin(diag.get_tau()), 0) += cor;
  if (diag.get_order() < 2) {Data(taumap.bin(diag.get_tau()), diag.get_order()+1) += cor;}
  if (diag.get_order() ==2){
	if(vsq(diag.get_q(2)) < 0.00000000000001) {
	  Data(taumap.bin(diag.get_tau()), 3) += cor;
	}
	Data(taumap.bin(diag.get_tau()), 4) += cor;
  }
#endif


  if (diag.get_order() < orderstat.size()){
	orderstat(diag.get_order()) +=1;
  }
  return 1;
}