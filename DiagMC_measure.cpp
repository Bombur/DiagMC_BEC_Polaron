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
	Data((int)(diag.get_tau()/taumax*static_cast<double>(taubin)), 1) += cor;
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
	Data((int)(diag.get_tinit(2)/taumax*static_cast<double>(taubin)), 2) += cor;
#else
	Data((int)(diag.get_tau()/taumax*static_cast<double>(taubin)), 2) += corbefore;
#endif
  }
  
  //second order G0SE
  if (diag.get_order() == 2)  {
	if (vsq(diag.get_p(1)-diag.get_p(3)) < 0.00000000000001){
	  //Diagram 2a
	  Data((int)(diag.get_tinit(2*diag.get_order())/taumax*static_cast<double>(taubin)), 3) += cor;
	} else {
	  //Diagram 2b
	  Data((int)(diag.get_tinit(2*diag.get_order())/taumax*static_cast<double>(taubin)), 4) += cor;
	}
  }
  
#ifndef SECUMUL
  //all G0SE
  if (diag.get_order() > 1 ){
	Data((int)(diag.get_tinit(2*diag.get_order())/taumax*static_cast<double>(taubin)), 0) += cor;
  }
  
#else //if it's cumulative we do not need to measure all diagrams at once
    //Norm Diagram
  if (diag.get_order() == minord) {
	Norms((int)(diag.get_tinit(2*diag.get_order())/taumax*static_cast<double>(taubin))) += cor;
  }
  
  //End Diagram
  if (diag.get_order() == maxord) {
	Ends((int)(diag.get_tinit(2*diag.get_order())/taumax*static_cast<double>(taubin))) += cor;
  }
  
  //Inbetween Diagram
  if (diag.get_order()>1){
	if (diag.get_order() > minord && diag.get_order() <= maxord) {
	  SEib((int)(diag.get_tinit(2*diag.get_order())/taumax*static_cast<double>(taubin))) += cor;
	}
  }
#endif
  
#else
  Data((int)(diag.get_tau()/taumax*static_cast<double>(taubin)), 0) += cor;
  if (diag.get_order() < 2) {Data((int)(diag.get_tau()/taumax*static_cast<double>(taubin)), diag.get_order()+1) += cor;}
  if (diag.get_order() ==2){
	if(vsq(diag.get_q(2)) < 0.00000000000001) {
	  Data((int)(diag.get_tau()/taumax*static_cast<double>(taubin)), 3) += cor;
	}
	Data((int)(diag.get_tau()/taumax*static_cast<double>(taubin)), 4) += cor;
  }
#endif


  if (diag.get_order() < orderstat.size()){
	orderstat(diag.get_order()) +=1;
  }
  return 1;
}