#include "DiagMC.h"
  
int DiagMC::measure(const int & ordstp) {
#ifdef SELFENERGY
  if (diag.get_order() > 1) {
	if (diag.is_reducible()) {return 0;}
  }
  //zero order (fake check)
  if (diag.get_order() == 0) {
	Data((int)(diag.get_tau()/taumax*static_cast<double>(taubin)), 1) += 1;
  }
  //first order G0SeG0
  if (diag.get_order() == 1) {
	Data((int)(diag.get_tau()/taumax*static_cast<double>(taubin)), 2) += 1;
  }
  //second order G0SE
  if (diag.get_order() == 2)  {
	if (vsq(diag.get_p(1)-diag.get_p(3)) < 0.00000000000001){
	  //Diagram 2a
	  Data((int)(diag.get_tinit(2*diag.get_order())/taumax*static_cast<double>(taubin)), 3) += 1;
	} else {
	  //Diagram 2b
	  Data((int)(diag.get_tinit(2*diag.get_order())/taumax*static_cast<double>(taubin)), 4) += 1;
	}
  }
#ifndef SECUMUL
  //all G0SE
  if (diag.get_order() > 1 ){
	Data((int)(diag.get_tinit(2*diag.get_order())/taumax*static_cast<double>(taubin)), 0) +=1;
  }
  
#else //if it's cumulative we do not need to measure all diagrams at once
  if (diag.get_order()>1){
  //Norm Diagram
  if (diag.get_order() == minord) {
	Norms((int)(diag.get_tinit(2*diag.get_order())/taumax*static_cast<double>(taubin))) += 1;
  }
  
  //End Diagram
  if (diag.get_order() == maxord) {
	Ends((int)(diag.get_tinit(2*diag.get_order())/taumax*static_cast<double>(taubin))) += 1;
  }
  
  //Inbetween Diagram
  if (diag.get_order() > minord && diag.get_order() <= maxord) {
	SEib((int)(diag.get_tinit(2*diag.get_order())/taumax*static_cast<double>(taubin))) += 1;
  }
  }
#endif
  
#else
  Data((int)(diag.get_tau()/taumax*static_cast<double>(taubin)), 0) +=1;
  if (diag.get_order() < 2) {Data((int)(diag.get_tau()/taumax*static_cast<double>(taubin)), diag.get_order()+1) += 1;}
  if (diag.get_order() ==2){
	if(vsq(diag.get_q(2)) < 0.00000000000001) {
	  Data((int)(diag.get_tau()/taumax*static_cast<double>(taubin)), 3) += 1;
	}
	Data((int)(diag.get_tau()/taumax*static_cast<double>(taubin)), 4) += 1;
  }
#endif
  if (diag.get_order() < orderstat.size()){
	orderstat(diag.get_order()) +=1;
  }
  return 1;
}