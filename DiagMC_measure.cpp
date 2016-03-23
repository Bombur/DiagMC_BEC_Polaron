#include "DiagMC.h"
  
int DiagMC::measure(const int & ordstp) {
//Fake weight reweighting
  double cor = 1./pow(fw, diag.get_order());  

#ifdef SELFENERGY
  if (diag.get_order() > 1) {
	if (diag.is_reducible()) {
	  std::cout << 1<< std::endl;
	  return 0;}
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
	Data(taumap.bin(diag.get_tinit(2)), 2) += cor;
#else
	Data(taumap.bin(diag.get_tau()), 2) += corbefore;
#endif
  }
  
  
//Tau to be measured
  double tau =0.;
  if (diag.get_order() > 1) {
#ifdef SIGMA
	tau = diag.get_tinit(2*diag.get_order())- diag.get_tinit(1); //to measur SE on its own
#else
	tau = diag.get_tinit(2*diag.get_order());
#endif
  } else {
	tau = diag.get_tau();
  }
  
  //second order G0SE
  if (diag.get_order() == 2)  {
	if (vsq(diag.get_p(1)-diag.get_p(3)) < 0.00000000000001){
	  //Diagram 2a
	  Data(taumap.bin(tau), 3) += cor;
	} else {
	  //Diagram 2b
	  Data(taumap.bin(tau), 4) += cor;
	}
  }
  
#ifndef SECUMUL
  //all G0SE
  if (diag.get_order() > 1 ){
	Data(taumap.bin(tau), 0) += cor;
  }
  
#else //if it's cumulative we do not need to measure all diagrams at once
    //Inbetween Diagram
  if ((diag.get_order() > 1) && (diag.get_order() > minord) && (diag.get_order() <= maxord)) {
	SEib(taumap.bin(tau)) += cor;
  }
	
	//Norm Diagram
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
	if(vsq(diag.get_p(2)-diag.get_p(0)) < 0.00000000000001) {
	  Data(taumap.bin(diag.get_tau()), 3) += cor;
	}
	Data(taumap.bin(diag.get_tau()), 4) += cor;
  }
#endif

//Estimator measurement. We need usually Self Energy sampling for it
  if (diag.get_order() == 1 || (diag.get_order() > 1  && !diag.is_reducible())) {
	meas_Epol(cor);
  }

//q Histogram
  for (auto vert=1; vert < 2*diag.get_order(); vert++){
	if (diag.get_link(vert) > vert) {
		double qtmp = sqrt(vsq(diag.get_p(vert-1)-diag.get_p(vert)));
		if (qtmp < qc) {
			qstat(static_cast<int>((qtmp/qc*static_cast<double>(qstat.size()))-0.5))+=1;
		}
	}
  }

  if (diag.get_order() < orderstat.size()){
	orderstat(diag.get_order()) +=1;
  }
  return 1;
}
