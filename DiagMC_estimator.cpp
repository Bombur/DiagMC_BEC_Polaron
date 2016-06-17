#include "DiagMC.h"
 
//Estimators
//Polaron Energy
void DiagMC::fill_Epcont(const double & cor){
  double Eptau = diag.get_tinit(2*diag.get_order()) - diag.get_tinit(1);
  for (int i=0; i<ws.size(); ++i){
	  Eptmp[i] = (expfun(ws(i)*Eptau - mu*Eptau) *cor);
  }
}

void DiagMC::meas_Epol(){
  if (((diag.get_order() == 1) && fst_Ep_meas) || (diag.get_order()>1)) {
	Epol += Eptmp;
  }
}

//Self Energy
void DiagMC::meas_SE(const double & cor){
  double SEtaubin = taumap.bin(diag.get_tinit(2*diag.get_order()) - diag.get_tinit(1));
  if (diag.get_order() == 1) {
	SE(SEtaubin, 0) += cor; //all Sigma
	SE(SEtaubin, 1) += cor; //Sigma order ==1

  } else if (diag.get_order() > 1){
	SE(SEtaubin, 0) += cor; //all Sigma
	SE(SEtaubin, 2) += cor; //Sigma Order >1
  }
}

//G0SE in Matsubara
void DiagMC::meas_G0SEiw(const double & cor){
  double G0SEtau = diag.get_tinit(2*diag.get_order());
  const std::complex<double> i(0, 1);
  if (diag.get_order() ==1) {
	for (int it = 0; it < wbin; it ++) {
	  G0SEiw(it, 0) += cor * std::exp((i * static_cast<double>(it)/static_cast<double>(wbin)*wmax) * G0SEtau);
	  G0SEiw(it, 1) += cor * std::exp((i * static_cast<double>(it)/static_cast<double>(wbin)*wmax)* G0SEtau);
	}
  } else if (diag.get_order() >1){
	for (int it = 0; it < wbin; it ++) {
	  G0SEiw(it, 0) += cor * std::exp((i * static_cast<double>(it)/static_cast<double>(wbin)*wmax) * G0SEtau);
	}
  }
}


//binning
void DiagMC::bin_Epol(){
  	if (((diag.get_order() == 1) && fst_Ep_meas) || (diag.get_order()>1)) {
		Epbin["step"+std::to_string(ordstep)] << Eptmp;
	}
}

//measure Ep as integral
void DiagMC::meas_Ep_intv() {
#ifndef SECUMUL
  double denominator = (Data.col(1).sum()-last_g0_count);
#else
  double denominator = (Norms.sum()-last_g0_count);
#endif
  overflowstat(13,0)+=1;
  if (denominator < 1e-8) {overflowstat(13,1)+=1;}
  else {
	Ep_intv["step"+std::to_string(ordstep)] << (Epol-last_measured_Epol)*get_Ep_norm()/denominator;
  }
  last_g0_count += denominator;
  last_measured_Epol = Epol;
}

void DiagMC::meas_G0tau0() {
  double denominator = (Data.col(1).sum()-last_g0_count);
  if (!(denominator < 1e-8)) {
	ArrayXd tmp= taumap.norm_table();
	G0tau0["tau0"] << (Data(0,1)-last_count_g0attau0)*G0ptmima/denominator*tmp(0);
  }
  last_count_g0attau0 = Data(0,1);
}

void DiagMC::meas_G0tau0_each() {
  if (taumap.bin(diag.get_tau()) == 0){
	G0tau0_each["tau0"]<<1.;
  } else {
	G0tau0_each["tau0"]<<0.;
  }
}
 
void DiagMC::meas_ordesti(){
	ordesti["step"+std::to_string(ordstep)] << diag.get_order();
}
