#include "DiagMC.h"
 
//Estimators
//Polaron Energy
void DiagMC::meas_Epol(const double & cor){
  double Eptau = diag.get_tinit(2*diag.get_order()) - diag.get_tinit(1);
  if (diag.get_order() == 1) {
	Epol.row(taumap.bin(Eptau)) += ((ws*Eptau - mu*Eptau).unaryExpr(std::ptr_fun(expfun)) *cor);
  } else {
	Epol.row(taumap.bin(Eptau)) += ((ws*Eptau - mu*Eptau).unaryExpr(std::ptr_fun(expfun)) *cor);
  }
}

//Self Energy
void DiagMC::meas_SE(const double & cor){
  double SEtaubin = taumap.bin(diag.get_tinit(2*diag.get_order()) - diag.get_tinit(1));
  if (diag.get_order() == 1) {
	SE(SEtaubin, 0) += (cor); //all Sigma
	SE(SEtaubin, 1) += (cor); //Sigma order ==1

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
