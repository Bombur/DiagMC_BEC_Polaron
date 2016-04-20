#include "DiagMC.h"

//Estimators
//Polaron Energy
void DiagMC::meas_Epol(const double & cor){
  double Eptau = diag.get_tinit(2*diag.get_order()) - diag.get_tinit(1);
  if (diag.get_order() == 1) {
	//Epol.row(taumap.bin(Eptau)) += ((ws*Eptau - mu*Eptau).unaryExpr(std::ptr_fun(expfun)) *cor/fwone);
  } else {
	//std::cout << Eptau << '\t' << taumap.bin(Eptau) << '\n' << Epol << std::endl;
	Epol.row(taumap.bin(Eptau)) += ((ws*Eptau - mu*Eptau).unaryExpr(std::ptr_fun(expfun)) *cor);
  }
}

//Self Energy
void DiagMC::meas_SE(const double & cor){
  double SEtau = diag.get_tinit(2*diag.get_order()) - diag.get_tinit(1);
  if (diag.get_order() == 1) {
	SE(taumap.bin(SEtau), 0) += (cor/fwone); //all Sigma
	SE(taumap.bin(SEtau), 1) += (cor/fwone); //Sigma order ==1

  } else if (diag.get_order() > 1){
	SE(taumap.bin(SEtau), 0) += cor; //all Sigma
	SE(taumap.bin(SEtau), 2) += cor; //Sigma Order >1
  }
}

//G0SE in Matsubara
void DiagMC::meas_G0SEiw(const double & cor){
  double G0SEtau = diag.get_tinit(2*diag.get_order());
  const std::complex<double> i(0, 1);
  if (diag.get_order() ==1) {
	for (int it = 0; it < wbin; it ++) {
	  G0SEiw(it, 0) += cor/fwone * std::exp((i * static_cast<double>(it)/static_cast<double>(wbin)*wmax) * G0SEtau);
	  G0SEiw(it, 1) += cor/fwone * std::exp((i * static_cast<double>(it)/static_cast<double>(wbin)*wmax)* G0SEtau);
	}
  } else if (diag.get_order() >1){
	for (int it = 0; it < wbin; it ++) {
	  G0SEiw(it, 0) += cor * std::exp((i * static_cast<double>(it)/static_cast<double>(wbin)*wmax) * G0SEtau);
	}
  }
}