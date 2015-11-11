#include "DiagMC.h"
 
void DiagMC::write() {
  try {
	ArrayXXd  output(taubin, 5);
	
	int CG0p = 0;
	for (int i=0; i< taubin; i++) {
	  CG0p += Data(i, 1);
	  output(i, 0) = (static_cast<double>(i)+0.5)*taumax/static_cast<double>(taubin);
	}
	output.rightCols(4)=Data.cast<double>()*(G0p/static_cast<double>(CG0p))*static_cast<double>(taubin)/taumax;
  }
  catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

  
ArrayXXd DiagMC::get_Data() {
  ArrayXXd  output(taubin, 6);
	
  int CG0p = 0;
  for (int i=0; i< taubin; i++) {
	CG0p += Data(i, 1);
	output(i, 0) = (static_cast<double>(i)+0.5)*taumax/static_cast<double>(taubin);
#ifndef SELFENERGY
	output(i,1) = static_cast<double>(Data(i, 0))-(1.-1./fw)*static_cast<double>(Data(i,1));
#endif
  }
 

#ifdef SELFENERGY
  output.rightCols(5)=Data.cast<double>()/static_cast<double>(CG0p)*fw*static_cast<double>(taubin)/taumax;
  //zero order (fake check)
  output.col(2)*= G0p/fw;
  output.col(3)*= G0p; //first order  G0SEG0
#else
  output.col(1) *= (G0p/static_cast<double>(CG0p))*fw*static_cast<double>(taubin)/taumax;
  output.rightCols(4)=Data.rightCols(4).cast<double>()*(G0p/static_cast<double>(CG0p))*fw*static_cast<double>(taubin)/taumax;
//order fake
  output.col(2)/= fw;
#endif	
  return output;	
}
