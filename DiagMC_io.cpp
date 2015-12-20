#include "DiagMC.h"
 
void DiagMC::write() {
  try {
	ArrayXXd  output(taubin, 5);
	
	double CG0p = 0;
	for (int i=0; i< taubin; i++) {
	  CG0p += Data(i, 1);
	  output(i, 0) = (static_cast<double>(i)+0.5)*taumax/static_cast<double>(taubin);
	}
	output.rightCols(4)=Data*(G0p/CG0p)*static_cast<double>(taubin)/taumax;
  }
  catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

  
ArrayXXd DiagMC::get_Data() {
  ArrayXXd  output(taubin, 6);
	
  double CG0p = 0;
  for (int i=0; i< taubin; i++) {
	CG0p += Data(i, 1);
	output(i, 0) = (static_cast<double>(i)+0.5)*taumax/static_cast<double>(taubin);
  }
 

#ifdef SELFENERGY
//if for how to reweight the data
#ifdef MEASREWEIGHT
  output.rightCols(5)=Data*((G0p/CG0p)*static_cast<double>(taubin)/taumax);
  std::cout << testg0p/count << std::endl;
#else
  output.rightCols(5)=Data*(1./CG0p*static_cast<double>(taubin)/taumax);
  //zero order (fake check)
  output.col(2)*= G0p;
#ifndef FOG0SE
  //first order  G0SEG0 
  output.col(3)*= G0p; 
#endif
#endif

#else
  output.rightCols(5)=Data*(G0p/CG0p)*static_cast<double>(taubin)/taumax;
#endif	
  return output;	
}
