#include "DiagMC.h"

/*
void DiagMC::write() {
  try {
	ArrayXXd  output(taumap.taubin, 5);
	
	double CG0p = 0;
	for (int i=0; i< taumap.taubin; i++) {
	  CG0p += Data(i, 1);
	  output(i, 0) = (static_cast<double>(i)+0.5)*taumap.taumax/static_cast<double>(taumap.taubin);
	}
	output.rightCols(4)=Data*(G0p/CG0p)*static_cast<double>(taumap.taubin)/taumap.taumax;
  }
  catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}
*/

  
ArrayXXd DiagMC::get_Data() {
  ArrayXXd  output(taumap.taubin, 6);
	
  double CG0p = Data.col(1).sum()/fwzero/fwone;
  output.col(0) = taumap.print();
  output.rightCols(5) = taumap.norm_table().replicate<1,5>();
  output.col(2) /= (fwzero*fwone);
  output.col(3) /= fwone;

#ifdef SELFENERGY
//if for how to reweight the data
#ifdef MEASREWEIGHT
  output.rightCols(5)*=Data*(G0p/CG0p);
#else
//Copy everything  
  output.rightCols(5)*=Data*(1./CG0p);
  
//Reweights and Norms  
  //zero order (fake check)
  output.col(2)*= G0p;
#ifndef FOG0SE
  //first order  G0SEG0
  output.col(3)*= G0p; 
#endif
  
#ifdef SIGMA
  //Orders >2 we measure just SE
  //all_orders
  output.col(1) /= G0p;
  //2nd orders
  output.rightCols(2) /= G0p;
#endif
    
#endif

#else
  output.rightCols(5)*=Data*(G0p/CG0p);
#endif	
  return output;	
}
