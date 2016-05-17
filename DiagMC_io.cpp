#include "DiagMC.h"

ArrayXXd DiagMC::get_Data() {
  ArrayXXd output(taumap.taubin, 6);
  output.col(0) = taumap.print();    // Taus
	
  double CG0p = Data.col(1).sum();  //counts in 0 Order
  output.rightCols(5) = taumap.norm_table().replicate<1,5>();   // delta T pro bin
  output.rightCols(5) *= Data/CG0p *G0ptmima; 
#ifdef SELFENERGY
  output.rightCols(5)/=G0p1arm;   //In the case of Green Sampling we add the last G0 to the data
#endif
   
  //Zero Order
#ifdef SELFENERGY
  output.col(2)*= G0p1arm;
#endif

   //First Order
#ifdef SELFENERGY
  if (fog0seg0) { // if you want to sample g0seg0 in Selfenergy sampling
	output.col(3)*= G0p1arm; 
  }
#endif

  return output;	
}

//Polaron Energy Estimator
ArrayXXd DiagMC::get_Eptest() {
  ArrayXXd output = Epol;
  for (int i=0 ; i< Epol.cols() ; i++) {
	output.col(i+1) *= taumap.norm_table();
  }
  return output *(G0ptmima/G0p2arms) / Data.col(1).sum();
}

ArrayXd DiagMC::get_Ep() {
  ArrayXXd output = Epol;
#ifdef SECUMUL
  return output.colwise().sum() *(G0ptmima/G0p2arms) *pref_calc();
#else
  return output.colwise().sum() *(G0ptmima/G0p2arms) /Data.col(1).sum();
#endif
}

//Self Energy Estimator
ArrayXXd DiagMC::get_SE() {
  ArrayXXd output(taumap.taubin, SE.cols()+1);
  output.col(0) = taumap.print();
#ifdef SECUMUL
  output.rightCols(SE.cols()) = SE *(G0ptmima/G0p2arms) *pref_calc();
#else
  output.rightCols(SE.cols()) = SE *(G0ptmima/G0p2arms) /Data.col(1).sum();
#endif
  for (int i=0 ; i< SE.cols() ; i++) {
	output.col(i+1) *= taumap.norm_table();
  }
  return output;
}

//G0SE in Matsubara
ArrayXXcd DiagMC::get_G0SEiw() {
  ArrayXXcd output(wbin, G0SEiw.cols()+1);
  for (int it =0; it<wbin; it++) {
	output(it,0) = static_cast<double>(it)/static_cast<double>(wbin)*wmax;
  }
#ifdef SECUMUL
  output.rightCols(G0SEiw.cols()) = G0SEiw /G0p1arm *pref_calc();
#else
  output.rightCols(G0SEiw.cols()) = G0SEiw /G0p1arm /Data.col(1).sum();
#endif
  return output;  
}



//tau Histogram
ArrayXXd DiagMC::get_taus() {
  ArrayXXd output(taumap.taubin, tstat.cols()+1);
  output.col(0) = taumap.print();
  output.rightCols(tstat.cols()) = tstat;
  for (int i=0 ; i< tstat.cols() ; i++) {
	output.col(i+1) *= taumap.norm_table();
  }
  return output;
}

ArrayXXd DiagMC::get_testhisto() {
  return testhisto;
}