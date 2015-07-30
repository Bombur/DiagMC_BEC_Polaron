#include "DiagMC.h"

void DiagMC::status() {
  //mean(Data);
  //bin_ana(Data);
  
  //printmean();
  printstats();
}

double DiagMC::mean(const VectorXd & measdata) {
  int N=measdata.size();
  anabuffer(0) = measdata.dot(VectorXd::Ones(N)) / N;
  return anabuffer(0);
}

void DiagMC::bin_ana(const VectorXd & measdata) {
  int N = measdata.size();
  double nB= N / binl;
  if (N%binl !=0) {std::cout<< "Warning! Measurements and Binning Length do not fit!"<< std::endl;}
  
  double Error=0;
  double Stddev=0;
	
  //binning error
  for (int n=0; n<nB; n++) {
	double binmean= measdata.segment(n*binl, binl).dot(VectorXd::Ones(binl)) / double(binl);
	Error += pow(binmean-anabuffer(0), 2)/double(nB*(nB-1));
  }
  anabuffer(1) = sqrt(Error);

	
  //integrated correlation time
  for (int j=0; j<N; j++){
	Stddev += pow(double(measdata(j)-anabuffer(0)), 2)/double(N-1);
  }
  anabuffer(2) = double(N)/2 * Error/Stddev;
}