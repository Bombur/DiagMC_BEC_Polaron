#include "DiagMC.h"

class openwritefile: public std::exception {
  virtual const char* what() const throw()
  {
    return "Write File not open!";
  }
};


class filepos: public std::exception {
  virtual const char* what() const throw()
  {
    return "Positions in Write File do not fit!";
  }
};

void DiagMC::write() {
  try {
	if (!(Datafile.is_open())) {throw openwritefile();}
	//long long tmp=Datafile.tellp();
	if (Datafile_pos!=0) {throw filepos();}
	
	MatrixXd  output(taubin, 5);
	
	int CG0p = 0;
	for (int i=0; i< taubin; i++) {
	  CG0p += Data(i, 1);
	  output(i, 0) = (i+0.5)*taumax/taubin;
	}
	output.rightCols(4)=Data.cast<double>()*(G0p/CG0p);
	
	Datafile.seekp(Datafile_pos);
	Datafile << output << '\n';
	Datafile.flush();
  }
  catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

void DiagMC::Stattofile() {
  try {
	if (!(Statsfile.is_open())) {throw openwritefile();}
	long long tmp=Statsfile.tellp();
	if (Statfile_pos!=tmp) {throw filepos();}
	
	Statsfile << stats << '\n';
	Statsfile.flush();
	Statfile_pos=Statsfile.tellp(); 
  }
  catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

void DiagMC::status() {
  for (int i = 0; i < 8; i++) {
	stats(i, 4) = stats(i, 3)/stats(i, 1);
	stats(i, 5) = stats(i, 3)/stats(i, 0);
  }
  
  printstats();  
}

void DiagMC::printmean() {
  std::cout << "G(p)" << '\t' << "Binning Error" << '\t' << "Integrated Correlation Time" << '\n';
  std::cout << anabuffer.transpose() << std::endl;
}

void DiagMC::printstats() {
  std::cout <<'\n'<< "*************************************" << std::endl;
  std::cout << "COLUMNS" << '\n' << "1:ATTEMPTED" << '\n' << "2:POSSIBLE" << '\n' << "3:REJECTED" << '\n' << "4:ACCEPTED" << '\n' << "5:ACCEPTANCE RATIO POSSIBLE" << '\n' << "6:ACCEPTANCE RATIO TOTAL" << '\n';
  std::cout << "*************************************" << std::endl;
  std::cout << "CHANGE OF TAU: " << '\t' << stats.topRows(1) << '\n' << "CT IN HO: " << '\t' << stats.block(7, 0, 1, 6) << '\n' << "INSERT: " << '\t' << stats.block(1, 0, 1, 6) << '\n' << "REMOVE: " << '\t' << stats.block(2, 0, 1, 6) << '\n' << "SWAP:   " << '\t' << stats.block(3, 0, 1, 6) << '\n'<< "SWAPoocc: " << '\t' << stats.block(4, 0, 1, 6) << '\n'<< "SWAPoc: " << '\t' << stats.block(5, 0, 1, 6) << '\n'<< "SWAPco: " << '\t' << stats.block(6, 0, 1, 6) << '\n' <<std::endl;
}

/*
for (int i=0; i< taubin; i++) {
	  for (int i2 = 0; i2 < Data.size() ; i2++) {
		output(i, i2+1)= Data(i, i2)*(G0p/CG0p);
	  }
	}
	
	*/