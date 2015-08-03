#include "DiagMC.h"

void DiagMC::write() {
  try {
	if (!(Datafile.is_open())) {throw openwritefile();}
	//long long tmp=Datafile.tellp();
	if (Datafile_pos!=0) {throw filepos();}
	
	MatrixXd output(taubin, 5);
	
	int CG0p = 0;
	for (int i=0; i< taubin; i++) {
	  CG0p += Data(i, 2);
	  output(i, 0) = i*taumax/taubin;
	}
	output.rightCols(4)=Data*(G0p/CG0p);
	
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



void DiagMC::printmean() {
  std::cout << "G(p)" << '\t' << "Binning Error" << '\t' << "Integrated Correlation Time" << '\n';
  std::cout << anabuffer.transpose() << std::endl;
}

void DiagMC::printstats() {
  std::cout <<'\n'<< "*************************************" << std::endl;
  std::cout << "1:ATTEMPTED" << '\n' << "2:POSSIBLE" << '\n' << "3:REJECTED" << '\n' << "4:ACCEPTED" << '\n' << "5:ACCEPTANCE RATIO POSSIBLE" << '\n' << "6:ACCEPTANCE RATIO TOTAL" << '\n';
  std::cout << stats << '\n'<< std::endl;
}
/*
void DiagMC::printallvariables() {
  std::cout<< N << '\n' << E << ' '<< m << '\n'<< beta <<'\n' << "no sigma" <<'\n'<< "no Nbr" <<'\n' << '\n';
  std::cout << Meas_its << ' ' << Test_its << ' '<< Write_its << '\n' << RunTime << '\n'<< irnd() <<'\n' << drnd()<<'\n'<<'\n' << binl<< '\n'<< meanbuffer <<'\n' << stats <<'\n'<< Datafile_pos <<'\n' << Statfile_pos << '\n';
}
*/