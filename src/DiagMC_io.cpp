#include "DiagMC.h"

void DiagMC::write() {
  try {
	if (!(Datafile.is_open())) {throw openwritefile();}
	//long long tmp=Datafile.tellp();
	if (Datafile_pos!=0) {throw filepos();}
	
	Datafile.seekp(Datafile_pos);
	Datafile << Data*G0p/stats(0,1) << '\n';
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



/*
void Ising_markov_method::meantofile(std::ofstream &file, const long long writing_pos) {
  try {
	if (!(file.is_open())) {throw openwritefile();}
	long long tmp=file.tellp();
	if (writing_pos!=tmp) {throw filepos();}
	
	file << meanbuffer.transpose() << '\n';
	file.flush();
  }
  catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}
*/

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