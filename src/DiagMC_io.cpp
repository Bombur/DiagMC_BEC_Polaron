#include "DiagMC.h"

void DiagMC::write() {
  try {
	if (!(Datafile.is_open())) {throw openwritefile();}
	long long tmp=Datafile.tellp();
	if (Datafile_pos!=tmp) {throw filepos();}
	
	Datafile << Data << '\n';
	Datafile.flush();
	Datafile_pos=Datafile.tellp(); 
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

void Ising_markov_method::printmean(const std::string Observable) {
  std::cout << Observable << '\t' << "Binning Error" << '\t' << "Integrated Correlation Time" << '\n';
  std::cout << meanbuffer.transpose() << std::endl;
}

void Ising_markov_method::printstats() {
  std::cout <<'\n'<< "*************************************" << std::endl;
  std::cout << "1:ATTEMPTED" << '\n' << "2:POSSIBLE" << '\n' << "3:REJECTED" << '\n' << "4:ACCEPTED" << '\n' << "5:ACCEPTANCE RATIO POSSIBLE" << '\n' << "6:ACCEPTANCE RATIO TOTAL" << '\n';
  std::cout << stats << '\n'<< std::endl;
}

void Ising_markov_method::printallvariables() {
  std::cout<< N << '\n' << E << ' '<< m << '\n'<< beta <<'\n' << "no sigma" <<'\n'<< "no Nbr" <<'\n' << '\n';
  std::cout << Meas_its << ' ' << Test_its << ' '<< Write_its << '\n' << RunTime << '\n'<< irnd() <<'\n' << drnd()<<'\n'<<'\n' << binl<< '\n'<< meanbuffer <<'\n' << stats <<'\n'<< Datafile_pos <<'\n' << Statfile_pos << '\n';
}

*/