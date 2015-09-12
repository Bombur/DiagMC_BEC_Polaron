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

void DiagMC::Stattofile(const VectorXd & timestat) {
  try {
	if (!(udsfile.is_open() && osfile.is_open() && tsfile.is_open())) {throw openwritefile();}
	if (uds_pos!=0 && os_pos!=0 && ts_pos!=0) {throw filepos();}
	
	udsfile << "Update Statistics" << '\n';
    udsfile << "*************************************" << '\n';
    udsfile << "COLUMNS" << '\n' << "1:ATTEMPTED" << '\n' << "2:POSSIBLE" << '\n' << "3:REJECTED" << '\n' << "4:ACCEPTED" << '\n' << "5:ACCEPTANCE RATIO POSSIBLE" << '\n' << "6:ACCEPTANCE RATIO TOTAL" << '\n';
    udsfile << "*************************************" << '\n';
    udsfile << "CHANGE TAU: " << '\t' << updatestat.topRows(1) << '\n' << "CT IN HO: " << '\t' << updatestat.block(7, 0, 1, 6) << '\n' << "INSERT: " << '\t' << updatestat.block(1, 0, 1, 6) << '\n' << "REMOVE: " << '\t' << updatestat.block(2, 0, 1, 6) << '\n' << "SWAP:   " << '\t' << updatestat.block(3, 0, 1, 6) << '\n'<< "SWAPoocc: " << '\t' << updatestat.block(4, 0, 1, 6) << '\n'<< "SWAPoc: " << '\t' << updatestat.block(5, 0, 1, 6) << '\n'<< "SWAPco: " << '\t' << updatestat.block(6, 0, 1, 6) << '\n' << "DQ:      " << '\t' << updatestat.block(8, 0, 1, 6)<< '\n' << '\n';
	udsfile.flush();
	udsfile.seekp(ts_pos);
	
	MatrixXi orderprint(orderstat.size(),2);
	orderprint << VectorXi::LinSpaced(orderstat.size(), 0, orderstat.size()-1), orderstat;
	
	osfile << "Order Statistics" << '\n';
	osfile << orderprint  << '\n';
	osfile.flush();
	osfile.seekp(ts_pos);
	
	tsfile << "Time Statistics [%]" << std:: endl;
	tsfile << "CHANGE TAU:" << '\t' << timestat(0) << '\n' << "CT HO:     " << '\t' <<timestat(1)<< '\n' <<  "INSERT: " << '\t' << timestat(2) << '\n' << "REMOVE: " << '\t' << timestat(3) << '\n' <<  "SWAP:     " << '\t' << timestat(4) << '\n' << "DQ:      " << '\t' << timestat(5) << '\n';		
	tsfile.flush();
	tsfile.seekp(ts_pos); 
  }
  catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

void DiagMC::status() {
  for (int i = 0; i < 9; i++) {
	updatestat(i, 4) = updatestat(i, 3)/updatestat(i, 1);
	updatestat(i, 5) = updatestat(i, 3)/updatestat(i, 0);
  }
  
  updatestat();
  orderstats();
}

void DiagMC::updatestat() {
  std::cout << "Update Statistics" << '\n';
  std::cout <<'\n'<< "*************************************" << std::endl;
  std::cout << "COLUMNS" << '\n' << "1:ATTEMPTED" << '\n' << "2:POSSIBLE" << '\n' << "3:REJECTED" << '\n' << "4:ACCEPTED" << '\n' << "5:ACCEPTANCE RATIO POSSIBLE" << '\n' << "6:ACCEPTANCE RATIO TOTAL" << '\n';
  std::cout << "*************************************" << std::endl;
  std::cout << "CHANGE OF TAU: " << '\t' << updatestat.topRows(1) << '\n' << "CT IN HO: " << '\t' << updatestat.block(7, 0, 1, 6) << '\n' << "INSERT: " << '\t' << updatestat.block(1, 0, 1, 6) << '\n' << "REMOVE: " << '\t' << updatestat.block(2, 0, 1, 6) << '\n' << "SWAP:   " << '\t' << updatestat.block(3, 0, 1, 6) << '\n'<< "SWAPoocc: " << '\t' << updatestat.block(4, 0, 1, 6) << '\n'<< "SWAPoc: " << '\t' << updatestat.block(5, 0, 1, 6) << '\n'<< "SWAPco: " << '\t' << updatestat.block(6, 0, 1, 6) << '\n' << "DQ:      " << '\t' << updatestat.block(8, 0, 1, 6)<< '\n' <<std::endl;
}

void DiagMC::orderstats() {
  MatrixXi orderprint(orderstat.size(),2);
  orderprint << VectorXi::LinSpaced(orderstat.size(), 0, orderstat.size()-1), orderstat;

  std::cout << "Order Statistics" << '\n';
  std::cout<< orderprint  << '\n' << std::endl;
}

void DiagMC::timestats(const VectorXd & timestat) {
  std::cout << "Time Statistics [%]" << std:: endl;
  std::cout << "CHANGE TAU:" << '\t' << timestat(0) << '\n' << "CT HO:     " << '\t' << timestat(1) << '\n' <<  "INSERT: " << '\t' << timestat(2) << '\n' << "REMOVE: " << '\t' << timestat(3) << '\n' <<  "SWAP:     " << '\t' << timestat(4) << '\n' << "DQ:      " << '\t' << timestat(5) << std::endl;
}
  
const MatrixXi & DiagMC::get_Data() {
  MatrixXd  output(taubin, 5);
	
  int CG0p = 0;
  for (int i=0; i< taubin; i++) {
	CG0p += Data(i, 1);
	output(i, 0) = (i+0.5)*taumax/taubin;
  }
  output.rightCols(4)=Data.cast<double>()*(G0p/CG0p);
	
  return output;	
}
