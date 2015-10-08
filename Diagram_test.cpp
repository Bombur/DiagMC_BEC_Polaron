#include "Diagram.h"

//Errors for Testing

class dnf_diagvec: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Diagram vector sizes do not fit!";
  }
};

class timeserr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Error in Diagram::times!";
  }
};

class phproperr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Error in Diagram::phprop!";
  }
};

class elproperr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Error in Diagram::elprop!";
  }
};

class momerr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Momentum is not conserved!";
  }
};

class propopen: public std::exception
{
  virtual const char* what() const throw()
  {
    return "A propergator is not closed";
  }
};


class swpos: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Number of possible vertices for swap do not fit number of zero loops!";
  }
};


void Diagram::test() {
  try {
	//printall();
	if (times.size() != (2*order)+2) {throw dnf_diagvec();}
	if (phprop.size() != times.size()-1) {throw dnf_diagvec();}
	if (elprop.size() != phprop.size()) {throw dnf_diagvec();}
	
	if (times[0].link != (2*order)+1) {throw timeserr();}
	if (times[(2*order)+1].link != 0) {throw timeserr();}
	
	for (int i = 0; i< 3; i++) {
	  if (fabs(phprop[0][i] - phprop[2*order][i]) > 0.0000001) {throw phproperr();}
	  if (fabs(elprop[0][i] - elprop[2*order][i]) > 0.0000001) {throw elproperr();}
	}
	
	VectorXi tmp=VectorXi::Zero((2*order)+2);
	for (int i=0; i<(2*order)+2; i++) {
	  tmp(times[i].link) +=1;
	  if (times[i].link == i) {throw timeserr();}
	}
	for (int i=0; i<(2*order)+1; i++) {

	  if (times[i].t > times[i+1].t) {throw timeserr();}
	  for (int i2 = 0; i2< 3; i2++) {
		if (fabs((get_p(i).at(i2)+get_q(i).at(i2)) - get_p(0).at(i2)) > 0.0000001)  {throw momerr();}

	  }


	}
	
	for (int i=0; i<(2*order)+2; i++) {
	  if (tmp(i) != 1) {throw propopen();}
	}
  } catch (std::exception& e) {
      std::cerr << e.what() << std::endl;
      exit(EXIT_FAILURE);
  }

}


void Diagram::printall() {
  
  std::cout <<  order << '\t' << drnd()  << std::endl;
  std::cout <<'\n';
  for (vertex i: times){
	  std::cout << i.t <<'\t';
	  std::cout << i.link <<'\t';
	
	std::cout <<'\n';
  }
  std::cout <<'\n';
  for (auto i : phprop){
	  std::cout << i <<'\n';
  }
  std::cout <<'\n';
  for (auto i : elprop){
	  std::cout<< i << '\n';
  }
  std::cout <<'\n';

  std::cout << pr_arc <<'\t' << pr_tauin << '\t' << pr_taufin  << std::endl;
  std::cout <<'\n';
	std::cout << pr_tau1.t <<'\t';
	std::cout << pr_tau1.link <<'\t';
	std::cout <<'\n';
	std::cout <<'\n';
	std::cout << pr_tau2.t <<'\t';
	std::cout << pr_tau2.link <<'\t';
	std::cout <<'\n';
	std::cout <<'\n';
  std::cout << pr_q <<'\n';
  std::cout <<'\n';
  std::cout << pr_p <<'\n';
  std::cout <<std::endl;
  

}
