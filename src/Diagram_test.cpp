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
	if (times.size() != (2*get_order())+2) {throw dnf_vecsize();}
	if (phprop.size() != times.size()-1) {throw dnf_vecsize();}
	if (elprop.size() != phprop.size()) {throw dnf_vecsize();}
	
	if ((int)(times[0][1]+0.5) != (2*get_order())+1) {throw timeserr();}
	if ((int)(times[(2*get_order())+1][1]+0.5) != 0) {throw timeserr();}
	
	for (int i = 0; i< 3; i++) {
	  if (fabs(get_q(0).at(i) - get_q(2*get_order()).at(i)) > 0.0000001) {throw phproperr();}
	  if (fabs(get_p(0).at(i) - get_p(2*get_order()).at(i)) > 0.0000001) {throw elproperr();}
	}
	
	
	int zl = 0; // number of zero loops for possible swap
	  
	VectorXd tmp=VectorXd::Zero((2*get_order())+2);
	for (int i=0; i<(2*get_order())+2; i++) {
	  tmp((int)(times[i][1]+0.5)) +=1;
	  if ((int)(times[i][1]+0.5) == i) {throw timeserr();}
	}
	for (int i=0; i<(2*get_order())+1; i++) {
	  if (times[i] > times[i+1]) {throw timeserr();}
	  for (int i2 = 0; i2< 3; i2++) {
		if (fabs((get_p(i).at(i2)+get_q(i).at(i2)) - get_p(0).at(i2)) > 0.0000001)  {throw momerr();}

	  }
	  if (get_link(i) == (i+1)) {zl += 1;}	  
	}
	
	if ( ((order != 0) && (order != 1)) && (((2*order)-1) != (sw_pos +zl)))  {
	  std::cout << zl <<'\t' << sw_pos << '\t' << ((2*order)-1) << std::endl;
	  throw swpos();
	  
	}
	
	for (int i=0; i<(2*get_order())+2; i++) {
	  if ((int)(tmp(i)+0.5) != 1) {throw propopen();}
	}
  } catch (std::exception& e) {
      std::cerr << e.what() << std::endl;
      exit(EXIT_FAILURE);
  }

}



void Diagram::printall() {
  std::cout <<  order << '\t' << drnd()  << '\t' << sw_pos << std::endl;
  std::cout <<'\n';
  for (int i=0 ;  i<times.size() ; i++){
	for (int i2=0 ;  i2<times[0].size() ; i2++) {
	  std::cout << times[i][i2] <<'\t';
	}
	std::cout <<'\n';
  }
  std::cout <<'\n';
  for (int i=0 ;  i<phprop.size() ; i++){
	for (int i2=0 ;  i2<phprop[0].size() ; i2++) {
	  std::cout << phprop[i][i2] <<'\t';
	}
	std::cout <<'\n';
  }
  std::cout <<'\n';
  for (int i=0 ;  i<elprop.size() ; i++){
	for (int i2=0 ;  i2<elprop[0].size() ; i2++) {
	  std::cout << elprop[i][i2] <<'\t';
	}
	std::cout <<'\n';
  }
  std::cout <<'\n';
  
  std::cout << pr_arc <<'\t' << pr_tauin << '\t' << pr_taufin  << std::endl;
  std::cout <<'\n';
  for (int i=0 ;  i<pr_tau1.size() ; i++){
	std::cout << pr_tau1[i] <<'\t';
	std::cout <<'\n';
  }
  std::cout <<'\n';
  for (int i=0 ;  i<pr_tau2.size() ; i++){
	std::cout << pr_tau2[i] <<'\t';
	std::cout <<'\n';
  }
  std::cout <<'\n';
  for (int i=0 ;  i<pr_q.size() ; i++){
	std::cout << pr_q[i] <<'\t';
	std::cout <<'\n';
  }
  std::cout <<'\n';
  for (int i=0 ;  i<pr_p.size() ; i++){
	std::cout << pr_p[i] <<'\t';
	std::cout <<'\n';
  }
  std::cout <<std::endl;
  

}