#include "Diagram.h"

void Diagram::test(const double & qc) {
  try {
	if (times.size() != (2*order)+2) {throw dnf_diagvec();}
	if (elprop.size() != times.size()-1) {throw dnf_diagvec();}
#ifndef NCHECK
	if (elprop.size() != phprop.size()) {throw dnf_diagvec();}
#endif
	
	if (times[0].link != (2*order)+1) {throw timeserr();}
	if (times[(2*order)+1].link != 0) {throw timeserr();}
	
	for (int i = 0; i< 3; i++) {
#ifndef NCHECK
	  if (fabs(phprop[0][i] - phprop[2*order][i]) > 0.0000001) {throw phproperr();}
#endif
	  if (fabs(elprop[0][i] - elprop[2*order][i]) > 0.0000001) {throw elproperr();}
	}
	  
	
	VectorXi tmp=VectorXi::Zero((2*order)+2);
	for (int i=0; i<(2*order)+2; i++) {
	  tmp(times[i].link) +=1;
	  if (times[i].link == i) {throw timeserr();}
	  if (times[i].t <0) {throw timeserr();}
	}
	for (int i=0; i<(2*order)+1; i++) {

	  if (times[i].t >= times[i+1].t) {throw timeserr();}
#ifndef NCHECK
	  if (vsq(get_p(i)+ get_q(i) - get_p(0)) > 0.0000001)  {throw momerr();}
#endif
	  if(times[i].link > i && i>0){
		if (vsq(elprop[i-1]-elprop[i]-elprop[times[i].link]+elprop[times[i].link-1]) > 0.0000001){throw qerr();}
	  }
	}
	
#ifdef BEC	
	//Cut Off Check
	for (auto i = 1; i< (elprop.size()-1); i++) {
	  if (times[i].link > i && vsq(elprop[i-1]- elprop[i]) > qc*qc){throw qoor();}
	}
#endif	
	
	for (int i=0; i<(2*order)+2; i++) {
	  if (tmp(i) != 1) {throw propopen();}
	}
  } catch (std::exception& e) {
      std::cerr << e.what() << std::endl;
      exit(EXIT_FAILURE);
  }

}


void Diagram::printall() {
  
  //std::cout << drnd()  << std::endl;
  std::cout <<  order << '\n';
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

bool Diagram::is_reducible() {
  auto p0 = elprop[0];
  for (int i = 2; i < 2*order; i+=2 ) {
	if (vsq(elprop[i]-p0) < 1e-8) {return true;}
  }
  return false;  
}

int Diagram::arch_num(const int & ver = 0) {
  int lim;
  int num = 0;
  if (ver == 0) {lim = pr_arc;}
  else {lim = ver;} 
  
  if (times[lim].link < lim) {return arch_num(times[lim].link);}
  for (int i = 2; i<  lim+1; i++){
	if (times[i].link > i) { num +=1;}
  }
  return num;
}

void Diagram::capacity_check() {
  if (times.capacity() > befcap) {
	std::cout  << "Capacity of Diagram containers changed from " << befcap << " to " << times.capacity() << "!" <<std::endl;
	befcap = times.capacity();
  }
}
