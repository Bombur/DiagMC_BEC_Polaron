#include <exception>
#include <stdexcept>
#include "DiagMC.h"
   
//tests
class oor_Gp: public std::exception
{
  virtual const char* what() const throw()
  {
    return "G(p) is out of range!";
  }
};

class oor_tau: public std::exception
{
  virtual const char* what() const throw()
  {
    return "tau is out of range!";
  }
};

class data_empty: public std::exception
{
  virtual const char* what() const throw()
  {
    return "G(p, tau_i) is empty!";
  }
};

class greenerr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Error in the Data of the Green's function!";
  }
};

class weight_check: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Weight Check failed!";
  }
};

class fake_check: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Fake Check failed!";
  }
};

class staterr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Statistics do not match!";
  }
};

class ins_rem: public std::exception
{
  virtual const char* what() const throw()
  {
    return "More Removes than Inserts!";
  }
};

class maxoerr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Maximum Order does not hold!";
  }
};



void DiagMC::test() {
  try{
	//Diagram
	diag.test(); 
	if (diag.get_order() > maxord && maxord != -1) {throw maxoerr();}
	
	//tau
	if (diag.get_tau() > taumax || diag.get_tau() < 0) {throw oor_tau();}
	
	//Weight Check and Fake Check
	//if (Data((int)(drnd()*taubin), 0) == 0) {throw data_empty();}
	double CG0p = 0;
	for (int i=0; i< taubin; i++) {
#ifdef SELFENERGY
	  if (Data(i, 0)<  Data(i, 3)) {throw  greenerr();}
#else
	  if (Data(i, 0)< (Data(i, 1) + Data(i, 2) + Data(i, 3))) {throw  greenerr();}
#endif
	  CG0p += static_cast<double>(Data(i, 1));	  
	}
	/*
	for (int i=0; i< taubin; i++) {
	  double G0deltatau = (exp(-(E*static_cast<double>(i)*taumax/static_cast<double>(taubin)))/E) *(1 - exp(-E*(taumax/static_cast<double>(taubin))));
	  //std::cout << (G0deltatau/G0p * CG0p)<<'\t' << static_cast<double>(Data(i,1)) << '\t' << fabs((G0deltatau/G0p * CG0p) - static_cast<double>(Data(i,1))) << std::endl;
	  if (((G0deltatau/G0p) * CG0p) > 100 && fabs((G0deltatau/G0p * CG0p) - static_cast<double>(Data(i,1))) > 0.5*(G0deltatau/G0p) * CG0p )  {throw fake_check();}
	  if (((G0deltatau/G0p) * CG0p) > 1000 && fabs((G0deltatau/G0p * CG0p) - static_cast<double>(Data(i,1))) > 0.16*(G0deltatau/G0p) * CG0p )  {throw fake_check();}
	  if (((G0deltatau/G0p) * CG0p) > 10000 && fabs((G0deltatau/G0p * CG0p) - static_cast<double>(Data(i,1))) > 0.05*(G0deltatau/G0p) * CG0p )  {throw fake_check();}
	  if (((G0deltatau/G0p) * CG0p) > 100000 && fabs((G0deltatau/G0p * CG0p) - static_cast<double>(Data(i,1))) > 0.016*(G0deltatau/G0p) * CG0p )  {throw fake_check();}
	}
	*/
	double weight_diff = fabs(weight_calc() - global_weight);
	if (weight_diff > (0.0000001*global_weight)) {
	  std::cout << weight_diff << '\t' << weight_calc() << '\t' << global_weight << std::endl;  
	  throw weight_check();
	}
	
	//Statistics
	for (int i = 0; i< 4; i++) {
	  if (fabs(updatestat(i,1)- updatestat(i,2) -updatestat(i, 3)) >  0.001) {throw staterr();}
	}
	if ((updatestat(1,3)+updatestat(9,3)) < (updatestat(2,3)+updatestat(10,3))) {throw ins_rem();}
	
  }
  catch (std::exception& e){
	std::cout << lu <<std::endl;
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

double DiagMC::weight_calc() {
  double weight = G0el(diag.get_p(0), diag.get_tfin(0), 0.);  //G0(p, t1, 0)
  
  for (int i=1; i< 1+2*diag.get_order(); i++) {
	weight *= G0el(diag.get_p(i), diag.get_tfin(i), diag.get_tinit(i));   //G0(pi, ti+1, ti)
	if (diag.get_link(i) > i) { 											//opening of an arc
	  weight *= Dph(diag.get_tinit(diag.get_link(i)), diag.get_tinit(i));	//Dph
	  weight *= alpha / vsq(diag.get_p(i-1) - diag.get_p(i));			//alpha/q^2
	  weight /= pow(2*M_PI,3);
	}
  }
  
  return weight;
}
  

	
void DiagMC::printall(){
  
  std::cout << p <<'\t'<< mu <<'\n';
  std::cout << taumax <<'\n';
  std::cout << taubin <<'\n';
  std::cout << alpha <<'\n';
  std::cout << omegap <<'\n';
  std::cout << E <<'\n';
  std::cout << G0p <<'\n';
  std::cout  <<'\n';
  std::cout << ctcor <<'\n';
  std::cout << qcor <<'\n';
  std::cout << global_weight <<'\n' <<std::endl;
  
   int CG0p = 0;
  for (int i=0; i< taubin; i++) {
	CG0p += Data(i, 1);
	std::cout  << Data(i,1) <<'\n';
  }
  std::cout<<'\n' << CG0p <<std::endl;
  
}
