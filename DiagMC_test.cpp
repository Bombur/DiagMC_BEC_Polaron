#include "DiagMC.h"
   
void DiagMC::test() {
  try{
	//SECUMUL
#ifdef SECUMUL
	if(ordstsz < 2) {throw ossoor();}
#endif
	
	//Diagram
	diag.test(qc); 
	if (diag.get_order() > maxord && maxord != -1) {throw maxoerr();}
	
	//tau
	if (diag.get_tau() > taumax || diag.get_tau() < 0) {throw oor_tau();}
	
	//Weight Check
	//if (Data((int)(drnd()*taubin), 0) == 0) {throw data_empty();}
	double CG0p = 0;
	for (int i=0; i< taubin; i++) {
#ifdef SELFENERGY
#ifndef SECUMUL
	  if (Data(i, 0)<  Data(i,3)+ Data(i, 4)) {throw  greenerr();}
#endif	  
#else
	  if (Data(i, 0)< (Data(i, 1) + Data(i, 2) + Data(i, 4))) {throw  greenerr();}
	  
#endif
	  CG0p += static_cast<double>(Data(i, 1));	  
	}
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
	if (diag.get_link(i) > i) { //opening of an arc
	  std::array<double, 3> q = diag.get_p(i-1)-diag.get_p(i);
	  weight *= Dph(q, diag.get_tinit(diag.get_link(i)), diag.get_tinit(i));	//Dph
	  weight *= Vq2(q);			//
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
  std::cout << relm <<'\n';
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

void DiagMC::printdiag() {
 diag.printall(); 
}
