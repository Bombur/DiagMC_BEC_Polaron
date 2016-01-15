#include "DiagMC.h"
   
void DiagMC::test() {
  try{
	
	//SECUMUL
#ifdef SECUMUL
	if(ordstsz < 1) {throw ossoor();}
#endif
	//Diagram
	diag.test(qc); 
	if (diag.get_order() > maxord && maxord != -1) {throw maxoerr();}
	
	//tau
	if (diag.get_tau() > taumap.taumax || diag.get_tau() < 0) {throw oor_tau();}
	
	//Weight Check
	//if (Data((int)(drnd()*taumap.taubin), 0) == 0) {throw data_empty();}
	double CG0p = 0;
	for (int i=0; i< taumap.taubin; i++) {
#ifdef SELFENERGY
#ifndef SECUMUL
	  if (Data(i, 0) + 0.0000001<  Data(i,3)+ Data(i, 4)) {
		throw  greenerr();}
#endif	  
#else
	  if (Data(i, 0)< (Data(i, 1) + Data(i, 2) + Data(i, 4)- 0.0000001)) {
		  std::cout << Data.topRows(20)<< std::endl;
		throw  greenerr();}
	  
#endif
	  CG0p += Data(i, 1);	  
	}
	double weight_diff = fabs(weight_calc() - global_weight);
	if (weight_diff > (0.0000001*global_weight)){
	  std::cerr << weight_diff << '\t' << weight_calc() << '\t' << global_weight << std::endl;  
	  if (weight_diff > 1.0e-250) {
		throw weight_check();
	  } else {
		std::cerr << "Warning! Weight Check would fail in double limits!" <<std::endl;
	  }
	}

	//Statistics
	for (int i = 0; i< 4; i++) {
	  if (fabs(updatestat(i,1)- updatestat(i,2) -updatestat(i, 3)) >  0.001) {throw staterr();}
	}
	if ((updatestat(1,3)+updatestat(9,3)) < (updatestat(2,3)+updatestat(10,3))) {throw ins_rem();}
	
  }
  catch (std::exception& e){
	std::cerr << lu <<std::endl;
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
	  //weight /= pow(2*M_PI,3);
	}
  }
  
  return weight;
}
  

	
void DiagMC::printall(){
  
  std::cout << p <<'\t'<< mu <<'\n';
  std::cout << qc <<'\n';
  std::cout << taumap.taumax <<'\n';
  std::cout << taumap.taubin <<'\n';
  std::cout << alpha <<'\n';
  std::cout << relm <<'\n';
  std::cout << E <<'\n';
  std::cout << G0p <<'\n';
  std::cout << maxord <<'\n';
  std::cout  <<'\n';
#ifdef FP
  std::cout << wp <<'\n';
  std::cout  <<'\n';
#endif
  std::cout << ctcor <<'\n';
  std::cout << qcor <<'\n';
  std::cout << dtins <<'\n';
  std::cout << dqins <<'\n';
  std::cout << fw <<'\n';
  std::cout << sigfac <<'\n';
  std::cout <<'\n';
  
  std::cout << minord <<'\n';
  std::cout << ordstsz <<'\n';
  std::cout << ordstep <<'\n';
  std::cout <<'\n';
#ifdef SECUMUL
  std::cout << nnorms <<'\n';
  std::cout << nends <<'\n';
  std::cout <<'\n';
#endif

  std::cout << updatestat <<'\n';
  std::cout << orderstat <<'\n';
  std::cout <<'\n';
  
  std::cout << path <<'\n' <<std::endl;
  std::cout << global_weight <<'\n' <<std::endl;
  std::cout << lu <<'\n' <<std::endl;
  std::cout <<'\n';
  
/*  
   int CG0p = 0;
  for (int i=0; i< taumap.taubin; i++) {
	CG0p += Data(i, 1);
	std::cout  << Data(i,1) <<'\n';
  }
  std::cout<<'\n' << CG0p <<std::endl;
*/
  std::cout << Prem << '\t' << Pins << '\t'<< Pct << '\t'<< Pctho << '\t'<< Psw << '\t'<< Pdq << '\t'<< Piae << '\t'<< Prae<<'\n' <<std::endl;
  std::cout << drnd() <<'\n' <<std::endl;
  std::cout << Meas_its << '\t' << Test_its << '\t'<< Write_its << '\n'; 
  std::cout << RunTime <<'\n' <<std::endl;
  std::cout <<'\n';
  
  std::cout << normmin <<'\n' <<std::endl;
  std::cout << TotRunTime <<'\n' <<std::endl;
  std::cout <<'\n';
  std::cout << std::endl;
  
  diag.printall();
  
}



void DiagMC::printdiag() {
 diag.printall(); 
}
