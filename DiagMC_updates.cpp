#include "DiagMC.h"

//Propagators
#ifdef FP
double DiagMC::Vq2(const std::array< double, 3> & q) {
  double pref =  alpha*2.*sqrt(2.)*M_PI;
  return pref/vsq(q);  
}

double DiagMC::disprel(const std::array< double, 3> & q){
  return wp;
}
#endif

#ifdef BEC
double DiagMC::Vq2(const std::array< double, 3> & q) {
  return alpha *M_PI * (1.+ 1./relm)*(1. + 1./relm)* sqrt(vsq(q)/(vsq(q)+2.));   
}

double DiagMC::disprel(const std::array< double, 3> & q){
  return sqrt(vsq(q) *(1.+(vsq(q)/2.)));
}
#endif

double DiagMC::G0el(const std::array< double, 3 > & p, const double & tfin, const double & tinit) {
  return expfun(-((vsq(p)/2./mass) - mu)*(tfin-tinit));
}

double DiagMC::Dph(const std::array< double, 3> & q, const double & tfin, const double & tinit) {
  return expfun(-disprel(q)*(tfin-tinit));
}

double DiagMC::irwsinV(const std::array< double, 3 > & p, const std::array< double, 3 > & q, const double & tfin, const double & tinit){
  double pterm = (vsq(p-q)-vsq(p))/2./mass;
  return expfun(-(pterm + disprel(q))*(tfin-tinit));
}

double DiagMC::tau2pref(const std::array< double, 3 > & p, const std::array< double, 3 > & q) {
  double pterm = (vsq(p-q)-vsq(p))/2./mass;
  return (pterm + disprel(q));
}

 
int DiagMC::change_tau() {
  lu = "change_tau()";
  updatestat(0,0) +=1;						//attempted
  double ntau = diag.get_tinit(2*diag.get_order()) - log(drnd())/E;		//new tau
  double otau = diag.get_tau();
 
  if (ntau < taumap.taumax) {
	if (diag.set_tau(ntau) == 0) {
	  updatestat(0,1) +=1;						//possible
	  updatestat(0,3) +=1;						//accepted
	  global_weight *= G0el(diag.get_p(0), ntau, 0.)/ G0el(diag.get_p(0), otau,0.);
	  return 0;
	}
  }
  return -1;
}

int DiagMC::ct(){
  lu = "ct()";
  updatestat(7,0) +=1;		//attempted

  //---------------Proposing and fast Rejections
  if (diag.get_order() == 0) {return -1;}
  diag.random_arc();
  if (diag.pr_arc == 0){return -1;}
  
  //check Phonon arch
  int open = (diag.get_link(diag.pr_arc) > diag.pr_arc)? 1 :0;
  
  //Getting momentas
  std::array<double, 3> pleft = diag.get_p(diag.pr_arc-1);
  std::array<double, 3> pright = diag.get_p(diag.pr_arc);
  std::array<double, 3> q;
  if (open) {
	q = diag.get_p(diag.pr_arc -1) - diag.get_p(diag.pr_arc);
  } else {
	q = diag.get_p(diag.pr_arc) - diag.get_p(diag.pr_arc - 1);
  }
  
  //Energyfactor for Sampling
  double Efac = (vsq(pleft) - vsq(pright))/2./mass;
  if (open) {Efac -= disprel(q);}
  else {Efac += disprel(q);}
    
  //Getting tau limits
  diag.pr_tauin = diag.get_tinit(diag.pr_arc - 1); //upper limit
  diag.pr_taufin = diag.get_tfin(diag.pr_arc);	//lower limit
  diag.pr_tau2.t = diag.get_tinit(diag.pr_arc); //old tau
  
 
  //Direct Sampling
  if (ct_taucut){
	do{
	  //select direction
	  int dir = (drnd() > 0.5)? 1 : 0; //0: to the left, 1:to the right
	  
	  //Calculate tmax and norm
	  double tmax = 0.;
	  double expEtmax = 0.;
	  if (dir) {tmax = diag.pr_taufin-diag.pr_tau2.t;}
	  else {tmax = diag.pr_tauin-diag.pr_tau2.t;}
	  expEtmax = expfun(- Efac *tmax);
		
	  diag.pr_tau1.t = diag.pr_tau2.t;
	  
	  //create random number
	  double rndnumb = 0.;
	  if ((dir && Efac > 0) ||(!dir && Efac<0)){
		assert(expEtmax<1);
		rndnumb = expEtmax + drnd() * (1-expEtmax);
	  }
	  else {
		assert(expEtmax>1);
		rndnumb = 1 + drnd() * (expEtmax-1);
	  }
	  
	  //Sample tau
	  diag.pr_tau1.t -= log(rndnumb)/Efac;
  
	  
	  int idx = (Efac>0)? 1:0;
	  if (diag.pr_tau1.t <  diag.pr_tau2.t) {testhisto(-taumap.bin(fabs(diag.pr_tau1.t - diag.pr_tau2.t))+taumap.taubin-1, idx) +=1;}
	  else {testhisto(taumap.bin(diag.pr_tau1.t - diag.pr_tau2.t)+taumap.taubin, idx) +=1;}
	  
	  if (!((diag.pr_tau1.t < diag.pr_tau2.t && !dir) || (diag.pr_tau1.t > diag.pr_tau2.t && dir))) {
	  std::cout << testhisto.sum() <<std::endl;
	  std::cout << (diag.pr_tau1.t < diag.pr_tau2.t && !dir) <<'\t' << (diag.pr_tau1.t > diag.pr_tau2.t && dir) << '\t' << ((diag.pr_tau1.t < diag.pr_tau2.t && !dir) || (diag.pr_tau1.t > diag.pr_tau2.t && dir)) <<std::endl;
	  std::cout << dir  << '\t' << open << '\t' << Efac << '\t' << expEtmax <<'\t' << rndnumb <<std::endl;
	  std::cout << diag.pr_tauin << '\t' << diag. pr_tau1.t << '\t' << diag.pr_taufin<< '\t' <<  diag.pr_tau1.t - diag.pr_tau2.t <<'\t'<<std::endl;
	  exit(EXIT_FAILURE);
	  }
	} while(diag.pr_tau1.t > diag.pr_taufin || diag.pr_tau1.t < diag.pr_tauin);
	
  } else {
	//Constant Sampling
    diag.pr_tau1.t = diag.pr_tau2.t + 2.*ctcor*(drnd()-0.5);
	if (diag.pr_tau1.t <  diag.pr_tau2.t) {testhisto(-taumap.bin(fabs(diag.pr_tau1.t - diag.pr_tau2.t))+taumap.taubin-1, 0) +=1;}
	else {testhisto(taumap.bin(diag.pr_tau1.t - diag.pr_tau2.t)+taumap.taubin, 0) +=1;}
	if (diag.pr_tau1.t > diag.pr_taufin || diag.pr_tau1.t < diag.pr_tauin){return -1;}
  }  
  
  double weight = expfun(- Efac*(diag.pr_tau1.t - diag.pr_tau2.t)) ;
  
  updatestat(7,1) +=1;		//possible
  if (drnd() < weight) {
    updatestat(7,3) +=1;		//accepted
    diag.ct();
	global_weight*=weight;
  }
  else {
    updatestat(7,2) +=1;		//rejected
  }  
  return 0; 
}


int DiagMC::insert() {
  lu = "insert()";
  updatestat(1,0) +=1;		//attempted
  
//Proposing and fast rejections
  if (diag.get_order() > maxord-1 && maxord != -1) {return -1;}
  
  //choose vertex
  diag.random_arc();										//arc to insert
#ifdef SELFENERGY
  if (diag.get_order() > 0 && (diag.pr_arc == 0 || diag.pr_arc == 2*diag.get_order())) {return -1;} 
#endif	
  
  //choose q
  if (qsigma != 0) {diag.gauss_q(static_cast<double>(qsigma), dqins);}
  else {diag.linear_q(dqins);}
  if (vsq(diag.pr_q) > dqins*dqins) {return -1;}
#ifdef BEC
  //Q Cut off
  if (vsq(diag.pr_q) > qc*qc){return -2;}
#endif
  diag.pr_p = diag.get_p(diag.pr_arc) - diag.pr_q;
	
  //tau limits
  diag.pr_tauin = diag.get_tinit(diag.pr_arc);
  diag.pr_taufin = diag.get_tfin(diag.pr_arc);
  
  //choose tau_1
  diag.pr_tau1.t = (drnd()*(diag.pr_taufin - diag.pr_tauin))+ diag.pr_tauin;			//tau1
  diag.pr_tau1.link = diag.pr_arc+2;
  
  //choose tau_2
  diag.pr_tau2.link = diag.pr_arc+1;
  if (ins_taucut){
	do{
	  //Cut Distribution of Tau at next Tau
	  diag.pr_tau2.t = - log(drnd())/ tau2pref(diag.get_p(diag.pr_arc), diag.pr_q) + diag.pr_tau1.t;
	} while (diag.pr_tau2.t > diag.pr_taufin || diag.pr_tau2.t < diag.pr_tauin);
	
  } else if (ins_tau2lin){
	//Linear choosen tau2
	diag.pr_tau2.t = (drnd()*(diag.pr_taufin - diag.pr_tau1.t))+ diag.pr_tau1.t;
	  testhisto(taumap.bin(diag.pr_tau2.t - diag.pr_tau1.t),0) += 1;
  } else {
	//Fast Rejection instead
	diag.pr_tau2.t = - log(drnd())/ tau2pref(diag.get_p(diag.pr_arc), diag.pr_q) + diag.pr_tau1.t;
	if (diag.pr_tau2.t > diag.pr_taufin || diag.pr_tau2.t < diag.pr_tau1.t) {return -1;}
  }
	
  updatestat(1,1) +=1;		//possible

  // weight part	
  double weight = irwsinV(diag.get_p(diag.pr_arc), diag.pr_q, diag.pr_tau2.t, diag.pr_tau1.t);
  /*weight *= G0el(diag.pr_p, diag.pr_tau2.t, diag.pr_tau1.t);	//G0(new)
  weight /= G0el(diag.get_p(diag.pr_arc), diag.pr_tau2.t, diag.pr_tau1.t);	//G0(old)
  weight *= Dph(diag.pr_q, diag.pr_tau2.t, diag.pr_tau1.t);	//Phonon(new) */
  weight *= Vq2(diag.pr_q);
  double weight_diff = weight;
  weight /= pow((2.*M_PI), 3);
  // a priori part
  
  //order
  double n = static_cast<double>(diag.get_order());
  weight *= fw;  //Fake-Weight
  if (diag.get_order() == 0) {weight /= fwzero;}
  if (diag.get_order() == 1) {weight /= fwone;}

  weight *= Prem/Pins;	
  weight /= (1. + (n+1.)*2.); 	//remove selecting vertex 
  weight *= (1. + n*2.); 		//insert selecting vertex
  weight *= diag.pr_taufin-diag.pr_tauin; //sample first vertex to insert
  
  //select second vertex to insert
  if (ins_tau2lin) {
	weight *= diag.pr_taufin-diag.pr_tau1.t;
  } else{
  weight /= tau2pref(diag.get_p(diag.pr_arc), diag.pr_q);
  weight /= irwsinV(diag.get_p(diag.pr_arc), diag.pr_q, diag.pr_tau2.t, diag.pr_tau1.t);
  if (ins_taucut) {
	weight *= (1 - expfun(-tau2pref(diag.get_p(diag.pr_arc), diag.pr_q)*(diag.pr_taufin-diag.pr_tau1.t)));
  } 
  }

  //select phonon momentum
  //linear
  if (qsigma != 0) {weight *= 4./3.*M_PI*pow(dqins,3);}
  else {weight *= pow(2.*dqins,3);}
  
  if (drnd() < weight) {
    updatestat(1,3) +=1;		//accepted
    diag.insert();
	global_weight *= weight_diff;
  }
  else {
    updatestat(1,2) +=1;		//rejected
  }  
  return 0;    
}

int DiagMC::remove() {
  lu = "remove()";
  updatestat(2,0) += 1;		//attempted
  
//Proposing and Fast Rejection
  if (diag.get_order()==0) {return -1;}
#ifdef SECUMUL
  if (diag.get_order() < minord + 1) {return -1;} 
#endif
  //choose vertex
  diag.random_arc();
  if (diag.pr_arc != diag.get_link(diag.pr_arc+1)) {return -1;}
  if (vsq(diag.get_p(diag.pr_arc-1) - diag.get_p(diag.pr_arc)) > dqins*dqins) {return -1;}

  diag.pr_tauin = diag.get_tinit(diag.pr_arc-1);
  diag.pr_taufin = diag.get_tfin(diag.pr_arc+1);
	
  diag.pr_tau1.t = diag.get_tinit(diag.pr_arc);
  diag.pr_tau1.link  = -1;
  diag.pr_tau2.t  = diag.get_tfin(diag.pr_arc);
  diag.pr_tau2.link = -1;
	
  diag.pr_q.fill(0);		
	
  diag.pr_p = diag.get_p(diag.pr_arc-1);
  
  
//possible
  updatestat(2,1) +=1;		

  // weight part	
  std::array<double, 3> q = diag.get_p(diag.pr_arc-1)-diag.get_p(diag.pr_arc); 
  
  double weight = 1/irwsinV(diag.pr_p, q, diag.pr_tau2.t, diag.pr_tau1.t);
  /*weight *= G0el(diag.pr_p, diag.pr_tau2.t, diag.pr_tau1.t);	//G0(new)
  weight /= G0el(diag.get_p(diag.pr_arc), diag.pr_tau2.t, diag.pr_tau1.t);	//G0(old)
  weight /= Dph(q, diag.pr_tau2.t, diag.pr_tau1.t);				  //Phonon(old)*/
  weight /= Vq2(q);
  double weight_diff = weight;
  weight *= pow((2.*M_PI), 3);

  // a priori part
  
  //order
  double n = static_cast<double>(diag.get_order());
  weight /= fw;   // for fake function
  if (diag.get_order() == 1) {weight *= fwzero;}
  if (diag.get_order() == 2) {weight *= fwone;}
  
  weight *= Pins/Prem;	
  weight /= (1. + (n-1.)*2.); 	//insert selecting vertex 
  weight *= (1. + n*2.); 		//remove selecting vertex
  
  weight /= diag.pr_taufin-diag.pr_tauin; //sample first vertex to insert
  //select second vertex to insert
  if (ins_tau2lin){
	weight /= diag.pr_taufin-diag.pr_tau1.t;
  } else {
	weight *= irwsinV(diag.pr_p, q, diag.pr_tau2.t, diag.pr_tau1.t);
	weight *= tau2pref(diag.pr_p, q);
	if (ins_taucut) {
	  weight /= (1 - expfun(-tau2pref(diag.get_p(diag.pr_arc), diag.pr_q)*(diag.pr_taufin-diag.pr_tau1.t)));
	} 
  }

  //select phonon momentum
  //linear
  if (qsigma != 0) {weight /= 4./3.*M_PI*pow(dqins,3);}
  else {weight /= pow(2.*dqins,3);}
 
  if (drnd() < weight) {
    updatestat(2,3) +=1;		//accepted
    diag.remove();
	global_weight *= weight_diff;
  }
  else {
    updatestat(2,2) +=1;		//rejected
  }
    
  return 0;
}


int DiagMC::swap() {
  lu = "swap()";
  updatestat(3,0) += 1; //attempted
  updatestat(4,0) += 1;
  updatestat(5,0) += 1;
  updatestat(6,0) += 1;
  int which = diag.propose_swap();
  if (which == -1) {return -1;} //not possible
  
  double weight = G0el(diag.pr_p, diag.pr_taufin, diag.pr_tauin); // G0(new)
  weight /= G0el(diag.get_p(diag.pr_arc), diag.pr_taufin, diag.pr_tauin);  //G0(old)
  
  // both vertices are open or closed
  if (which == 0) {
	updatestat(4,1) +=1; //possible
	if (diag.get_link(diag.pr_arc) > diag.pr_arc) { // both opening
	  weight /= Dph(diag.get_p(diag.pr_arc-1)-diag.get_p(diag.pr_arc),diag.pr_taufin, diag.pr_tauin);
	  weight *= Dph(diag.get_p(diag.pr_arc)-diag.get_p(diag.pr_arc+1),diag.pr_taufin, diag.pr_tauin);
	} else { //both closing
	  weight *= Dph(diag.get_p(diag.pr_arc)-diag.get_p(diag.pr_arc-1),diag.pr_taufin, diag.pr_tauin);
	  weight /= Dph(diag.get_p(diag.pr_arc+1)-diag.get_p(diag.pr_arc),diag.pr_taufin, diag.pr_tauin);
	}	  
	if (drnd() < weight) {
	  updatestat(4,3) +=1; //accepted
	  diag.swap();
	  global_weight *= weight;
	}
	else {
	  updatestat(4,2)+= 1; 	// rejected
	}
  }
  
  // first vertex is start point, second one is end point 
  if (which == -2) {
	updatestat(5,1) +=1; //possible
	
	// 2 phonons are destroyed between tauin and taufin
	weight /= Dph(diag.get_p(diag.pr_arc-1)-diag.get_p(diag.pr_arc),diag.pr_taufin, diag.pr_tauin); // D(q1)
	weight /= Dph(diag.get_p(diag.pr_arc+1)-diag.get_p(diag.pr_arc),diag.pr_taufin, diag.pr_tauin);	// D(q2)
	
	if (drnd() < weight) {
	  updatestat(5,3) +=1; //accepted
	  diag.swap();
	  global_weight *= weight;
	}
	else {
	  updatestat(5,2)+= 1; 	// rejected
	}
  }
  
   // first vertex is end point, second one is start point 
  if (which == 2) {
	updatestat(6,1) +=1; //possible
	
	// 2 phonons are created between tauin and taufin
	weight *= Dph(diag.get_p(diag.pr_arc)-diag.get_p(diag.pr_arc-1), diag.pr_taufin, diag.pr_tauin);
	weight *= Dph(diag.get_p(diag.pr_arc)-diag.get_p(diag.pr_arc+1), diag.pr_taufin, diag.pr_tauin);
	
	if (drnd() < weight) {
	  updatestat(6,3) +=1; //accepted
	  diag.swap();
	  global_weight *= weight;
	}
	else {
	  updatestat(6,2)+= 1; 	// rejected
	}
  } 
  
  
  for (int i =1 ; i<4 ; i++){
	updatestat(3,i) = updatestat(4,i) + updatestat(5,i) + updatestat(6,i);
  }

  return 0;
}

int DiagMC::dq() {

  lu = "dq()";
  updatestat(8,0) +=1;		//attempted
  if(diag.propose_dq(qcor)!=0) {return -1;}
#ifdef BEC
  //Q Cut off
  if (vsq(diag.get_p(diag.pr_arc-1)-diag.get_p(diag.pr_arc)+diag.pr_q) > qc*qc){return -2;}
#endif
  updatestat(8,1) +=1;		//possible
  
  double weight = 1.;
  double arg = 0.;
  for (int i=0 ; i < (diag.get_link(diag.pr_arc) - diag.pr_arc); i++) {
	arg -= vsq(diag.get_p(diag.pr_arc + i) - diag.pr_q)/2./mass * (diag.get_tfin(diag.pr_arc + i), diag.get_tinit(diag.pr_arc + i)); //new
	arg += vsq(diag.get_p(diag.pr_arc + i))/2./mass *(diag.get_tfin(diag.pr_arc + i), diag.get_tinit(diag.pr_arc + i));
    //weight *= G0el(diag.get_p(diag.pr_arc + i) - diag.pr_q, diag.get_tfin(diag.pr_arc + i), diag.get_tinit(diag.pr_arc + i)); //G0el(neu)
    //weight /= G0el(diag.get_p(diag.pr_arc + i), diag.get_tfin(diag.pr_arc + i), diag.get_tinit(diag.pr_arc + i));					//G0el(alt)
  }
  std:: array<double,3>  qold = diag.get_p(diag.pr_arc - 1)- diag.get_p(diag.pr_arc); 
  //Phonon term
  arg -= disprel(qold + diag.pr_q)*(diag.pr_taufin - diag.pr_tauin);
  arg += disprel(qold)*(diag.pr_taufin - diag.pr_tauin);
  //weight *= Dph(qold + diag.pr_q, diag.pr_taufin, diag.pr_tauin);
  //weight /= Dph(qold, diag.pr_taufin, diag.pr_tauin); 
  weight *= expfun(arg);
  
  //alpha term
  weight *= Vq2(qold + diag.pr_q); 
  weight /= Vq2(qold);
   

  if (drnd() < weight) {
    updatestat(8,3) +=1;		//accepted
    diag.dq();
	global_weight *= weight;
	//diag.printall();
  }
  else {
    updatestat(8,2) +=1;		//rejected
  }  
  return 0;    
}



int DiagMC::inscrossed() {
  lu = "inscrossed()";
  updatestat(9,0) +=1;		//attempted
  
//Proposing and fast rejections
  if (diag.get_order() == 0) {return -1;}
  if (diag.get_order() > maxord-1 && maxord != -1) {return -1;}  
  diag.random_arc();										//arc to insert
  if (diag.pr_arc == 0) {return -1;} //We need a vertex in the middle
  
  // linear chosen q
  if (qsigma != 0) {diag.gauss_q(static_cast<double>(qsigma), dqins);}
  else {diag.linear_q(dqins);}
  if (vsq(diag.pr_q) > dqins*dqins) {return -1;}
#ifdef BEC
  //Q Cut off
  if (vsq(diag.pr_q) > qc*qc){return -2;}
#endif  
  
  std::array<double, 3>  pleft = diag.get_p(diag.pr_arc-1);
  std::array<double, 3>  pright = diag.get_p(diag.pr_arc);
  
  //tau
  diag.pr_tauin = diag.get_tinit(diag.pr_arc);
  diag.pr_taufin = diag.get_tfin(diag.pr_arc);
  
  //tau1
  diag.pr_tau1.t = diag.pr_tauin + log(drnd())/tau2pref(pleft, diag.pr_q);
  diag.pr_tau1.link = diag.pr_arc+2;

  if (diag.pr_tau1.t < diag.get_tinit(diag.pr_arc-1)) {return -1;} 
  if (diag.pr_tau1.t > diag.pr_tauin) {return -1;} 
  
  //tau2
  diag.pr_tau2.t = diag.pr_tauin - log(drnd())/tau2pref(pright, diag.pr_q);
  diag.pr_tau2.link = diag.pr_arc;
  if (diag.pr_tau2.t > diag.pr_taufin) {return -1;} 
  if (diag.pr_tau2.t < diag.pr_tauin) {return -1;} 

  updatestat(9,1) +=1;		//possible

  // weight part	
  double weight = irwsinV(pleft, diag.pr_q, diag.pr_tauin, diag.pr_tau1.t);
  weight *= irwsinV(pright, diag.pr_q, diag.pr_tau2.t, diag.pr_tauin);
  weight *= Vq2(diag.pr_q);
  double weight_diff = weight;
  weight /= pow((2.*M_PI), 3);
    
  // a priori part
  
  //order
  double n = static_cast<double>(diag.get_order());
  weight *= fw;  //Fake-Weight
  if (diag.get_order() == 0) {weight /= fwzero;}
  if (diag.get_order() == 1) {weight /= fwone;}

  weight *= Prc/Pic;	
  weight /= (1. + (n+1.)*2.); 	//remove selecting vertex 
  weight *= (1. + n*2.); 		//insert selecting vertex
  
  //tau1
  weight /= tau2pref(pleft, diag.pr_q);
  weight /= irwsinV(pleft, diag.pr_q, diag.pr_tauin, diag.pr_tau1.t); 
  
  //tau2
  weight /= tau2pref(pright, diag.pr_q);
  weight /= irwsinV(pright, diag.pr_q, diag.pr_tau2.t, diag.pr_tauin); 

  //select phonon momentum
  //linear
  if (qsigma != 0) {weight *= 4./3.*M_PI*pow(dqins,3);}
  else {weight *= pow(2.*dqins,3);}
  
  if (drnd() < weight) {
    updatestat(9,3) +=1;		//accepted
    diag.inscrossed();
	global_weight *= weight_diff;
  }
  else {
    updatestat(9,2) +=1;		//rejected
  }  
  return 0;    
}


int DiagMC::remcrossed() {
  lu = "remcrossed()";
  updatestat(10,0) += 1;		//attempted
 
//Proposing and fast rejections
  if (diag.get_order() == 0 || diag.get_order() == 1) {return -1;}
#ifdef SECUMUL
  if (diag.get_order() < minord + 1) {return -1;} 
#endif
  diag.random_arc();	//arc to remove
  if (diag.get_link(diag.pr_arc + 1) != diag.pr_arc-1) {return -1;}
  
  std::array<double, 3> q = diag.get_p(diag.pr_arc-2)-diag.get_p(diag.pr_arc-1); 
  if (vsq(q) > dqins*dqins){return -1;} 

  diag.pr_tauin = diag.get_tinit(diag.pr_arc);
  diag.pr_taufin = diag.get_tfin(diag.pr_arc+1);
	
  diag.pr_tau1.t = diag.get_tinit(diag.pr_arc-1);
  diag.pr_tau1.link  = -1;
  diag.pr_tau2.t  = diag.get_tfin(diag.pr_arc);
  diag.pr_tau2.link = -1;
  
  //pleft
  std::array<double,3> pleft = diag.get_p(diag.pr_arc-2);	
  
  //pright
  std::array<double,3> pright = diag.get_p(diag.pr_arc+1);

  //possible
  updatestat(10,1) +=1;		

  // weight part	
  double weight = 1./irwsinV(pleft, q, diag.pr_tauin, diag.pr_tau1.t);
  weight /= irwsinV(pright, q, diag.pr_tau2.t, diag.pr_tauin);
  weight /= Vq2(q);
  double weight_diff = weight;
  weight *= pow((2.*M_PI), 3);

  // a priori part
  
  //order
  double n = static_cast<double>(diag.get_order());
  weight /= fw;   // for fake function
  if (diag.get_order() == 1) {weight *= fwzero;}
  if (diag.get_order() == 2) {weight *= fwone;}
  
  weight *= Pic/Prc;	
  weight /= (1. + (n-1.)*2.); 	//insert selecting vertex 
  weight *= (1. + n*2.); 		//remove selecting vertex
  
  //tau1
  weight *= tau2pref(pleft, q);
  weight *= irwsinV(pleft, q, diag.pr_tauin, diag.pr_tau1.t); 
  
  //tau2
  weight *= tau2pref(pright, q);
  weight *= irwsinV(pright, q, diag.pr_tau2.t, diag.pr_tauin); 
  
  //select phonon momentum
  //linear
  if (qsigma != 0) {weight /= 4./3.*M_PI*pow(dqins,3);}
  else {weight /= pow(2.*dqins,3);}
 
  if (drnd() < weight) {
    updatestat(10,3) +=1;		//accepted
    diag.remcrossed();
	global_weight *= weight_diff;
  }
  else {
    updatestat(10,2) +=1;		//rejected
  }
  return 0;
}

