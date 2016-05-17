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

double DiagMC::logVq2(const std::array< double, 3 > & q) {
  return std::log(Vq2(q));
}

double DiagMC::logG0el(const std::array< double, 3 > & p, const double & tfin, const double & tinit) {
  return -((vsq(p)/2./mass) - mu)*(tfin-tinit);
}

double DiagMC::G0el(const std::array< double, 3 > & p, const double & tfin, const double & tinit) {
  return expfun(logG0el(p, tfin, tinit));
}

double DiagMC::logDph(const std::array< double, 3 > & q, const double & tfin, const double & tinit) {
  return -disprel(q)*(tfin-tinit);
}

double DiagMC::Dph(const std::array< double, 3> & q, const double & tfin, const double & tinit) {
  return expfun(logDph(q, tfin, tinit));
}

double DiagMC::irwsinV(const std::array< double, 3 > & p, const std::array< double, 3 > & q, const double & tfin, const double & tinit){
  return expfun(- tau2pref(p, q)*(tfin - tinit));
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
 
  if ((ntau > taumap.taumin) && (ntau < taumap.taug0max)) {
	if (diag.set_tau(ntau) == 0) {
	  updatestat(0,1) +=1;						//possible
	  updatestat(0,3) +=1;						//accepted
#ifndef NCHECK
	  global_weight += logG0el(diag.get_p(0), ntau, otau);
#endif
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
  std::array<double, 3> q;
  if (open) {
	q = diag.get_p(diag.pr_arc -1) - diag.get_p(diag.pr_arc);
  } else {
	q = diag.get_p(diag.pr_arc) - diag.get_p(diag.pr_arc - 1);
  }
  
  //Energyfactor for Sampling
  double Efac = (vsq(diag.get_p(diag.pr_arc-1)) - vsq(diag.get_p(diag.pr_arc)))/2./mass;
  if (open) {Efac -= disprel(q);}
  else {Efac += disprel(q);}

  //Getting tau limits
  diag.pr_tauin = diag.get_tinit(diag.pr_arc - 1); //upper limit
  diag.pr_taufin = diag.get_tfin(diag.pr_arc);	//lower limit
  diag.pr_tau1.t = diag.get_tinit(diag.pr_arc); //old tau
  
  //Linear chosen tau
  if (ct_lin){
	diag.pr_tau2.t = diag.pr_tau1.t + 2*ctcor*(drnd()-0.5);
	if (diag.pr_tau2.t >= diag.pr_taufin || diag.pr_tau2.t <= diag.pr_tauin){return -1;}
  }	
  //direct Transformation
  else if (ct_taucut){
	double dir= (drnd() > 0.5)? 1. :-1.;
	//create rndnumb
	double Etmaxl = Efac*(diag.pr_tauin - diag.pr_tau1.t);
	double Etmaxr = Efac*(diag.pr_taufin - diag.pr_tau1.t);
	double rndnumb = 0.;
	
	//Overflow
	if ((Etmaxl < -explim)|| (Etmaxr < -explim)) {rndnumb = drnd()*expfun(explim);}
	else {
	  double expEtmaxl= (Etmaxl > explim)? 0: expfun(-Etmaxl);
	  double expEtmaxr= (Etmaxr > explim)? 0: expfun(-Etmaxr);
	  //normal cases
	  if (Efac > 0) {rndnumb = expEtmaxr + drnd()*(expEtmaxl - expEtmaxr);}
	  else {rndnumb = expEtmaxl + drnd()*(expEtmaxr - expEtmaxl);}
	}	
	//sample tau2
	diag.pr_tau2.t = diag.pr_tau1.t - (dir*log(rndnumb)/Efac);
	if ((diag.pr_tau2.t >= diag.pr_taufin) || (diag.pr_tau2.t <= diag.pr_tauin)){return -1;}
  }	
	//Hybrid Method
  else {
	//select direction
	double dir= (drnd() > 0.5)? 1. :-1.;
	diag.pr_tau2.t = diag.pr_tau1.t - (dir*log(drnd())/Efac);
	if ((diag.pr_tau2.t >= diag.pr_taufin) || (diag.pr_tau2.t <= diag.pr_tauin)){return -1;}
  } 
  
  //For Tau cut
  if (diag.get_order() > 0){
	double deltat = (taumap.taumax + taumap.taumin)/2.;
	if (diag.pr_arc == 1){ deltat =(diag.get_tinit(2*diag.get_order()) - diag.pr_tau2.t);}
	else if (diag.pr_arc == 2*diag.get_order()) {deltat = (diag.pr_tau2.t - diag.get_tinit(1));}
	if ((deltat < taumap.taumin) || (deltat > taumap.taumax)) {return -1;}
  }

  //weight variables
  int accept = 0; //-1:metropolis, 0:rejected explim, 1:accepted explim
  double weight =1.;
  double arg = (- Efac*(diag.pr_tau2.t - diag.pr_tau1.t));
  
  //Underflow
  if (arg < -explim) {accept = 0; overflowstat(7,1)+=1;}
  //Overflow
  else if (arg >explim) {accept =1; overflowstat(7,0)+=1;}
  //Normal Case
  else if (arg > 0.) {accept = 1;}
  else{accept =-1;}
  
 
  updatestat(7,1) +=1;		//possible
  
  if (accept == 0) {
	updatestat(7,2) +=1; //rejected
	overflowstat(7, 2) += 1; //rejected due to overflow
	return -1;
  }
  
  if (accept == 1) {
	updatestat(7,3) +=1;		//accepted
	overflowstat(7,3) +=1;    //accepted due to overflow
	diag.ct();
#ifndef NCHECK
	global_weight += arg;
#endif
	return 0;
  }
  
  if (drnd() < expfun(arg)) {
	assert(accept == -1);
    updatestat(7,3) +=1;		//accepted
    diag.ct();
#ifndef NCHECK
	global_weight += arg;
#endif
  }
  else {
    updatestat(7,2) +=1;		//rejected
  }  
  return 0; 
}


int DiagMC::insert() {
  if (diag.get_order() ==0 && fo_insrem){return fo_insert();}
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
  if (qsigma != 0) {
	diag.gauss_q(static_cast<double>(qsigma), dqins);
  }
  else {
	diag.linear_q(dqins);
	if (vsq(diag.pr_q) > dqins*dqins) {return -1;}
  }
  
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
  
  
  //weight variables
  int accept = 0; //-1:metropolis, 0:rejected explim, 1:accepted explim
  double weight = 1.;  
  
  //choose tau_2
  double Efac = tau2pref(diag.get_p(diag.pr_arc), diag.pr_q); //factor for tau2 distribution
  double tmax = diag.pr_taufin - diag.pr_tau1.t;
  double expEtmax=0.;
  diag.pr_tau2.link = diag.pr_arc+1;

  //Linear choosen tau2  
  if (ins_tau2lin){
	diag.pr_tau2.t = (drnd()*tmax)+ diag.pr_tau1.t;
	
	//Underflow
	if ((Efac*(diag.pr_tau2.t - diag.pr_tau1.t)) > explim) {accept = 0;  overflowstat(1,1)+=1;}
	//Overflow
	else if ((Efac*(diag.pr_tau2.t - diag.pr_tau1.t)) <  -explim) {accept = 1;  overflowstat(1,0)+=1;}
	//normal case
	else {
	  accept = -1;
	  weight *= expfun(-Efac*(diag.pr_tau2.t - diag.pr_tau1.t));
	  weight *= tmax;
	}
  }
  // Tau2 according to distribution
  else if (ins_tau2cut) {
	//with Cut off
	double rndnumb =0.;
	//Underflow
	if (Efac*tmax > explim) {
	  accept = -1;
	  diag.pr_tau2.t = diag.pr_tau1.t - log(drnd())/Efac;
	  weight /= Efac;
	  overflowstat(1,1)+=1;
	} 
	//Overflow
	else if (Efac*tmax < -explim) {
	  accept = 1;
	  expEtmax = expfun(explim);
	  rndnumb = 1 + drnd()*(expEtmax -1);
	  diag.pr_tau2.t = diag.pr_tau1.t - log(rndnumb)/Efac;
	  overflowstat(1,0)+=1;
	}
	//normal case
	else {
	  accept = -1;
	  expEtmax = expfun(-Efac*tmax);
	  if (Efac>0) {rndnumb = expEtmax + drnd()*(1-expEtmax);}
	  else {rndnumb = 1 + drnd()*(expEtmax -1);}
	  diag.pr_tau2.t = diag.pr_tau1.t - log(rndnumb)/Efac;
	  weight *= (1 - expEtmax);
	  weight /= Efac;
	}
  } else {
	//without cut off
	accept = -1;
	diag.pr_tau2.t = diag.pr_tau1.t - log(drnd())/Efac;
	if ((diag.pr_tau2.t >= diag.pr_taufin) || (diag.pr_tau2.t <= diag.pr_tau1.t)) {return -1;}
	weight /= Efac;
  }
  
    
  updatestat(1,1) +=1;		//possible

  if (accept == 0) {
	updatestat(1,2) +=1; //rejected
	overflowstat(1, 2) += 1; //rejected due to overflow
	return -1;
  }
  
  double arg = 0.;
#ifndef NCHECK
   //weight check
  arg = -Efac*(diag.pr_tau2.t - diag.pr_tau1.t);
#endif
  
  if (accept == 1) {
	updatestat(1,3) +=1;		//accepted
	overflowstat(1,3) +=1;    //accepted due to overflow
	diag.insert();
#ifndef NCHECK
	global_weight += (arg + logVq2(diag.pr_q));
#endif
	return 0;
  }
  
  //Rest of weight part
  weight /= pow((2.*M_PI), 3);
  weight *= Vq2(diag.pr_q);
  
  // a priori part
  //order
  double n = static_cast<double>(diag.get_order());
  weight *= fw;  //Fake-Weight
  if (n < fwtab.size()) {weight *= fwtab[n];} 

  weight *= Prem/Pins;	
  weight /= (1. + (n+1.)*2.); 	//remove selecting vertex 
  weight *= (1. + n*2.); 		//insert selecting vertex
  weight *= diag.pr_taufin-diag.pr_tauin; //sample first vertex to insert
  
  //select phonon momentum
  if (qsigma != 0) {weight *= 4./3.*M_PI*pow(dqins,3);}
  else {weight *= pow(2.*dqins,3);}
  
  //Actual metropolis
  if (drnd() < weight) {
	updatestat(1,3) +=1;		//accepted
	diag.insert();
#ifndef NCHECK
	global_weight += (arg + logVq2(diag.pr_q));
#endif
  }
  else {
	updatestat(1,2) +=1;		//rejected
  }  
  return 0;
}



int DiagMC::remove() {
  if (diag.get_order() == 1 && fo_insrem){return fo_remove();}
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
  std::array<double, 3> q = diag.get_p(diag.pr_arc-1)-diag.get_p(diag.pr_arc); 
  if (vsq(q) > dqins*dqins) {return -1;}

  diag.pr_tauin = diag.get_tinit(diag.pr_arc-1);
  diag.pr_taufin = diag.get_tfin(diag.pr_arc+1);
	
  diag.pr_tau1.t = diag.get_tinit(diag.pr_arc);
  diag.pr_tau1.link  = -1;
  diag.pr_tau2.t  = diag.get_tfin(diag.pr_arc);
  diag.pr_tau2.link = -1;
	
  diag.pr_q.fill(0);		
	
  diag.pr_p = diag.get_p(diag.pr_arc-1);
  
  // weight variables
  int accept = 0; //-1:metropolis, 0:rejected explim, 1:accepted explim
  double weight = 1.;
  
  //tau2
  double Efac = tau2pref(diag.pr_p, q); //factor for tau2 distribution
  double tmax = diag.pr_taufin - diag.pr_tau1.t;
  double expEtmax=0.;
  //Linear choosen tau2  
  if (ins_tau2lin){
	//Underflow
	if ((Efac*(diag.pr_tau2.t - diag.pr_tau1.t)) > explim) {accept = 1;  overflowstat(2,1)+=1;}
	//Overflow
	else if ((Efac*(diag.pr_tau2.t - diag.pr_tau1.t)) <  -explim) {accept = 0;  overflowstat(2,0)+=1;}
	//normal case
	else {
	  accept = -1;
	  weight /= expfun(-Efac*(diag.pr_tau2.t - diag.pr_tau1.t));
	  weight /= tmax;
	}

  // Tau2 according to distribution
  } else if (ins_tau2cut) {
	//with Cut off
	//Underflow
	if (Efac*tmax > explim) {
	  accept = -1;
	  weight *= Efac;
	  overflowstat(2,1)+=1;
	} 
	//Overflow
	else if (Efac*tmax < -explim) {accept = 0; overflowstat(1,0)+=1;}
	//normal case
	else {
	  accept = -1;
	  expEtmax = expfun(-Efac*tmax);
	  weight /= (1 - expEtmax);
	  weight *= Efac;
	}
  
  //without cut off
  } else {accept = -1;	weight *= Efac;}
  

  //possible
  updatestat(2,1) +=1;	
 
  //Acceptance and Rejection because of the Overflow
  if (accept == 0) {
	updatestat(2,2) +=1; //rejected
	overflowstat(2, 2) += 1; //rejected due to overflow
	return -1;
  }
  
  double arg = 0.;
#ifndef NCHECK
   //weight check
  arg = - tau2pref(diag.pr_p, q) *(diag.pr_tau2.t - diag.pr_tau1.t);
#endif
  
  if (accept == 1) {
	updatestat(2,3) +=1;		//accepted
	overflowstat(2,3) +=1;    //accepted due to overflow
	diag.remove();
#ifndef NCHECK
	global_weight -= (arg + logVq2(q));
#endif
	return 0;
  }
  
  //From here on we have the normal case
  //Rest of the weight part
  weight *= pow((2.*M_PI), 3);
  weight /= Vq2(q);
  
  // a priori part
  //order
  double n = static_cast<double>(diag.get_order());
  weight /= fw;   // for fake function
  if ((n-1) < fwtab.size()) {weight /= fwtab[n-1];} 
  
  weight *= Pins/Prem;	
  weight /= (1. + (n-1.)*2.); 	//insert selecting vertex 
  weight *= (1. + n*2.); 		//remove selecting vertex
  
  //select phonon momentum
  //linear
  if (qsigma != 0) {weight /= (4./3.*M_PI*pow(dqins,3));}
  else {weight /= pow(2.*dqins,3);}
  
  //tau1
  weight /= diag.pr_taufin-diag.pr_tauin;
  
  if (drnd() < weight) {
	assert(accept == -1);
	updatestat(2,3) +=1;		//accepted
	diag.remove();
#ifndef NCHECK
	global_weight -= (arg+logVq2(q));
#endif
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
  // 0: oc, 1: co, 2: oo, 3: cc
  int which = diag.propose_swap();
  if (which == -1) {return -1;} //not possible
  int idx = (which !=3) ? 4 + which: 6;
  updatestat(3,1) +=1; //possible
  updatestat(idx,1) +=1; //possible
  
  double Efac= (vsq(diag.pr_p) - vsq(diag.get_p(diag.pr_arc)))/2./mass;
  std::array<double,3> qleft;
  std::array<double,3> qright;
  switch(which){
	//oc
	case 0: {qleft = diag.get_p(diag.pr_arc-1) - diag.get_p(diag.pr_arc);
			qright = diag.get_p(diag.pr_arc+1) - diag.get_p(diag.pr_arc);
			Efac -= disprel(qleft);
			Efac -= disprel(qright); break;}
	//co
	case 1: {qleft = diag.get_p(diag.pr_arc) - diag.get_p(diag.pr_arc-1);
			qright = diag.get_p(diag.pr_arc) - diag.get_p(diag.pr_arc+1);
			Efac += disprel(qleft);
			Efac += disprel(qright); break;}
	//oo
	case 2: {qleft = diag.get_p(diag.pr_arc-1) - diag.get_p(diag.pr_arc);
			qright = diag.get_p(diag.pr_arc) - diag.get_p(diag.pr_arc+1);
			Efac -= disprel(qleft);
			Efac += disprel(qright); break;}
	//cc
	case 3: {qleft = diag.get_p(diag.pr_arc) - diag.get_p(diag.pr_arc-1);
			qright = diag.get_p(diag.pr_arc+1) - diag.get_p(diag.pr_arc);
			Efac += disprel(qleft);
			Efac -= disprel(qright); break;}
  }
  
  double arg = Efac*(diag.pr_taufin - diag.pr_tauin);
  
  //Underflow
  if (arg > explim){overflowstat(idx, 1)+=1;}
  //Overflow
  else if (arg < -explim){overflowstat(idx, 0)+=1;}
  //normal case
 
  if (log(drnd()) < (- arg)) {
	updatestat(3,3) +=1;
    updatestat(idx,3) +=1;		//accepted
    diag.swap();
#ifndef NCHECK
	global_weight -= arg;
#endif
  }
  else {
	updatestat(3,2) +=1;
    updatestat(idx,2) +=1;		//rejected
  }  
  return 0; 
}

int DiagMC::dq() {
  lu = "dq()";
  updatestat(8,0) +=1;		//attempted
  if (diag.get_order() == 0) {return -1;}
  diag.random_arc(); 
  if (diag.pr_arc == 0) {return -1;}
  if (diag.get_link(diag.pr_arc) < diag.pr_arc) {diag.pr_arc = diag.get_link(diag.pr_arc);}
	
  diag.pr_tauin = diag.get_tinit(diag.pr_arc);			
  diag.pr_taufin = diag.get_tinit(diag.get_link(diag.pr_arc));   //get_tinit(get_link(pr_arc));
  
  diag.linear_q(qcor);
#ifdef BEC
  //Q Cut off
  if (vsq(diag.get_p(diag.pr_arc-1)-diag.get_p(diag.pr_arc)+diag.pr_q) > qc*qc){return -2;}
#endif
  updatestat(8,1) +=1;		//possible

  double weight = 1.;
  double arg =0.;
  for (int i=0 ; i < (diag.get_link(diag.pr_arc) - diag.pr_arc); i++) {
    arg -= (vsq(diag.get_p(diag.pr_arc + i) - diag.pr_q)/2./mass) * (diag.get_tfin(diag.pr_arc + i)- diag.get_tinit(diag.pr_arc + i)); //new
    arg += (vsq(diag.get_p(diag.pr_arc + i))/2./mass) *(diag.get_tfin(diag.pr_arc + i) - diag.get_tinit(diag.pr_arc + i)); //old
  }
  std:: array<double,3>  qold = diag.get_p(diag.pr_arc - 1)- diag.get_p(diag.pr_arc); 
  //Phonon term
  arg -= disprel(qold + diag.pr_q)*(diag.pr_taufin - diag.pr_tauin);
  arg += disprel(qold)*(diag.pr_taufin - diag.pr_tauin);
  weight *= expfun(arg);
  //alpha term
  weight *= Vq2(qold + diag.pr_q); 
  weight /= Vq2(qold);
   

  if (drnd() < weight) {
    updatestat(8,3) +=1;		//accepted
    diag.dq();
#ifndef NCHECK
	global_weight += (arg + logVq2(qold + diag.pr_q) - logVq2(qold));
#endif
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
  if (qsigma != 0) {
	diag.gauss_q(static_cast<double>(qsigma), dqins);
  }
  else {
	diag.linear_q(dqins);
	if (vsq(diag.pr_q) > dqins*dqins) {return -1;}
  }
#ifdef BEC
  //Q Cut off
  if (vsq(diag.pr_q) > qc*qc){return -2;}
#endif  

  //tau borders and initial ta
  diag.pr_tauin = diag.get_tinit(diag.pr_arc);
  double tmaxl = diag.pr_tauin - diag.get_tinit(diag.pr_arc -1);
  double tmaxr = diag.get_tfin(diag.pr_arc) - diag.pr_tauin;
  
  //tau links
  diag.pr_tau1.link = diag.pr_arc +2;
  diag.pr_tau2.link = diag.pr_arc;
  
  double Eleft = tau2pref(diag.get_p(diag.pr_arc-1), diag.pr_q);
  double Eright = tau2pref(diag.get_p(diag.pr_arc), diag.pr_q);
  double expEtmaxl = 0.;
  double expEtmaxr = 0.;
  
  //weight variables
  int accept = -1; //-1:metropolis, 0:rejected explim, 1:accepted explim
  double weight = 1.;  
  
  //linear chosen taus
  if (ic_taulin) {
	diag.pr_tau1.t  = diag.pr_tauin - drnd()* tmaxl;
	diag.pr_tau2.t  = diag.pr_tauin + drnd()* tmaxr;
	double Earg = (Eleft*(diag.pr_tauin - diag.pr_tau1.t)) + (Eright*(diag.pr_tau2.t- diag.pr_tauin));
	//Underflow
	if (Earg > explim) {accept = 0; overflowstat(9,1)+=1;}
	//Overflow
	else if (Earg < -explim) {accept = 1; overflowstat(9,0)+=1;}
	//normal case
	else{
	  accept = -1;
	  weight *= expfun(-Earg);
	  weight *= (tmaxl*tmaxr);
	}
  }
  //taus according to distribution
  //with taucutoff
  else if (ic_taucut){
	//with Cut off
	double rndnumb =0.;
	
	//tau1 
	//Underflow
	if (Eleft*tmaxl > explim) {
	  accept = -1;
	  diag.pr_tau1.t = diag.pr_tauin + log(drnd())/Eleft;
	  weight /= Eleft;
	  overflowstat(9,1)+=1;
	} 
	//Overflow
	else if (Eleft*tmaxl < -explim) {
	  accept = 1;
	  expEtmaxl = expfun(explim);
	  rndnumb = 1. + drnd()*(expEtmaxl -1.);
	  diag.pr_tau1.t = diag.pr_tauin + log(rndnumb)/Eleft;
	  overflowstat(9,0)+=1;
	}
	//normal case
	else {
	  accept = -1;
	  expEtmaxl = expfun(-Eleft*tmaxl);
	  if (Eleft>0) {rndnumb = expEtmaxl + drnd()*(1-expEtmaxl);}
	  else {rndnumb = 1. + drnd()*(expEtmaxl -1.);}
	  diag.pr_tau1.t = diag.pr_tauin + log(rndnumb)/Eleft;
	  weight *= (1. - expEtmaxl);
	  weight /= Eleft;
	}
	
	//tau2
	//Underflow
	if (Eright*tmaxr > explim) {
	  accept = (accept == 1)? 1: -1;
	  diag.pr_tau2.t = diag.pr_tauin - log(drnd())/Eright;
	  weight /= Eright;
	  overflowstat(9,1)+=1;
	} 
	//Overflow
	else if (Eright*tmaxr < -explim) {
	  accept = 1;
	  expEtmaxr = expfun(explim);
	  rndnumb = 1. + drnd()*(expEtmaxr -1.);
	  diag.pr_tau2.t = diag.pr_tauin - log(rndnumb)/Eright;
	  overflowstat(9,0)+=1;
	}
	//normal case
	else {
	  accept = (accept == 1)? 1: -1;
	  expEtmaxr = expfun(-Eright*tmaxr);
	  if (Eright>0) {rndnumb = expEtmaxr + drnd()*(1.-expEtmaxr);}
	  else {rndnumb = 1. + drnd()*(expEtmaxr -1.);}
	  diag.pr_tau2.t = diag.pr_tauin - log(rndnumb)/Eright;
	  weight *= (1 - expEtmaxr);
	  weight /= Eright;
	}
  }
  //without cut off
  else{
	accept = -1;
	diag.pr_tau1.t = diag.pr_tauin + log(drnd())/Eleft;
	if ((diag.pr_tau1.t >= diag.pr_tauin) || (diag.pr_tau1.t <= (diag.pr_tauin - tmaxl))) {return -1;}
	diag.pr_tau2.t = diag.pr_tauin - log(drnd())/Eright;
	if ((diag.pr_tau2.t <= diag.pr_tauin) || (diag.pr_tau2.t >= (diag.pr_tauin + tmaxr))) {return -1;}
	weight /= Eleft;
	weight /= Eright;
  }
  
  //For Tau cut
  double deltat = (taumap.taumax + taumap.taumin)/2.;
  if (diag.pr_arc == 1){ deltat =(diag.get_tinit(2*diag.get_order()) - diag.pr_tau1.t);}
  else if (diag.pr_arc == 2*diag.get_order()) {deltat = (diag.pr_tau2.t - diag.get_tinit(1));}
  if (deltat > taumap.taumax) {return -1;}

	
  updatestat(9,1) +=1;		//possible
  
  if (accept == 0) {
	updatestat(9,2) +=1; //rejected
	overflowstat(9,2) += 1; //rejected due to overflow
	return -1;
  }
  
  
#ifndef NCHECK
   //weight check
  double arg = (-(Eleft*(diag.pr_tauin - diag.pr_tau1.t))- (Eright*(diag.pr_tau2.t-diag.pr_tauin)));
#endif
  
  if (accept == 1) {
	updatestat(9,3) +=1;		//accepted
	overflowstat(9,3) +=1;    //accepted due to overflow
	diag.inscrossed();
#ifndef NCHECK
	global_weight += (arg + logVq2(diag.pr_q));
#endif
	return 0;
  }
  
  //Rest of weight part
  weight /= pow((2.*M_PI), 3);
  weight *= Vq2(diag.pr_q);
  
  // a priori part
  //order
  double n = static_cast<double>(diag.get_order());
  weight *= fw;  //Fake-Weight
  if (n < fwtab.size()) {weight *= fwtab[n];} 

  weight *= Pic/Pic;	
  weight /= (1. + (n+1.)*2.); 	//remove selecting vertex 
  weight *= (1. + n*2.); 		//insert selecting vertex
  
  //select phonon momentum
  if (qsigma != 0) {weight *= 4./3.*M_PI*pow(dqins,3);}
  else {weight *= pow(2.*dqins,3);}
  
  //Actual metropolis
  if (drnd() < weight) {
	assert(accept == -1);
	updatestat(9,3) +=1;		//accepted
	diag.inscrossed();
#ifndef NCHECK
	global_weight += (arg + logVq2(diag.pr_q));
#endif
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
  if (diag.get_order() < 2) {return -1;}
#ifdef SECUMUL
  if (diag.get_order() < minord + 1) {return -1;} 
#endif
  diag.random_arc();	//arc to remove
  if (diag.pr_arc == 0) {return -1;}
  if (diag.get_link(diag.pr_arc + 1) != (diag.pr_arc-1)) {return -1;}
  
  std::array<double, 3> q = (diag.get_p(diag.pr_arc-2)-diag.get_p(diag.pr_arc-1)); 
  if (vsq(q) > dqins*dqins){return -1;} 

  //tau1 and tau2
  diag.pr_tau1.t = diag.get_tinit(diag.pr_arc-1);
  diag.pr_tau1.link  = -1;
  diag.pr_tau2.t  = diag.get_tfin(diag.pr_arc);
  diag.pr_tau2.link = -1;
  
  //tau borders and initial tau
  diag.pr_tauin = diag.get_tinit(diag.pr_arc);
  double deltat = (taumap.taumax + taumap.taumin)/2.;
  if (diag.pr_arc == 2){deltat =(diag.get_tinit(2*diag.get_order()) - diag.pr_tauin);}
  else if (diag.pr_arc == (2*diag.get_order()-1)) {deltat = (diag.pr_tauin - diag.get_tinit(1));}
  if (deltat < taumap.taumin) {return -1;}
  
  double tmaxl = diag.pr_tauin - diag.get_tinit(diag.pr_arc -2);
  double tmaxr = diag.get_tfin(diag.pr_arc+1) - diag.pr_tauin;
  
  double Eleft = tau2pref(diag.get_p(diag.pr_arc-2), q);
  double Eright = tau2pref(diag.get_p(diag.pr_arc+1), q);

  //weight variables
  int accept = -1; //-1:metropolis, 0:rejected explim, 1:accepted explim
  double weight = 1.;  
  
  //linear chosen taus
  if (ic_taulin){
	double Earg = (Eleft*(diag.pr_tauin - diag.pr_tau1.t)) + (Eright*(diag.pr_tau2.t- diag.pr_tauin));
	//Underflow
	if (Earg > explim){accept=1; overflowstat(10,1)+=1;}
	//Overflow
	if (Earg < -explim) {accept = 0; overflowstat(10,0) +=1;}
	//normal case
	else {
	  accept = -1;
	  weight /= expfun(-Earg);
	  weight /= tmaxl*tmaxr;
	}
  }
  //taus according to weights
  if (ic_taucut){
	//with tau cut offs
	
	//tau1 
	//Underflow
	if (Eleft*tmaxl > explim) {
	  accept = -1;
	  weight *= Eleft;
	  overflowstat(10,1)+=1;
	} 
	//Overflow
	else if (Eleft*tmaxl < -explim) {accept = 0;  overflowstat(10,0)+=1;}
	//normal case
	else {
	  accept = -1;
	  weight /= (1. - expfun(-Eleft*tmaxl));
	  weight *= Eleft;
	}
	
	//tau2
	if (accept !=0){
	  //Underflow
	  if (Eright*tmaxr > explim) {
		accept = (accept == 1)? 1: -1;
		weight *= Eright;
		overflowstat(10,1)+=1;
	  } 
	  //Overflow
	  else if (Eright*tmaxr < -explim) {accept = 0;  overflowstat(10,0)+=1;}
	  //normal case
	  else {
		accept = (accept == 1)? 1: -1;
		weight /= (1 - expfun(-Eright*tmaxr));
		weight *= Eright;
	  }
	}
  }
  //without tau cut offs
  else{
	accept = -1;
	weight *= Eleft;
	weight *= Eright;
  }

  //possible
  updatestat(10,1) +=1;		

  if (accept == 0) {
	updatestat(10,2) +=1; //rejected
	overflowstat(10,2) += 1; //rejected due to overflow
	return -1;
  }
  
#ifndef NCHECK
   //weight check
  double arg = (-(Eleft*(diag.pr_tauin - diag.pr_tau1.t))- (Eright*(diag.pr_tau2.t-diag.pr_tauin)));
#endif
  
  if (accept == 1) {
	updatestat(10,3) +=1;		//accepted
	overflowstat(10,3) +=1;    //accepted due to overflow
	diag.remcrossed();
#ifndef NCHECK
	global_weight -= (arg +logVq2(q)) ;
#endif
	return 0;
  }
  
  //Rest of weight part
  weight *= pow((2.*M_PI), 3);
  weight /= Vq2(q);
  
  // a priori part
  //order
  double n = static_cast<double>(diag.get_order());
  weight /= fw;   // for fake function
  if ((n-1) < fwtab.size()) {weight /= fwtab[n-1];} 
  
  weight *= Pic/Prc;	
  weight /= (1. + (n-1.)*2.); 	//insert selecting vertex 
  weight *= (1. + n*2.); 		//remove selecting vertex
  
  //select phonon momentum
  //linear
  if (qsigma != 0) {weight /= (4./3.*M_PI*pow(dqins,3));}
  else {weight /= pow(2.*dqins,3);}
  
  //Actual metropolis
  if (drnd() < weight) {
	assert(accept == -1);
	updatestat(10,3) +=1;		//accepted
	diag.remcrossed();
#ifndef NCHECK
	global_weight -= (arg +logVq2(q)) ;
#endif
  }
  else {
	updatestat(9,2) +=1;		//rejected
  }  
  return 0;
}

int DiagMC::fo_insert() {
  lu = "fo_insert()";
  updatestat(11,0) +=1;		//attempted
  
//Proposing and fast rejections
  if (diag.get_order() != 0) {return -1;}
//choose vertex
  diag.pr_arc =0;
 
//choose q
  if (qsigma != 0) {
	diag.gauss_q(static_cast<double>(qsigma), dqins);
  }
  else {
	diag.linear_q(dqins);
	if (vsq(diag.pr_q) > dqins*dqins) {return -1;}
  }
#ifdef BEC
  //Q Cut off
  if (vsq(diag.pr_q) > qc*qc){return -2;}
#endif
  diag.pr_p = diag.get_p(diag.pr_arc) - diag.pr_q;
  	
  //tau limits
  diag.pr_tauin = diag.get_tinit(diag.pr_arc);
  diag.pr_taufin = diag.get_tfin(diag.pr_arc);
  
  //choose tau_1
  diag.pr_tau1.t = (drnd()*(diag.pr_taufin - taumap.taumin - diag.pr_tauin))+ diag.pr_tauin;			//tau1
  diag.pr_tau1.link = diag.pr_arc+2;
  
  //choose tau_2
  diag.pr_tau2.t = (drnd()*(diag.pr_taufin - taumap.taumin - diag.pr_tau1.t))+ diag.pr_tau1.t + taumap.taumin;
  diag.pr_tau2.link = diag.pr_arc+1;
  
  if ((diag.pr_tau2.t - diag.pr_tau1.t) > taumap.taumax) {return -1;}
  
  //weight variables
  int accept = 0; //-1:metropolis, 0:rejected explim, 1:accepted explim
  double weight = 1.; 
  double arg = -tau2pref(diag.get_p(diag.pr_arc), diag.pr_q)*(diag.pr_tau2.t - diag.pr_tau1.t);
  
  //Underflow
  if (arg < -explim) {accept = 0;  overflowstat(11,1)+=1;}
  //Overflow
  else if (arg > explim) {accept = 1;  overflowstat(11,0)+=1;}
  //normal case
  else {accept = -1;
	weight *= expfun(arg);
  }
    
  updatestat(11,1) +=1;		//possible

  if (accept == 0) {
	updatestat(11,2) +=1; //rejected
	overflowstat(11, 2) += 1; //rejected due to overflow
	return -1;
  }
  
  if (accept == 1) {
	updatestat(11,3) +=1;		//accepted
	overflowstat(11,3) +=1;    //accepted due to overflow
	diag.insert();
#ifndef NCHECK
	global_weight += (arg + logVq2(diag.pr_q));
#endif
	return 0;
  }
  
  //Rest of weight part
  weight /= pow((2.*M_PI), 3);
  weight *= Vq2(diag.pr_q);
  
  // a priori part
  //order
  weight *= fw;  //Fake-Weight
  weight *= fwtab[0]; 

  weight *= Prem/Pins;	
  weight *= (diag.pr_taufin-taumap.taumin-diag.pr_tauin); //sample first vertex to insert
  weight *= (diag.pr_taufin-taumap.taumin-diag.pr_tau1.t); //sample first vertex to insert
  
  //select phonon momentum
  if (qsigma != 0) {weight *= 4./3.*M_PI*pow(dqins,3);}
  else {weight *= pow(2.*dqins,3);}
  
  //Actual metropolis
  if (drnd() < weight) {
	updatestat(11,3) +=1;		//accepted
	diag.insert();
#ifndef NCHECK
	global_weight += (arg + logVq2(diag.pr_q));
#endif
  }
  else {
	updatestat(11,2) +=1;		//rejected
  }  
  return 0;
}

int DiagMC::fo_remove() {
  lu = "fo_remove()";
  updatestat(12,0) += 1;		//attempted
  
//Proposing and Fast Rejection
  if (diag.get_order()!=1) {return -1;}

//choose vertex
  diag.pr_arc = 1;
  std::array<double, 3> q = diag.get_p(diag.pr_arc-1)-diag.get_p(diag.pr_arc); 
  if (vsq(q) > dqins*dqins) {return -1;}

  diag.pr_tauin = diag.get_tinit(diag.pr_arc-1);
  diag.pr_taufin = diag.get_tfin(diag.pr_arc+1);
	
  diag.pr_tau1.t = diag.get_tinit(diag.pr_arc);
  diag.pr_tau1.link  = -1;
  diag.pr_tau2.t  = diag.get_tfin(diag.pr_arc);
  diag.pr_tau2.link = -1;

  if ((diag.pr_tau2.t - diag.pr_tau1.t) > taumap.taumax) {return -1;}
  
  diag.pr_q.fill(0);		
	
  diag.pr_p = diag.get_p(diag.pr_arc-1);
  
  // weight variables
  int accept = 0; //-1:metropolis, 0:rejected explim, 1:accepted explim
  double weight = 1.;
  double arg = -tau2pref(diag.pr_p, q)*(diag.pr_tau2.t - diag.pr_tau1.t); //factor for tau2 distribution

  //Underflow
  if (arg < -explim) {accept = 1;  overflowstat(12,1)+=1;}
  //Overflow
  else if (arg > explim) {accept = 0;  overflowstat(12,0)+=1;}
  //normal case
  else {
	  accept = -1;
	  weight /= expfun(arg);
  }

  //possible
  updatestat(12,1) +=1;	
 
  //Acceptance and Rejection because of the Overflow
  if (accept == 0) {
	updatestat(12,2) +=1; //rejected
	overflowstat(12, 2) += 1; //rejected due to overflow
	return -1;
  }
  
  
  if (accept == 1) {
	updatestat(12,3) +=1;		//accepted
	overflowstat(12,3) +=1;    //accepted due to overflow
	diag.remove();
#ifndef NCHECK
	global_weight -= (arg + logVq2(q));
#endif
	return 0;
  }
  
  //From here on we have the normal case
  //Rest of the weight part
  weight *= pow((2.*M_PI), 3);
  weight /= Vq2(q);
  
  // a priori part
  //order
  weight /= fw;   // for fake function
  weight /= fwtab[0];
  
  weight *= Pins/Prem;	
 
  //select phonon momentum
  //linear
  if (qsigma != 0) {weight /= (4./3.*M_PI*pow(dqins,3));}
  else {weight /= pow(2.*dqins,3);}
  
  //tau1
  weight /= (diag.pr_taufin-taumap.taumin-diag.pr_tauin);
  //tau2
  weight /= (diag.pr_taufin-taumap.taumin-diag.pr_tau1.t);
  
  if (drnd() < weight) {
	updatestat(12,3) +=1;		//accepted
	diag.remove();
#ifndef NCHECK
	global_weight -= (arg+logVq2(q));
#endif
  }
  else {
	updatestat(12,2) +=1;		//rejected
  }
  return 0;
}
