#include "DiagMC.h"

 
int DiagMC::change_tau() {
  lu = "change_tau()";
  double ntau = - log(drnd())/E;		//new tau
  double otau = diag.get_tau();
  updatestat(0,0) +=1;						//attempted
  
  if (diag.get_order() != 0) {return -1;}
  //if (diag.get_order() != 0 || ntau<diag.get_tinit(2*diag.get_order()-1)) {return -1;}

    
  if (ntau < taumax) {
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
  int which = diag.propose_ct(taumax, ctcor);
  double weight = 0.;
  if(which == 0) {
	weight = G0el(diag.pr_p, diag.pr_taufin, diag.pr_tauin); //G0 (new) 
  } 
   else if (which==1){ // opening arc
	weight = G0el(diag.pr_p, diag.pr_taufin, diag.pr_tauin); //G0(new)
	weight /=G0el(diag.pr_q, diag.pr_taufin, diag.pr_tauin); //G0(old)
	
	std::array<double, 3> q = diag.get_p(diag.pr_arc-1)-diag.get_p(diag.pr_arc);
	weight /= Dph(q, diag.pr_taufin, diag.pr_tauin); // Dph(old)  or new depending on t_fin
  }
  else if (which==2){ // closing arc
	weight = G0el(diag.pr_p, diag.pr_taufin, diag.pr_tauin); //G0(new)
	weight /=G0el(diag.pr_q, diag.pr_taufin, diag.pr_tauin); //G0(old)
	
	std::array<double, 3> q = diag.get_p(diag.pr_arc)-diag.get_p(diag.pr_arc-1);
	weight *= Dph(q, diag.pr_taufin, diag.pr_tauin); // Dph(new)  or old depending on t_fin
  }
  else {return -1;}

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
  if (diag.get_order() > maxord-1 && maxord != -1) {return -1;}
#ifdef FP
  const double qctemp = dqins;
#elif BEC
  const double qctemp = qc;
#endif
  if(diag.propose_insert(qctemp, sigfac)!=0) {return -1;}
  updatestat(1,1) +=1;		//possible

  // weight part	
  double weight = G0el(diag.pr_p, diag.pr_tau2.t, diag.pr_tau1.t);	//G0(new)
  weight /= G0el(diag.get_p(diag.pr_arc), diag.pr_tau2.t, diag.pr_tau1.t);	//G0(old)
  weight *= Dph(diag.pr_q, diag.pr_tau2.t, diag.pr_tau1.t);	//Phonon(new)
  weight *= Vq2(diag.pr_q);
  weight /= pow((2.*M_PI), 3);
  
  double weight_diff = weight;
    
  // a priori part
  
  //order
  double n = static_cast<double>(diag.get_order());
  if (diag.get_order() == 0) {weight /= fw;}   // for fake function


  weight *= Prem/Pins;	
  weight /= (1. + (n+1.)*2.); 	//remove selecting vertex 
  weight *= (1. + n*2.); 		//insert selecting vertex
  
  weight *= diag.pr_taufin-diag.pr_tauin; //sample first vertex to insert
  weight *= diag.pr_taufin-diag.pr_tau1.t; //select second vertex to insert

  //select phonon momentum
#ifdef QGAUSS
  //Gauss
  weight *= pow(M_PI* 2.*qctemp/sigfac *qctemp/sigfac, 3./2.);
  weight /= exp(-vsq(diag.pr_q)/2./(qctemp/sigfac)/(qctemp/sigfac));
#else  
  //linear
  weight *= pow(2*qctemp, 3);
#endif
  
  
  //weight /= pow((diag.pr_tau2.t-diag.pr_tau1.t)/2./M_PI, 3./2.) ;
  //weight /= exp((-vsq(diag.pr_q)/2.) *(diag.pr_tau2.t-diag.pr_tau1.t));
  

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

#ifdef SECUMUL
  if (diag.get_order() < minord + 1) {return -1;} 
#endif

#ifdef FP
  const double qctemp = dqins;
#elif BEC
  const double qctemp = qc;
#endif

  if(diag.propose_remove(qctemp)!=0) {return -1;}
  updatestat(2,1) +=1;		//possible

  // weight part	
  std::array<double, 3> q = diag.get_p(diag.pr_arc-1)-diag.get_p(diag.pr_arc); 
  
  double weight = G0el(diag.pr_p, diag.pr_tau2.t, diag.pr_tau1.t);	//G0(new)
  weight /= G0el(diag.get_p(diag.pr_arc), diag.pr_tau2.t, diag.pr_tau1.t);	//G0(old)
  weight /= Dph(q, diag.pr_tau2.t, diag.pr_tau1.t);				  //Phonon(old)
  weight /= Vq2(q);
  weight *= pow((2.*M_PI), 3);

  double weight_diff = weight;

  // a priori part
  
  //order
  double n = static_cast<double>(diag.get_order());
  if (diag.get_order() == 1) {weight *= fw;}   // for fake function

  weight *= Pins/Prem;	
  weight /= (1. + (n-1.)*2.); 	//insert selecting vertex 
  weight *= (1. + n*2.); 		//remove selecting vertex
  
  weight /= diag.pr_taufin-diag.pr_tauin; //sample first vertex to insert
  weight /= diag.pr_taufin-diag.pr_tau1.t; //select second vertex to insert

  //select phonon momentum
#ifdef QGAUSS
  //Gauss
  weight /= pow(M_PI* 2.*qctemp/sigfac*qctemp/sigfac, 3./2.);
  weight *= exp(-vsq(q)/2./(qctemp/sigfac)/(qctemp/sigfac));
#else
  //linear
  weight /= pow(2*qctemp, 3);
#endif
  //weight *= pow((diag.pr_tau2.t-diag.pr_tau1.t)/2./M_PI, 3./2.) ;
  //weight*= exp((-vsq(diag.get_q(diag.pr_arc))/2.) *(diag.pr_tau2.t-diag.pr_tau1.t));
				  
  
  if (drnd() < weight) {
    updatestat(2,3) +=1;		//accepted
    diag.remove();
	global_weight *= weight_diff;
	//diag.printall();
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
#ifdef FP
  const double qcortemp = dqins;
  const double qctemp = dqins;
#elif BEC
  const double qcortemp = qcor;
  const double qctemp = qc;
#endif
  if(diag.propose_dq(qcortemp)!=0) {return -1;}
#ifdef BEC
  //Q Cut off
  if (vsq(diag.get_p(diag.pr_arc-1)-diag.get_p(diag.pr_arc)+diag.pr_q) > qctemp*qctemp){return -2;}
#endif
  updatestat(8,1) +=1;		//possible

  double weight = 1.;
  for (int i=0 ; i < (diag.get_link(diag.pr_arc) - diag.pr_arc); i++) {
    weight *= G0el(diag.get_p(diag.pr_arc + i) - diag.pr_q, diag.get_tfin(diag.pr_arc + i), diag.get_tinit(diag.pr_arc + i)); //G0el(neu)
    weight /= G0el(diag.get_p(diag.pr_arc + i), diag.get_tfin(diag.pr_arc + i), diag.get_tinit(diag.pr_arc + i));					//G0el(alt)
  }
  std:: array<double,3>  qold = diag.get_p(diag.pr_arc - 1)- diag.get_p(diag.pr_arc); 
  //Phonon term
  weight *= Dph(qold + diag.pr_q, diag.pr_taufin, diag.pr_tauin);
  weight /= Dph(qold, diag.pr_taufin, diag.pr_tauin);  
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

int DiagMC::insatend(){
  lu = "insatend()";
  updatestat(9,0) +=1;		//attempted
  if (diag.get_order() > maxord-1 && maxord != -1) {return -1;}
  if(diag.propose_insatend(taumax, dtins, dqins)!=0) {return -1;}
  updatestat(9,1) +=1;		//possible

  // weight part	
  double weight = G0el(diag.get_p(0), diag.pr_tau2.t, diag.pr_tau1.t);	//G0(new) second last and last vertex
  weight *= G0el(diag.pr_p, diag.pr_tau1.t, diag.pr_tauin);	//G0(new) 
  weight *= Dph(diag.pr_q, diag.pr_tau1.t, diag.pr_tauin);	//Phonon(new)
  weight *= Vq2(diag.pr_q);
  weight /= pow((2.*M_PI), 3);
  
  double weight_diff = weight;
    
  // a priori part
  if (diag.get_order() == 0) {weight /= fw;}   // for fake function
  
  weight *= Prae/Piae;	
  weight *= dtins; 					//sample second vertex to insert
  weight *= dtins; 					//select last tau
  weight *= pow(dqins,3);			//select phonon momentum
  
  if (drnd() < weight) {
  
    updatestat(9,3) +=1;		//accepted
    diag.insatend();
	global_weight *= weight_diff;
  }
  else {
    updatestat(9,2) +=1;		//rejected
  }  
  return 0;  
  
}

int DiagMC::rematend() {
  lu = "rematend()";
  updatestat(10,0) += 1;		//attempted
  if(diag.propose_rematend(dtins, dqins)!=0) {return -1;}
  updatestat(10,1) +=1;		//possible

  // weight part	
  std::array<double, 3> q = diag.get_p(diag.pr_arc-1)-diag.get_p(diag.pr_arc); 
  
  double weight = 1.;
  weight /= G0el(diag.get_p(diag.pr_arc), diag.pr_tau1.t, diag.pr_tauin);	//G0(old)
  weight /= G0el(diag.get_p(diag.pr_arc+1), diag.pr_tau2.t, diag.pr_tau1.t);	//G0(old)
  weight /= Dph(q, diag.pr_tau1.t, diag.pr_tauin);				  //Phonon(old)
  weight /= Vq2(q);
  weight *= pow((2.*M_PI), 3);
  double weight_diff = weight;

  // a priori part
  if (diag.get_order() == 1) {weight *= fw;}   // for fake function
  
  weight *= Piae/Prae;	
  weight /= dtins; 					//sample second vertex to insert
  weight /= dtins; 					//select last tau
  weight /= pow(dqins,3);			//select phonon momentum
	  
  
  if (drnd() < weight) {
    updatestat(10,3) +=1;		//accepted
    diag.rematend();
	global_weight *= weight_diff;
  }
  else {
    updatestat(10,2) +=1;		//rejected
  }
    
  return 0;
}