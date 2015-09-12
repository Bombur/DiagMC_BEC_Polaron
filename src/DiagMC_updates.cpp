#include "DiagMC.h"

 
int DiagMC::change_tau() {
  double ntau = - log(drnd())/E;		//new tau
  updatestat(0,0) +=1;						//attempted
  
  if (diag.get_order() != 0) {return -1;}
    
  if (ntau < taumax) {
	if (diag.set_tau(ntau) == 0) {
	  updatestat(0,1) +=1;						//possible
	  updatestat(0,3) +=1;						//accepted
	  //tau = ntau;
	  return 0;
	}
  }
  return -1;
}



int DiagMC::ct(){
  updatestat(7,0) +=1;		//attempted
  int which = diag.propose_ct(taumax, ctcor);
  double weight;
  if(which == 0) {
	weight = G0el(diag.pr_p, diag.pr_taufin, diag.pr_tauin);
  } else if (which==1){
	weight = G0el(diag.pr_p, diag.pr_taufin, diag.pr_tauin)/G0el(diag.pr_q, diag.pr_taufin, diag.pr_tauin) * Dph(diag.pr_taufin, diag.pr_tauin);
  }
  else {return -1;}
  
  updatestat(7,1) +=1;		//possible
  if (drnd() < weight) {
    updatestat(7,3) +=1;		//accepted
    diag.ct();
	//std::cout << diag.get_tau() << std:: endl;
  }
  else {
    updatestat(7,2) +=1;		//rejected
  }  
  
  return 0; 
}

int DiagMC::insert() {
  updatestat(1,0) +=1;		//attempted
  if(diag.propose_insert()!=0) {return -1;}
  updatestat(1,1) +=1;		//possible
  
  double weight = (alpha/vsq(diag.pr_q)*G0el(vsub(diag.pr_p,diag.pr_q), diag.pr_tau2[0], diag.pr_tau1[0])*Dph(diag.pr_tau2[0], diag.pr_tau1[0]))
				  *(Prem /double(diag.get_order()+1))
				  /(G0el(diag.pr_p, diag.pr_tau2[0], diag.pr_tau1[0]))
				  /(Pins/(2*double(diag.get_order())+1) /(diag.pr_taufin-diag.pr_tauin) /(diag.pr_taufin-diag.pr_tau1[0]) *pow((diag.pr_tau2[0]-diag.pr_tau1[0])/2/M_PI, 3/2) * exp(-vsq(diag.pr_q)/2 *(diag.pr_tau2[0]-diag.pr_tau1[0])));

  if (drnd() < weight) {
    updatestat(1,3) +=1;		//accepted
    diag.insert();
  }
  else {
    updatestat(1,2) +=1;		//rejected
  }  
  return 0;    
}

int DiagMC::remove() {
  updatestat(2,0) += 1;		//attempted
  if(diag.propose_remove()!=0) {return -1;}
  updatestat(2,1) +=1;		//possible
  
  double weight = 1/(alpha/vsq(diag.pr_q)*G0el(vsub(diag.pr_p,diag.pr_q), diag.pr_tau2[0], diag.pr_tau1[0])*Dph(diag.pr_tau2[0], diag.pr_tau1[0]))
				  /(Prem /double(diag.get_order()+1))
				  *(G0el(diag.pr_p, diag.pr_tau2[0], diag.pr_tau1[0]))
				  *(Pins/(2*double(diag.get_order())+1) /(diag.pr_taufin-diag.pr_tauin) /(diag.pr_taufin-diag.pr_tau1[0]) *pow((diag.pr_tau2[0]-diag.pr_tau1[0])/2/M_PI, 3/2) * exp(-vsq(diag.pr_q)/2 *(diag.pr_tau2[0]-diag.pr_tau1[0])));

  if (drnd() < weight) {
    updatestat(2,3) +=1;		//accepted
    diag.remove();
  }
  else {
    updatestat(2,2) +=1;		//rejected
  }
    
  return 0;
}

int DiagMC::swap() {
  updatestat(3,0) += 1; //attempted
  updatestat(4,0) += 1;
  updatestat(5,0) += 1;
  updatestat(6,0) += 1;
  int which = diag.propose_swap();
  if (which == -1) {return -1;} //not possible
  
  double weight = G0el(diag.pr_p, diag.pr_taufin, diag.pr_tauin)/G0el(diag.get_p(diag.pr_arc), diag.pr_taufin, diag.pr_tauin);
  
  // both vertices are open or closed
  if (which == 0) {
	updatestat(4,1) +=1; //possible
	if (drnd() < weight) {
	  updatestat(4,3) +=1; //accepted
	  diag.swap();
	}
	else {
	  updatestat(4,2)+= 1; 	// rejected
	}
  }
  
  // first vertex is start point, second one is end point 
  if (which == -2) {
	updatestat(5,1) +=1; //possible
	weight /= pow(Dph(diag.pr_taufin, diag.pr_tauin),2);
	if (drnd() < weight) {
	  updatestat(5,3) +=1; //accepted
	  diag.swap();
	}
	else {
	  updatestat(5,2)+= 1; 	// rejected
	}
  }
  
   // first vertex is end point, second one is start point 
  if (which == 2) {
	updatestat(6,1) +=1; //possible
	weight *= pow(Dph(diag.pr_taufin, diag.pr_tauin),2);
	if (drnd() < weight) {
	  updatestat(6,3) +=1; //accepted
	  diag.swap();
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
  updatestat(8,0) +=1;		//attempted
  if(diag.propose_dq(qcor)!=0) {return -1;}
  updatestat(8,1) +=1;		//possible
  
  double weight = 1;
  for (int i=0 ; i < (diag.get_link(diag.pr_arc) - diag.pr_arc); i++) {
    weight *= G0el(vsub(diag.get_p(diag.pr_arc + i), diag.pr_q), diag.get_tfin(diag.pr_arc + i), diag.get_tinit(diag.pr_arc + i));
    weight /= G0el(diag.get_p(diag.pr_arc + i), diag.get_tfin(diag.pr_arc + i), diag.get_tinit(diag.pr_arc + i));
  }
  

  if (drnd() < weight) {
    updatestat(8,3) +=1;		//accepted
    diag.dq();
  }
  else {
    updatestat(8,2) +=1;		//rejected
  }  
  return 0;    
}
	


void DiagMC::measure(const int & whichstep) {
  Data((int)(diag.get_tau()/taumax*taubin), 0) +=1;
  if (get_order() < 2) {Data((int)(diag.get_tau()/taumax*taubin), get_order()+1) += 1;}
  if (get_order() ==2 && vsq(diag.get_q(2)) < 0.00000000000001) {	
	Data((int)(diag.get_tau()/taumax*taubin), get_order()+1) += 1;
  }
  
  if (get_order() < orderstat.size()){
	orderstat(get_order()) +=1;
  }
   
} 

