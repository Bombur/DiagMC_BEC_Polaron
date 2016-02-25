#include "Diagram.h"

//Definitions of Diagramm class

Diagram::Diagram() {
  order=0;
  
  vertex fv;
  fv.t = 0.;
  fv.link = 0;
  
  times.reserve(1000);
  times.assign(2, fv);
  
  std::array<double, 3> zero;
  zero.fill(0.);
#ifndef NCHECK
  phprop.reserve(1000);
  phprop.assign(1, zero);
 // qs.reserve(500);
#endif 
  
  elprop.reserve(1000);
  elprop.assign(1, zero);
  
  befcap = times.capacity();
  
  pr_arc = 0;
  pr_tauin = 0.;
  pr_taufin = 0.;
  
  pr_tau1.t = 0.;
  pr_tau1.link = 0;
  pr_tau2 = pr_tau1;
  
  pr_q.fill(0.);
  pr_p.fill(0.);
}

void Diagram::set(const double & p, const double & tau, const std::function<double()> & rnd) {
  drnd = rnd;
  
  times[0].link = 1;
  times[1].t = tau;

  elprop[0][0] = p;
}

void Diagram::random_arc() {
  pr_arc = (int)(drnd()*static_cast<double>((2*order)+1));
}


int Diagram::propose_insert(const double & dqins, const double & sigfac) {
	random_arc();										//arc to insert
	//if (order >1) {return -2;}						//step 2 just 1st order allowed
	
#ifdef SELFENERGY
	if (order > 0 && (pr_arc == 0 || pr_arc == 2*order)) {return -1;} 
#endif	
	pr_tauin = times[pr_arc].t;	  // get_tauinit
	pr_taufin = times[pr_arc+1].t; // get_taufin
	
	pr_tau1.t = (drnd()*(pr_taufin - pr_tauin))+ pr_tauin;			//tau1
	pr_tau1.link = pr_arc+2;
	pr_tau2.t = (drnd()*(pr_taufin - pr_tau1.t))+pr_tau1.t;							//tau2
	pr_tau2.link = pr_arc+1;
	

#ifdef QGAUSS
	//Gauss q	
	const double qr= (dqins/sigfac)*sqrt(-2.* log(drnd()));
	const double qtheta = drnd() *M_PI;
	const double qphi = drnd() *2. *M_PI;
	pr_q[0] = qr * cos(qtheta)*cos(qphi);	//qx
	pr_q[1] = qr * cos(qtheta)*sin(qphi);		//qy
	pr_q[2] = qr * sin(qtheta);		//qz
#else
	// linear chosen q
	pr_q[0] = (drnd()-0.5)*2*dqins;		//qx
	pr_q[1] = (drnd()-0.5)*2*dqins;		//qy
	pr_q[2] = (drnd()-0.5)*2*dqins;		//qz
#endif  
	if (vsq(pr_q) > dqins*dqins) {return -1;}
	
	pr_p = elprop[pr_arc] - pr_q;
	 
	return 0;
}

int Diagram::propose_remove(const double & dqins) {
	random_arc();	//arc to remove
	if (pr_arc != times[pr_arc+1].link) {return -1;}
	if (order == 0) {return -1;}
	
	pr_tauin = times[pr_arc-1].t;   // get_tinit(pr_arc-1);
	pr_taufin = times[pr_arc+2].t ;  // get_tfin(pr_arc+1);
	
	pr_tau1.t = times[pr_arc].t ;   //get_tinit(pr_arc);			//tau1
	pr_tau1.link  = -1;
	pr_tau2.t  = times[pr_arc+1].t;   // get_tfin(pr_arc);							//tau2
	pr_tau2.link = -1;
	
	pr_q.fill(0);		
	
	if (vsq(elprop[pr_arc-1] - elprop[pr_arc]) > dqins*dqins){return -1;} 	//only remove in circle with diameter qins
	
	pr_p = elprop[pr_arc-1];
	return 0;
}

int Diagram::propose_swap() {
  random_arc();
  if (order == 0) {return -1;}
  if (pr_arc == 0 || times[pr_arc+1].link == 0) {return -1;} // first and second last
  
  if (times[pr_arc].link == (pr_arc+1)) {return -1;} //small loop
  
  pr_tauin = times[pr_arc].t;   // get_tinit(pr_arc);
  pr_taufin = times[pr_arc+1].t ;  // get_tfin(pr_arc);
  
  
  pr_p = elprop[pr_arc -1] + elprop[pr_arc+1] - elprop[pr_arc];  
  
  pr_q= elprop[0] - pr_p ; 
  
  
  //first opening second closing 
  if (((times[pr_arc].link) > pr_arc) && ((times[pr_arc+1].link) < (pr_arc + 1)) ) {
#ifdef SELFENERGY
	if (vsq(pr_p - elprop[0]) < 0.00000000001) {return -1;} // check to not open SE Diagram
#endif
	return -2;
  }
  //first closing second opening
  else if (((times[pr_arc].link) < pr_arc) && ((times[pr_arc+1].link) > (pr_arc + 1)) ) {
	return 2;
  }
  else {
	return 0;	
  }
}


int Diagram::propose_ct(const double & taumax, const double & cor) {
  random_arc();
    
  pr_tau1.t = times[pr_arc].t;  //get_tinit(pr_arc);
  pr_tau1.link = 100;
  
  pr_p = elprop[pr_arc];		//p in the beginning
  
  pr_arc += 1;
  
  pr_tauin = times[pr_arc].t; //tau old
  
  pr_tau2.link = 100;
  
  
  
  if (times[pr_arc].link == 0) {
	pr_tau2.t=taumax;
	pr_taufin = ((drnd()-0.5)*cor)+ pr_tauin;		//tau new
	if (pr_taufin < pr_tau1.t || pr_taufin > pr_tau2.t) {return -1;}
	return 0; // last tau
  } 

  
  pr_tau2.t = times[pr_arc+1].t;   //get_tfin(pr_arc);
  
  pr_taufin = ((drnd()-0.5)*cor) + pr_tauin;			//tau new
  
  if (pr_taufin < pr_tau1.t || pr_taufin > pr_tau2.t) {return -1;}
  
  pr_q = elprop[pr_arc];   // get_p(pr_arc);  // p in the middle
  
  if (times[pr_arc].link > pr_arc) {return 1;} // opening arc
  if (times[pr_arc].link < pr_arc) {return 2;} // closing arc
  
  return -1; 
}

int Diagram::propose_dq(const double & qcor) {
	if (order == 0) {return -1;}
	random_arc(); 
	if (pr_arc == 0) {return -1;}
	if (times[pr_arc].link < pr_arc) {
	  pr_arc = times[pr_arc].link;
	}
	
	pr_tauin = times[pr_arc].t;   //get_tinit(pr_arc);			
	pr_taufin = times[times[pr_arc].link].t;   //get_tinit(get_link(pr_arc));
	
	pr_q[0] = (drnd()-0.5)* qcor;		//qx	
	pr_q[1] = (drnd()-0.5)* qcor;		//qy
	pr_q[2] = (drnd()-0.5)* qcor;		//qz
	
	pr_p.fill(0.);
	
	return 0;
}

int Diagram::propose_insatend(const double & taumax, const double & dtins, const double & dqins) {
	pr_arc= times[0].link;				//choose last vertex
  
	pr_tauin = times[pr_arc].t;   //get_tinit(pr_arc);     tau1					
	pr_taufin = taumax;  
	
	pr_tau1.t = (drnd()*dtins)+ pr_tauin;			// tau2
	pr_tau1.link = pr_arc;
	pr_tau2.t = (drnd()*dtins)+pr_tau1.t;			// tau
	pr_tau2.link = 0;
	
	if (pr_tau2.t > taumax) {return -1;}
	
	pr_q[0] = (drnd()-0.5)*dqins;		//qx
	pr_q[1] = (drnd()-0.5)*dqins;		//qy
	pr_q[2] = (drnd()-0.5)*dqins;		//qz
  
	if (vsq(pr_q) > dqins*dqins*0.25){return -1;}
	
	pr_p = elprop[0] - pr_q;
	return 0;
}


int Diagram::propose_rematend(const double & dtins, const double & dqins) {
	if (order ==0) {return -1;}
	pr_arc= times[0].link-2;				//choose third last vertex
	
	if (times[pr_arc].link != (pr_arc +1)) {return -1;}
	  
	pr_tauin = times[pr_arc].t;   //get_tinit(pr_arc);					//tau1
	pr_taufin = -1.;
	
	pr_tau1.t = times[pr_arc+1].t;   //get_tinit(pr_arc+1);			//tau2
	pr_tau1.link = -1;
	pr_tau2.t = times[pr_arc+2].t;  //get_tfin(pr_arc+1);							//last tau
	pr_tau2.link = -1;
	
	if (((pr_tau1.t-pr_tauin) > dtins) || ((pr_tau2.t-pr_tau1.t) > dtins)) {return -1;}
	if (vsq(elprop[pr_arc-1] - elprop[pr_arc]) > dqins*dqins*0.25 ){return -1;}

	pr_q.fill(0.); 
	
	pr_p.fill(0.);
	 
	return 0;
}


void Diagram::insert() {
	for (int i = 0; i<((2*order)+1); i++) {
	  if (times[i].link > pr_arc) {times[i].link +=2;}
	}
	
	times.insert(times.begin() + pr_arc + 1, pr_tau1);										//insert tau1 and tau2
	times.insert(times.begin() + pr_arc + 2, pr_tau2);
		
	elprop.insert(elprop.begin() + pr_arc + 1, std::move(elprop[pr_arc]- pr_q));					//insert oldp-q in arc+1
	elprop.insert(elprop.begin() + pr_arc + 2, std::move(elprop[pr_arc]));								//insert oldp in arc+2	  
	
#ifndef NCHECK
	phprop.insert(phprop.begin() + pr_arc + 1, std::move(pr_q + phprop[pr_arc]));							//insert q+oldq in arc+1
	phprop.insert(phprop.begin() + pr_arc + 2, std::move(phprop[pr_arc]));								//insert oldq in arc+2
#endif
	
	order+=1;

}

void Diagram::remove() {
	std::vector< vertex >::iterator tit = times.begin() + pr_arc;		//iterators for different matrices
	std::vector< std::array<double,3> >::iterator eit = elprop.begin() + pr_arc;
	  
	times.erase(tit, tit+2);	  
	elprop.erase(eit, eit+2); 
	
#ifndef NCHECK
	std::vector< std::array<double,3> >::iterator pit = phprop.begin() + pr_arc;
	phprop.erase(pit, pit+2);

#endif
  
	order-=1;
	  
	for (int i = 0; i<((2*order)+1); i++) {
	  if (times[i].link > pr_arc) {times[i].link -=2;}
	}

}

int Diagram::set_tau(double tau) {
  if (tau < times[2*order].t) {return -1;}
  times[(2*order)+1].t = tau;
  
  return 0;
}


void Diagram::swap() {
  int tmp = times[times[pr_arc].link].link;
  times[times[pr_arc].link].link = times[times[pr_arc+1].link].link;
  times[times[pr_arc+1].link].link = tmp;
  
  tmp = times[pr_arc].link;			//swap links
  times[pr_arc].link = times[pr_arc+1].link;
  times[pr_arc+1].link = tmp;
    
  elprop[pr_arc] = pr_p;
#ifndef NCHECK
  phprop[pr_arc] = pr_q;
#endif 
}


void Diagram::ct() {
  times[pr_arc].t= pr_taufin;
}

void Diagram::dq() {
  for (int i=0 ; i < (times[pr_arc].link - pr_arc); i++) {
#ifndef NCHECK
    phprop[pr_arc + i][0] += pr_q[0];
    phprop[pr_arc + i][1] += pr_q[1];
    phprop[pr_arc + i][2] += pr_q[2];
#endif
    
    elprop[pr_arc + i][0] -= pr_q[0];
    elprop[pr_arc + i][1] -= pr_q[1];
    elprop[pr_arc + i][2] -= pr_q[2];
  }

}

void Diagram::insatend() {	
	times.push_back(pr_tau1);						//insert tau2 and last tau
	times.push_back(pr_tau2);
	
	times[0].link += 2;				//link to last tau
	times[pr_arc].link += pr_arc + 1;			// opening arc
	
#ifndef NCHECK
	phprop.push_back(pr_q);							
	phprop.push_back(phprop[0]);
#endif
	
	elprop.push_back(pr_p);
	elprop.push_back(elprop[0]);  
	
	order+=1;
}


void Diagram::rematend() {	
	times.pop_back();						
	times.pop_back();
	
	times[0].link -= 2;				//link to last tau
	times[pr_arc].link = 0;			// last tau
	
	elprop.pop_back();
	elprop.pop_back();  
	
#ifndef NCHECK
	phprop.pop_back();							
	phprop.pop_back();
#endif
	
	order-=1;
}











