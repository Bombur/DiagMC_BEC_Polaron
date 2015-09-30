#include "Diagram.h"

//Definitions of Diagramm class

Diagram::Diagram() {
  order=0;
  
  times.reserve(1000);
  times.assign(2, std::vector<double>(2, 0.));
  
  elprop.reserve(1000);
  elprop.assign(1, std::vector<double>(3, 0.));
  
  phprop.reserve(1000);
  phprop.assign(1, std::vector<double>(3, 0.));
  
  pr_arc = 0;
  pr_tauin = 0.;
  pr_taufin = 0.;
  pr_tau1.assign(2, 0.);
  pr_tau2.assign(2, 0.);
  pr_q.assign(3, 0.);
  pr_p.assign(3, 0.);
}

Diagram::Diagram(const double & p, const double & tau, const std::function<double()> & rnd): order(0), drnd(rnd) {
  times.reserve(50);
  times.assign(2, std::vector<double>(2, 0.));
  times[0][1] = 1.;
  times[1][0] = tau;
    
  elprop.reserve(50);
  elprop.assign(1, std::vector<double>(3, 0.));
  elprop[0][0] = p;
  
  phprop.reserve(50);
  phprop.assign(1, std::vector<double>(3, 0.));
  
  pr_arc = 0;
  pr_tauin = 0.;
  pr_taufin = 0.;
  pr_tau1.assign(2, 0.);
  pr_tau2.assign(2, 0.);
  pr_q.assign(3, 0.);
  pr_p.assign(3, 0.);
}

void Diagram::set(const double & p, const double & tau, const std::function<double()> & rnd) {
  drnd = rnd;
  
  times[0][1] = 1.;
  times[1][0] = tau;
    
  elprop[0][0] = p;
}

void Diagram::random_arc() {
  pr_arc = (int)(drnd()*static_cast<double>((2*order)+1));
}


int Diagram::propose_insert(const double & dqins) {
  try{
	random_arc();										//arc to insert
	//if (order >1) {return -2;}						//step 2 just 1st order allowed
	
	pr_tauin = get_tinit(pr_arc);			
	pr_taufin = get_tfin(pr_arc);
	
	pr_tau1[0] = (drnd()*(pr_taufin - pr_tauin))+ pr_tauin;			//tau1
	pr_tau1[1] = static_cast<double>(pr_arc)+2.;
	pr_tau2[0] = (drnd()*(pr_taufin - pr_tau1[0]))+pr_tau1[0];							//tau2
	pr_tau2[1] = static_cast<double>(pr_arc)+1.;
	
	/*
	pr_q[0] = sqrt(-2.*log(drnd()))*cos(2.*M_PI*drnd())/sqrt(pr_tau2[0]-pr_tau1[0]);		//qx
	pr_q[1] = sqrt(-2.*log(drnd()))*sin(2.*M_PI*drnd())/sqrt(pr_tau2[0]-pr_tau1[0]);		//qy
	pr_q[2] = sqrt(-2.*log(drnd()))*cos(2.*M_PI*drnd())/sqrt(pr_tau2[0]-pr_tau1[0]);		//qz
	*/
	
	pr_q[0] = (drnd()-0.5)*dqins;		//qx
	pr_q[1] = (drnd()-0.5)*dqins;		//qy
	pr_q[2] = (drnd()-0.5)*dqins;		//qz
  
	if (vsq(pr_q) > pow(dqins*0.5,2)) {return -1;}
	
	pr_p = vsub(elprop[pr_arc], pr_q);
	 
	return 0;
  }	catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

int Diagram::propose_remove(const double & dqins) {
  try{
	random_arc();	//arc to remove
	if (pr_arc != (int)(times[pr_arc+1][1]+0.5)) {return -1;}
	if (order == 0) {return -1;}
	
	pr_tauin = get_tinit(pr_arc-1);
	pr_taufin = get_tfin(pr_arc+1);
	
	pr_tau1[0] = get_tinit(pr_arc);			//tau1
	pr_tau1[1] = -1.;
	pr_tau2[0] = get_tfin(pr_arc);							//tau2
	pr_tau2[1] = -1.;
	
	pr_q.assign(3,0.);		
	
	if (vsq(vsub(get_p(pr_arc-1), get_p(pr_arc))) > pow(dqins*0.5 , 2)){return -1;} 	//only remove in circle with diameter qins
	
	pr_p = get_p(pr_arc-1);
	return 0;
  }	catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

int Diagram::propose_swap() {
  random_arc();
  if (order == 0) {return -1;}
  if (pr_arc == 0 || (int)(times[pr_arc+1][1]+0.5) == 0) {return -1;}
  if ((int)(times[pr_arc][1]+0.5) == (pr_arc+1)) {return -1;}
  
  pr_tauin = get_tinit(pr_arc);
  pr_taufin = get_tfin(pr_arc);
  
  pr_tau1.assign(2,0.);
  pr_tau2.assign(2,0.);
  
  pr_q= vsub( vadd(get_q(pr_arc-1), get_q(pr_arc+1)), get_q(pr_arc));	
  
  pr_p =vsub(vadd(get_p(pr_arc -1), get_p(pr_arc+1)), get_p(pr_arc)); //pnew
  
  if ((((int)(times[pr_arc][1]+0.5)) > pr_arc) && (((int)(times[pr_arc+1][1]+0.5)) < pr_arc + 1) ) {
	return -2;
  }
  else if ((((int)(times[pr_arc][1]+0.5)) < pr_arc) && (((int)(times[pr_arc+1][1]+0.5)) > pr_arc + 1) ) {
	return 2;
  }
  else {
	return 0;		//wp tilde
  }
}


int Diagram::propose_ct(const double & taumax, const double & cor) {
  random_arc();
    
  pr_tau1[0] = get_tinit(pr_arc);
  pr_tau1[1] = 100.;
  
  pr_p = get_p(pr_arc);		//p in the beginning
  
  pr_arc += 1;
  
  pr_tauin = get_tinit(pr_arc); //tau old
  
  pr_tau2[1] = 100.;
  
  
  
  if ((int)(times[pr_arc][1]+0.5) == 0) {
	pr_tau2[0]=taumax;
	pr_taufin = ((drnd()-0.5)*cor)+ pr_tauin;		//tau new
	if (pr_taufin < pr_tau1[0] || pr_taufin > pr_tau2[0]) {return -1;}
	return 0; // last tau
  } 

  
  pr_tau2[0] = get_tfin(pr_arc);
  
  pr_taufin = ((drnd()-0.5)*cor) + pr_tauin;			//tau new
  
  if (pr_taufin < pr_tau1[0] || pr_taufin > pr_tau2[0]) {return -1;}
  
  pr_q = get_p(pr_arc);  // p in the middle
  
  if (get_link(pr_arc) > pr_arc) {return 1;} // opening arc
  if (get_link(pr_arc) < pr_arc) {return 2;} // closing arc
  
  return -1; 
}

int Diagram::propose_dq(const double & qcor) {
  try{
	if (order == 0) {return -1;}
	random_arc(); 
	if (pr_arc == 0) {return -1;}
	if (times[pr_arc][1] < static_cast<double>(pr_arc)) {
	  pr_arc = get_link(pr_arc);
	}
	
	pr_tauin = get_tinit(pr_arc);			
	pr_taufin = get_tinit(get_link(pr_arc));
	
	pr_tau1[0] = 9.;
	pr_tau1[1] = 9.;
	pr_tau2[0] = 9.;				
	pr_tau2[1] = 9.;
	

	pr_q[0] = (drnd()-0.5)* qcor;		//qx	
	pr_q[1] = (drnd()-0.5)* qcor;		//qy
	pr_q[2] = (drnd()-0.5)* qcor;		//qz
	
	pr_p.assign(3,0.);
	 
	return 0;
  }	catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

int Diagram::propose_insatend(const double & taumax, const double & dtins, const double & dqins) {
	pr_arc= get_link(0);				//choose last vertex
  
	pr_tauin = get_tinit(pr_arc);					//tau1
	pr_taufin = taumax;
	
	pr_tau1[0] = (drnd()*dtins)+ pr_tauin;			//tau2
	pr_tau1[1] = static_cast<double>(pr_arc);
	pr_tau2[0] = (drnd()*dtins)+pr_tau1[0];							//tau2
	pr_tau2[1] = 0;
	
	if (pr_tau2[0] > taumax) {return -1;}
	
	pr_q[0] = (drnd()-0.5)*dqins;		//qx
	pr_q[1] = (drnd()-0.5)*dqins;		//qy
	pr_q[2] = (drnd()-0.5)*dqins;		//qz
  
	if (vsq(pr_q) > pow(dqins*0.5,2)){return -1;}
	
	pr_p = vsub(elprop[0], pr_q);
	return 0;
}


int Diagram::propose_rematend(const double & dtins, const double & dqins) {
	if (order ==0) {return -1;}
	pr_arc= get_link(0)-2;				//choose second last vertex
	
	if (get_link(pr_arc) != pr_arc +1) {return -1;}
	  
	pr_tauin = get_tinit(pr_arc);					//tau1
	pr_taufin = -1.;
	
	pr_tau1[0] = get_tinit(pr_arc+1);			//tau2
	pr_tau1[1] = -1.;
	pr_tau2[0] = get_tfin(pr_arc+1);							//last tau
	pr_tau2[1] = -1.;
	
	if (((pr_tau1[0]-pr_tauin) > dtins) || ((pr_tau2[0]-pr_tau1[0]) > dtins)) {return -1;}
	if (vsq(vsub(get_p(pr_arc-1), get_p(pr_arc))) > pow(dqins*0.5 , 2)){return -1;}

	pr_q.assign(3, 0); 
	
	pr_p.assign(3, 0);
	 
	return 0;
}


void Diagram::insert() {
  try{
    if (pr_tau1.empty() || pr_tau2.empty() || pr_q.empty()) {throw no_propose();};
	std::vector< std::vector<double> >::iterator tit = times.begin() + pr_arc;		//iterators for different matrices
	std::vector< std::vector<double> >::iterator pit = phprop.begin() + pr_arc;
	std::vector< std::vector<double> >::iterator eit = elprop.begin() + pr_arc;
	
	for (int i = 0; i<((2*order)+1); i++) {
	  if ((int)(times[i][1]+0.5) > pr_arc) {times[i][1] +=2.;}
	}
	
	times.insert(tit+1, pr_tau1);										//insert tau1 and tau2
	times.insert(tit+2, pr_tau2);
	
	
	phprop.insert(pit+1, std::move(vadd(pr_q, get_q(pr_arc))));							//insert q+oldq in arc+1
	phprop.insert(pit+2, std::move(get_q(pr_arc)));								//insert oldq in arc+2
	
	elprop.insert(eit+1, std::move(vsub(get_p(pr_arc), pr_q)));					//insert oldp-q in arc+1
	elprop.insert(eit+2, std::move(get_p(pr_arc)));								//insert oldp in arc+2	  
 
	
	order+=1;

  }	catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }

}

void Diagram::remove() {
   try{
    if (pr_tau1.empty() || pr_tau2.empty() || pr_q.empty()) {throw no_propose();}
	std::vector< std::vector<double> >::iterator tit = times.begin() + pr_arc;		//iterators for different matrices
    std::vector< std::vector<double> >::iterator pit = phprop.begin() + pr_arc;
	std::vector< std::vector<double> >::iterator eit = elprop.begin() + pr_arc;
	  
	times.erase(tit, tit+2);	  
	
	phprop.erase(pit, pit+2);
      
	elprop.erase(eit, eit+2);  
  
	order-=1;
	  
	for (int i = 0; i<((2*order)+1); i++) {
	  if ((int)(times[i][1]+0.5) > pr_arc) {times[i][1] -=2.;}
	}

  } catch (std::exception& e) {
      std::cerr << e.what() << std::endl;
      exit(EXIT_FAILURE);
  }

}

int Diagram::set_tau(double tau) {
  if (tau < times[2*order][0]) {return -1;}
  times[(2*order)+1][0] = tau;
  
  return 0;
}


void Diagram::swap() {
  double tmp = times[(int)(times[pr_arc][1]+0.5)][1];
  times[(int)(times[pr_arc][1]+0.5)][1] = times[(int)(times[pr_arc+1][1]+0.5)][1];
  times[(int)(times[pr_arc+1][1]+0.5)][1] = tmp;
  
  tmp = times[pr_arc][1];			//swap links
  times[pr_arc][1] = times[pr_arc+1][1];
  times[pr_arc+1][1] = tmp;
    
  phprop[pr_arc] = pr_q;
  
  elprop[pr_arc] = pr_p;
}


void Diagram::ct() {
  times[pr_arc][0]= pr_taufin;
}

void Diagram::dq() {
  for (int i=0 ; i < (get_link(pr_arc) - pr_arc); i++) {
    phprop[pr_arc + i][0] += pr_q[0];
    phprop[pr_arc + i][1] += pr_q[1];
    phprop[pr_arc + i][2] += pr_q[2];
    
    elprop[pr_arc + i][0] -= pr_q[0];
    elprop[pr_arc + i][1] -= pr_q[1];
    elprop[pr_arc + i][2] -= pr_q[2];
  }

}

void Diagram::insatend() {	
	times.push_back(pr_tau1);						//insert tau2 and last tau
	times.push_back(pr_tau2);
	
	times[0][1] += 2.;				//link to last tau
	times[pr_arc][1] += pr_arc + 1.;			// opening arc
	
	
	phprop.push_back(pr_q);							
	phprop.push_back(get_q(0));
	
	elprop.push_back(pr_p);
	elprop.push_back(get_p(0));  
 
	
	order+=1;
}


void Diagram::rematend() {	
	times.pop_back();						
	times.pop_back();
	
	times[0][1] -= 2.;				//link to last tau
	times[pr_arc][1] = 0.;			// last tau
	
	phprop.pop_back();							
	phprop.pop_back();
	
	elprop.pop_back();
	elprop.pop_back();  
 
	
	order-=1;
}












