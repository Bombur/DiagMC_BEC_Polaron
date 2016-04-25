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

void Diagram::linear_q(const double & dqins){
  // linear chosen q
  for (int i=0; i<3; i++){
	pr_q[i] = (drnd()-0.5)*2*dqins;		
  }
}


void Diagram::gauss_q(const double& sigma, const double & dqins) {
  // Box Mueller from Krauth
  double Gamma = 0.;
  double x = 0.;
  double y = 0.;
  do{
	x = drnd()*2. - 1.;
	y = drnd()*2. - 1.;
	Gamma = x*x + y*y;
  } while (Gamma > 1);
  double Gamma2 = sigma*sqrt(2*(-log(Gamma))/Gamma);
  pr_q[0] = Gamma2*x;
  pr_q[1] = Gamma2*y;
  
  do{
	x = drnd()*2. - 1.;
	y = drnd()*2. - 1.;
	Gamma = x*x + y*y;
  } while (Gamma > 1.);
  Gamma2 = sigma*sqrt(2.*(-log(Gamma))/Gamma);
  pr_q[2] = Gamma2*x;
  
  //Now we revert it to get linear distributed qs also Krauth
  double sum = vsq(pr_q);
  Gamma = pow(drnd(),1./3.);
  for (int i = 0; i<3; i++) {
	pr_q[i] *= dqins*Gamma/sqrt(sum);
  }	
}



int Diagram::propose_swap() {
  //Fast Rejections
  if (order == 0) {return -1;}
  random_arc();
  if (pr_arc == 0 || times[pr_arc+1].link == 0) {return -1;} // first and second last
  if (times[pr_arc].link == (pr_arc+1)) {return -1;} //small loop
  
  //Taus
  pr_tauin = times[pr_arc].t;   // get_tinit(pr_arc);
  pr_taufin = times[pr_arc+1].t ;  // get_tfin(pr_arc);
  
  
  pr_p = elprop[pr_arc -1] + elprop[pr_arc+1] - elprop[pr_arc];  
  pr_q= elprop[0] - pr_p ; 
  
  
  //first opening second closing 
  if (((times[pr_arc].link) > pr_arc) && ((times[pr_arc+1].link) < (pr_arc + 1)) ) {
#ifdef SELFENERGY
	if (vsq(pr_p - elprop[0]) < 0.00000000001) {return -1;} // check to not open SE Diagram
#endif
	return 0;
  }
  //first closing second opening
  else if (((times[pr_arc].link) < pr_arc) && ((times[pr_arc+1].link) > (pr_arc + 1)) ) {
	return 1;
  }
  //opening opening
  else if (((times[pr_arc].link) > pr_arc) && ((times[pr_arc+1].link) > (pr_arc + 1)) ) {
	return 2;	
  }
  //closing closing
  else if (((times[pr_arc].link) < pr_arc) && ((times[pr_arc+1].link) < (pr_arc + 1)) ) {
	return 3;	
  }
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
  times[pr_arc].t= pr_tau2.t;
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

void Diagram::inscrossed() {
	for (int i = 0; i<((2*order)+1); i++) {
	  if (times[i].link > pr_arc) {times[i].link +=2;} 
	  else if (times[i].link == pr_arc) {times[i].link +=1;}
	}
	//insert tau1 and tau2 around pr_arc
	times.insert(times.begin() + pr_arc, pr_tau1);
	times.insert(times.begin() + pr_arc + 2, pr_tau2);
		
	//insert first new pright than pleft
	elprop.insert(elprop.begin() + pr_arc, std::move(elprop[pr_arc]-pr_q));	
	elprop.insert(elprop.begin() + pr_arc, std::move(elprop[pr_arc-1]- pr_q));
	
#ifndef NCHECK
	phprop.insert(phprop.begin() + pr_arc, std::move(pr_q + phprop[pr_arc]));	
	phprop.insert(phprop.begin() + pr_arc, std::move(pr_q + phprop[pr_arc-1]));
#endif
	order+=1;
}

void Diagram::remcrossed() {
	times.erase(times.begin() + pr_arc + 1);
	times.erase(times.begin() + pr_arc - 1);	
	elprop.erase(elprop.begin() + pr_arc - 1, elprop.begin() + pr_arc + 1); 
	
#ifndef NCHECK
	phprop.erase(phprop.begin() + pr_arc -1, phprop.begin() + pr_arc + 1);
#endif
  
	order-=1;
	  
	for (int i = 0; i<((2*order)+1); i++) {
	  if (times[i].link > pr_arc) {times[i].link -=2;}
	  else if (times[i].link == pr_arc) {times[i].link -= 1;}
	}

}











