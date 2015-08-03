#include "Diagram.h"

class dnf_vecsize: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Vector sizes do not fit!";
  }
};

class no_propose: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Nothing proposed!";
  }
};

std::vector<double> operator+(const std::vector<double> & vec1, const std::vector<double> & vec2){
  try {
	if (vec1.size() != vec2.size()) {throw dnf_vecsize();}
	std::vector<double> tmp(vec1.size());
	for (int i =0; i<tmp.size; i++){
	  tmp[i]= vec1[i] + vec2[i];
	}
	return tmp;
  }
  catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

std::vector<double> operator-(const std::vector<double> & vec1, const std::vector<double> & vec2){
  try {
	if (vec1.size() != vec2.size()) {throw dnf_vecsize();}
	std::vector<double> tmp(vec1.size());
	for (int i =0; i<tmp.size; i++){
	  tmp[i]= vec1[i] - vec2[i];
	}
	return tmp;
  }
  catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

double square(const std::vector<double> & vec1){
  try {
	double tmp=0;
	for (int i =0; i<vec1.size; i++){
	  tmp += pow(vec1[i],2);
	}
	return tmp;
  }
  catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

Diagram::Diagram(double p, double tau, double my, double omega, double, alph, int ord, std::function<double()> rnd): mu(my), wp(omega), alpha(alph), order(ord), drnd(rnd) {
  times.reserve(50);
  times.assign(2, std::vector<double>(2, 0));
  times[0][1] = 1;
  times[1][0] = tau;
    
  elprop.reserve(50);
  elprop.assign(1, std::vector<double>(3, 0));
  elprop[0][0] = p;
  
  phprop.reserve(50);
  phprop.assign(1, std::vector<double>(3, 0));
  
}

void Diagram::random_arc() {
  pr_arc = (int)(drnd()*(2*old.get_order()+1));
}


void Diagram::propose_insert() {
  try{
	random_arc();										//arc to insert
	
	pr_tau1[0] = (drnd()*(get_tfin(pr_arc)- get_tinit(pr_arc)))+ get_tinit(pr_arc);			//tau1
	pr_tau1[1] = pr_arc+2;
	pr_tau2[0] = (drnd()*(get_tfin(pr_arc)-pr_tau1[0]))+pr_tau1[0];							//tau2
	pr_tau2[1] = pr_arc+1;
	
	pr_q[0] = sqrt(-2*log(drnd()))*cos(2*M_PI*drnd())/sqrt(pr_tau2[0]-pr_tau1[0]);		//qx
	pr_q[1] = sqrt(-2*log(drnd()))*sin(2*M_PI*drnd())/sqrt(pr_tau2[0]-pr_tau1[0]);		//qy
	pr_q[2] = sqrt(-2*log(drnd()))*cos(2*M_PI*drnd())/sqrt(pr_tau2[0]-pr_tau1[0]);		//qz
	
  }	catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

int Diagram::propose_remove() {
  try{
	random_arc();	//arc to insert
	
	if (pr_arc != times[pr_arc+1][1]) {return 0;}
	
	pr_tau1[0] = get_tfin(pr_arc);			//tau1
	pr_tau1[1] = -1;
	pr_tau2[0] = get_tfin(pr_arc);							//tau2
	pr_tau2[1] = -1;
	
	pr_q = get_q(pr_arc);		//qx
		
	return 1;
  }	catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

}


double Diagram::high_weigth() {
  return alpha*exp(((square(get_p(pr_arc)-pr_q)/2) - mu)(pr_tau1[0]-pr_tau2[0]) - wp*(pr_tau2[0]-pr_tau1[0])) / square(pr_q);
}

double Diagram::low_weigth() {
  return exp(((square(get_p(pr_arc))/2) - mu)(pr_tau1[0]-pr_tau2[0]));
}

double Diagram::P_hilo(double Prem){
  return Prem/(get_order()+1);
}

double Diagram::P_lohi(double Pins){
  return Pins/(2*get_order()+1) /(get_tfin(pr_arc)-get_tinit(pr_arc)) /(get_tfin(pr_arc)-pr_tau1[0]) *pow((pr_tau2[0]-pr_tau1[0])/2/M_PI, 3/2) * exp((square(get_p(pr_arc)-pr_q)/2) (pr_tau1[0]-pr_tau2[0]));
}

void Diagram::insert() {
  try{
    if (pr_tau1.empty() || pr_tau2.empty() || pr_q.empty()) {throw no_propose();};
	std::vector< std::vector<double> >::iterator tit = times.begin() + pr_arc;		//iterators for different matrices
	std::vector< std::vector<double> >::iterator pit = phprop.begin() + pr_arc;
	std::vector< std::vector<double> >::iterator eit = elprop.begin() + pr_arc;
	
	times.insert(tit+1, std::move(pr_tau1));										//insert tau1 and tau2
	times.insert(tit+2, std::move(pr_tau2));
	
	phprop.insert(pit+1, std::move(pr_q + get_q(pr_arc)));							//insert q+oldq in arc+1
	phprop.insert(pit+2, std::move(get_q(pr_arc)));								//insert oldq in arc+2
	
	elprop.insert(eit+1, std::move(get_p(pr_arc) - get_q(pr_arc+1)));					//insert oldp-q in arc+1
	elprop.insert(eit+2, std::move(get_p(pr_arc)));								//insert oldp in arc+2	  
  
	order+=1;
	
  }	catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }

}

void Diagram::remove() {
   try{
    if (pr_tau1.empty() || pr_tau2.empty() || pr_q.empty()) {throw no_propose();};
      std::vector< std::vector<double> >::iterator tit = times.begin() + pr_arc;		//iterators for different matrices
      std::vector< std::vector<double> >::iterator pit = phprop.begin() + pr_arc;
      std::vector< std::vector<double> >::iterator eit = elprop.begin() + pr_arc;
	
      times.erase(tit, tit+1);
	
      phprop.erase(pit, pit+1);
      
      elprop.erase(eit, eit+1);  
  
      order-=1;
	
  } catch (std::exception& e) {
      std::cerr << e.what() << std::endl;
      exit(EXIT_FAILURE);
  }

}

