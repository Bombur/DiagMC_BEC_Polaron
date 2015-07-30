#include "Diagram.h"

class dnf_vecsize: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Vector sizes do not fit!";
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

void Diagram::propose() {
  try{
	arc = (int)(drnd()*(2*old.get_order()+1));								//arc to insert
	
	tau1[0] = (drnd()*(get_tfin(arc)- get_tinit(arc)))+ get_tinit(arc);			//tau1
	tau1[1] = arc+2;
	tau2[0] = (drnd()*(get_tfin(arc)-tau1[0]))+tau1[0];							//tau2
	tau2[1] = arc+1;
	
	q[0] = sqrt(-2*log(drnd()))*cos(2*M_PI*drnd())/sqrt(tau2[0]-tau1[0]);		//qx
	q[1] = sqrt(-2*log(drnd()))*sin(2*M_PI*drnd())/sqrt(tau2[0]-tau1[0]);		//qy
	q[2] = sqrt(-2*log(drnd()))*cos(2*M_PI*drnd())/sqrt(tau2[0]-tau1[0]);		//qz
	
  }	catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

double Diagram::high_weigth() {
  return alpha*exp(((square(get_p(arc)-q)/2) - mu)(tau1[0]-tau2[0]) - wp*(tau2[0]-tau1[0])) / square(q);
}

double Diagram::low_weigth() {
  return exp(((square(get_p())/2) - mu)(tau1[0]-tau2[0]));
}

double Diagram::P_hilo(double Prem){
  return Prem/(get_order()+1);
}

double Diagram::P_lohi(double Pins){
  return Pins/(2*get_order()+1) /(get_tfin(arc)-get_tinit(arc)) /(get_tfin(arc)-tau1[0]) *pow((tau2[0]-tau1[0])/2/M_PI, 3/2) * exp((square(get_p(arc)-q)/2) (tau1[0]-tau2[0]));
}

void Diagram::insert() {
  try{
	std::vector< std::vector<double> >::iterator tit = times.begin() + arc;		//iterators for different matrices
	std::vector< std::vector<double> >::iterator eit = elprop.begin() + arc;
	std::vector< std::vector<double> >::iterator pit = phprop.begin() + arc;
	
	times.insert(tit+1, std::move(tau1));										//insert tau1 and tau2
	times.insert(tit+2, std::move(tau2));
	
	phprop.insert(pit+1, std::move(q + get_q(arc)));							//insert q+oldq in arc+1
	phprop.insert(pit+2, std::move(get_q(arc)));								//insert oldq in arc+2
	
	elprop.insert(pit+1, std::move(get_p(arc) - get_q(arc+1)));					//insert oldp-q in arc+1
	elprop.insert(pit+2, std::move(get_p(arc)));								//insert oldp in arc+2	  
  
	order+=1;
	
  }	catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }

}
