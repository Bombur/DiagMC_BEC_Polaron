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


//Basic Vector operations

std::vector<double> operator+(const std::vector<double> & vec1, const std::vector<double> & vec2){
  try {
	if (vec1.size() != vec2.size()) {throw dnf_vecsize();}
	std::vector<double> tmp(vec1.size());
	for (int i =0; i<tmp.size(); i++){
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
	for (int i =0; i<tmp.size(); i++){
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
	for (int i =0; i<vec1.size(); i++){
	  tmp += pow(vec1[i],2);
	}
	return tmp;
  }
  catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}




//Definitions of Diagramm class

Diagram::Diagram() {
  order=0;
  
  times.reserve(10000);
  times.assign(2, std::vector<double>(2, 0));
  
  elprop.reserve(10000);
  elprop.assign(1, std::vector<double>(3, 0));
  
  phprop.reserve(100);
  phprop.assign(1, std::vector<double>(3, 0));
  
  pr_arc = 0;
  pr_tauin = 0;
  pr_taufin = 0;
  pr_tau1.assign(2, 0);
  pr_tau2.assign(2, 0);
  pr_q.assign(3, 0);
  pr_p.assign(3, 0);
}

Diagram::Diagram(const double & p, const double & tau, const double & my, const double & omega, const double & alph, const std::function<double()> & rnd): mu(my), wp(omega), alpha(alph), order(0), drnd(rnd) {
  times.reserve(50);
  times.assign(2, std::vector<double>(2, 0));
  times[0][1] = 1;
  times[1][0] = tau;
    
  elprop.reserve(50);
  elprop.assign(1, std::vector<double>(3, 0));
  elprop[0][0] = p;
  
  phprop.reserve(50);
  phprop.assign(1, std::vector<double>(3, 0));
  
  pr_arc = 0;
  pr_tauin = 0;
  pr_taufin = 0;
  pr_tau1.assign(2, 0);
  pr_tau2.assign(2, 0);
  pr_q.assign(3, 0);
  pr_p.assign(3, 0);
}

void Diagram::set(const double & p, const double & tau, const double & my, const double & omega, const double & alph, const std::function<double()> & rnd) {
  mu = my;
  wp = omega;
  alpha = alph;
  drnd = rnd;
  
  times[0][1] = 1;
  times[1][0] = tau;
    
  elprop[0][0] = p;
}

void Diagram::random_arc() {
  pr_arc = (int)(drnd()*(2*get_order()+1));
}


int Diagram::propose_insert() {
  try{
	random_arc();										//arc to insert
	//if (order != 0 && order !=1) {return -2;}						//step 2 just 1st order allowed
	
	pr_tauin = get_tinit(pr_arc);			
	pr_taufin = get_tfin(pr_arc);
	
	pr_tau1[0] = (drnd()*(pr_taufin - pr_tauin))+ pr_tauin;			//tau1
	pr_tau1[1] = pr_arc+2;
	pr_tau2[0] = (drnd()*(pr_taufin - pr_tau1[0]))+pr_tau1[0];							//tau2
	pr_tau2[1] = pr_arc+1;
	
	pr_q[0] = sqrt(-2*log(drnd()))*cos(2*M_PI*drnd())/sqrt(pr_tau2[0]-pr_tau1[0]);		//qx
	pr_q[1] = sqrt(-2*log(drnd()))*sin(2*M_PI*drnd())/sqrt(pr_tau2[0]-pr_tau1[0]);		//qy
	pr_q[2] = sqrt(-2*log(drnd()))*cos(2*M_PI*drnd())/sqrt(pr_tau2[0]-pr_tau1[0]);		//qz
	
	pr_p = get_p(pr_arc);
	 
	return 0;
  }	catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

int Diagram::propose_remove() {
  try{
	random_arc();	//arc to remove
	if (pr_arc != (int)(times[pr_arc+1][1]+0.5)) {
	  return -1;}
	if (get_order() == 0) {return -1;}
	
	pr_tauin = get_tinit(pr_arc-1);
	pr_taufin = get_tfin(pr_arc+1);
	
	pr_tau1[0] = get_tinit(pr_arc);			//tau1
	pr_tau1[1] = -1;
	pr_tau2[0] = get_tfin(pr_arc);							//tau2
	pr_tau2[1] = -1;
	
	pr_q = get_q(pr_arc);		
	
	pr_p = get_p(pr_arc-1);
	return 0;
  }	catch (std::exception& e) {
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

int Diagram::propose_swap() {
  random_arc();
  if (get_order() == 0) {return -1;}
  if (pr_arc == 0 || (int)(times[pr_arc+1][1]+0.5) == 0) {return -1;}
  if ((int)(times[pr_arc][1]+0.5) == (pr_arc+1)) {return -1;}
  
  pr_tauin = get_tinit(pr_arc);
  pr_taufin = get_tfin(pr_arc);
  
  pr_tau1.assign(2,0);
  pr_tau2.assign(2,0);
  pr_q.assign(3,0);	
  
  //std::cout << get_p(pr_arc) << '\t' << get_q(pr_arc) << '\t' << get_q(pr_arc-1) << '\t' << get_p(pr_arc) << '\t' <<
  
  pr_p = get_p(pr_arc)+ get_q(pr_arc) + get_q(pr_arc) - get_q(pr_arc-1) -get_q(pr_arc +1); //pnew
  
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

int Diagram::propose_ct() {
  if (get_order() == 0) {return -1;}				
  
  random_arc();
  
  
  
  pr_tau1[0] = get_tinit(pr_arc);
  pr_tau1[1] = double(pr_arc);
  pr_arc += 1;
  
  if ((int)(times[pr_arc][1]+0.5) == 0) {return 1;} // last tau
  
  pr_tau2[0] = get_tfin(pr_arc);
  pr_tau2[1] = double(pr_arc+1);
  
  
  pr_tauin = (drnd()*(pr_tau2[0] - pr_tau1[0]))+ pr_tau1[0];		//tau new
  pr_taufin = get_tinit(pr_arc);									//tau old
  
  pr_q.assign(3,0);
  
  pr_p = get_p(pr_arc-1);
  
  return 0; 
}



double Diagram::high_weight() {
  return (alpha*exp(((square(pr_p-pr_q)/2) - mu)*(pr_tau1[0]-pr_tau2[0]) - wp*(pr_tau2[0]-pr_tau1[0])) / square(pr_q));
}

double Diagram::low_weight() {
  return (exp(((square(pr_p)/2) - mu)*(pr_tau1[0]-pr_tau2[0])));
}

double Diagram::P_hilo(){
  return 1/(double(get_order())+1);
}

double Diagram::P_lohi(){
  return (1/(2*double(get_order())+1) /(pr_taufin-pr_tauin) /(pr_taufin-pr_tau1[0]) *pow((pr_tau2[0]-pr_tau1[0])/2/M_PI, 3/2) * exp((square(pr_q)/2) *(pr_tau1[0]-pr_tau2[0])));
}

double Diagram::G0el(const std::vector< double > & p) {
  return exp(-(square(p)/2 - mu)*(pr_taufin-pr_tauin));
}

double Diagram::Dph(const double & factor) {
  return exp((-wp)*factor*(pr_taufin-pr_tauin));
}

void Diagram::insert() {
  try{
    if (pr_tau1.empty() || pr_tau2.empty() || pr_q.empty()) {throw no_propose();};
	std::vector< std::vector<double> >::iterator tit = times.begin() + pr_arc;		//iterators for different matrices
	std::vector< std::vector<double> >::iterator pit = phprop.begin() + pr_arc;
	std::vector< std::vector<double> >::iterator eit = elprop.begin() + pr_arc;
	
	for (int i = 0; i<((2*get_order())+1); i++) {
	  if ((int)(times[i][1]+0.5) > pr_arc) {times[i][1] +=2;}
	}
	times.insert(tit+1, pr_tau1);										//insert tau1 and tau2
	times.insert(tit+2, pr_tau2);
	
	
	phprop.insert(pit+1, std::move(pr_q + get_q(pr_arc)));							//insert q+oldq in arc+1
	phprop.insert(pit+2, std::move(get_q(pr_arc)));								//insert oldq in arc+2
	
	elprop.insert(eit+1, std::move(get_p(pr_arc) - pr_q));					//insert oldp-q in arc+1
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
	  
	for (int i = 0; i<((2*get_order())+1); i++) {
	  if ((int)(times[i][1]+0.5) > pr_arc) {times[i][1] -=2;}
	}
	
  } catch (std::exception& e) {
      std::cerr << e.what() << std::endl;
      exit(EXIT_FAILURE);
  }

}

int Diagram::set_tau(double tau) {
  if (tau < times[2*get_order()][0]) {return -1;}
  times[2*get_order()+1][0] = tau;
  
  return 0;
}

void Diagram::swap() {
  double tmp = times[(int)(times[pr_arc][1]+0.5)][1];
  times[(int)(times[pr_arc][1]+0.5)][1] = times[(int)(times[pr_arc+1][1]+0.5)][1];
  times[(int)(times[pr_arc+1][1]+0.5)][1] = tmp;
  
  tmp = times[pr_arc][1];			//swap links
  times[pr_arc][1] = times[pr_arc+1][1];
  times[pr_arc+1][1] = tmp;
  
  
  
  phprop[pr_arc] = (get_q(pr_arc-1) + get_q(pr_arc+1) - get_q(pr_arc));
  
  elprop[pr_arc] = pr_p;

}

double Diagram::ct() {
  times[pr_arc][0]= pr_taufin;
  
  return pr_taufin;
}








//Errors for Testing

class dnf_diagvec: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Diagram vector sizes do not fit!";
  }
};

class timeserr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Error in Diagram::times!";
  }
};

class phproperr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Error in Diagram::phprop!";
  }
};

class elproperr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Error in Diagram::elprop!";
  }
};

class momerr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Momentum is not conserved!";
  }
};

class propopen: public std::exception
{
  virtual const char* what() const throw()
  {
    return "A propergator is not closed";
  }
};



void Diagram::test() {
  try {
	printall();
	if (times.size() != (2*get_order())+2) {throw dnf_vecsize();}
	if (phprop.size() != times.size()-1) {throw dnf_vecsize();}
	if (elprop.size() != phprop.size()) {throw dnf_vecsize();}
	
	if ((int)(times[0][1]+0.5) != (2*get_order())+1) {throw timeserr();}
	if ((int)(times[(2*get_order())+1][1]+0.5) != 0) {throw timeserr();}
	
	for (int i = 0; i< 3; i++) {
	  if (fabs(get_q(0).at(i) - get_q(2*get_order()).at(i)) > 0.0000001) {throw phproperr();}
	  if (fabs(get_p(0).at(i) - get_p(2*get_order()).at(i)) > 0.0000001) {throw elproperr();}
	}
	
	
	  
	  
	VectorXd tmp=VectorXd::Zero((2*get_order())+2);
	for (int i=0; i<(2*get_order())+2; i++) {
	  tmp((int)(times[i][1]+0.5)) +=1;
	  if ((int)(times[i][1]+0.5) == i) {throw timeserr();}
	}
	for (int i=0; i<(2*get_order())+1; i++) {
	  if (times[i] > times[i+1]) {throw timeserr();}
	  for (int i2 = 0; i2< 3; i2++) {
		if (fabs((get_p(i).at(i2)+get_q(i).at(i2)) - get_p(0).at(i2)) > 0.0000001)  {throw momerr();}

	  }
	}
	
	for (int i=0; i<(2*get_order())+2; i++) {
	  if ((int)(tmp(i)+0.5) != 1) {throw propopen();}
	}
  } catch (std::exception& e) {
      std::cerr << e.what() << std::endl;
      exit(EXIT_FAILURE);
  }

}

void Diagram::printall() {
  std::cout << mu <<'\t' << wp << '\t' << alpha <<'\t' << order << '\t' << drnd()  << std::endl;
  std::cout <<'\n';
  for (int i=0 ;  i<times.size() ; i++){
	for (int i2=0 ;  i2<times[0].size() ; i2++) {
	  std::cout << times[i][i2] <<'\t';
	}
	std::cout <<'\n';
  }
  std::cout <<'\n';
  for (int i=0 ;  i<phprop.size() ; i++){
	for (int i2=0 ;  i2<phprop[0].size() ; i2++) {
	  std::cout << phprop[i][i2] <<'\t';
	}
	std::cout <<'\n';
  }
  std::cout <<'\n';
  for (int i=0 ;  i<elprop.size() ; i++){
	for (int i2=0 ;  i2<elprop[0].size() ; i2++) {
	  std::cout << elprop[i][i2] <<'\t';
	}
	std::cout <<'\n';
  }
  std::cout <<'\n';
  
  std::cout << pr_arc <<'\t' << pr_tauin << '\t' << pr_taufin  << std::endl;
  std::cout <<'\n';
  for (int i=0 ;  i<pr_tau1.size() ; i++){
	std::cout << pr_tau1[i] <<'\t';
	std::cout <<'\n';
  }
  std::cout <<'\n';
  for (int i=0 ;  i<pr_tau2.size() ; i++){
	std::cout << pr_tau2[i] <<'\t';
	std::cout <<'\n';
  }
  std::cout <<'\n';
  for (int i=0 ;  i<pr_q.size() ; i++){
	std::cout << pr_q[i] <<'\t';
	std::cout <<'\n';
  }
  std::cout <<'\n';
  for (int i=0 ;  i<pr_p.size() ; i++){
	std::cout << pr_p[i] <<'\t';
	std::cout <<'\n';
  }
  std::cout <<std::endl;
  

}



