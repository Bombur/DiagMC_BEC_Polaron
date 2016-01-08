#include "tmap.h"

tmap::tmap(const std::vector<std::function<double(int)> > & fvec, const std::vector<int> bins, const std::vector<double> taus):taubin(bins[bins.size()-1]) {
  try{
	//Tests on Arguments
	std::size_t sizetest=fvec.size();
	if ((bins.size()!= sizetest+1) || (taus.size()!= sizetest+1)) {throw std::out_of_range("Input vector sizes do not fit!");}
	if (bins[0] != 0) {throw std::invalid_argument("bins[0] != 0");}
	for (std::size_t it = 0; it< fvec.size(); it++) {
	  if ((bins[it] > bins[it+1]) || (taus[it] > taus[it+1])) {throw std::invalid_argument("Bins or Taus vectors are not monoton increasing!");}
	}
  
	//creating total vector
	std::vector<double> out;
	int pos= 0; //current position in out
	out.reserve(bins[bins.size()-1]+100);
	for (std::size_t it = 0; it< fvec.size(); it++) {
	  std::vector<double> tmp = discretize(fvec[it], bins[it], bins[it+1]);
	  for (const auto & i : tmp) {
		out.push_back(i);
	  }
	  if (fabs(taus[it] - out[pos]) > 0.01 * taus[it]) {throw std::range_error("Taumin and created vector min do not fit!");}
	  pos = bins[it+1];
	  if (taus[it+1] < out[pos]) {throw std::range_error("Taumax is smaller than created vector max!");}
	}

	//Test on Vector
	for (std::size_t it = 0; it< out.size()-1; it++) {
	  std::cout << out[it] << '\t' <<out[it+1] <<std::endl;
	  if (out[it] > out[it+1]) {throw std::range_error("Created Vector is not monoton increasing!");}
	}
	
	//converting to map
	int val_bin = 0; // Value bin for map
	for (const auto & i : out) {
	  mymap[i] = val_bin;
	  val_bin += 1;
	}
	mymap[taus[taus.size()-1]] = bins[bins.size()-1]; 
	
	
		 print_all();
	//Tests on converted map
	if (fabs(mymap.cbegin()->first - 0.) > 0.000000001 || mymap.cbegin()->second != 0) {throw std::logic_error("tau[0] != 0!");}
	

	for (auto mapit = mymap.cbegin(); mapit!= std::prev(mymap.cend(), 1); ++mapit){
	  if (mapit->first > std::next(mapit, 1)-> first) {throw std::range_error("Map Keys are not monoton increasing!");}
	  if (mapit->second > std::next(mapit, 1)-> second) {throw std::range_error("Map Values are not monoton increasing!");}
	}
  
		
  } catch (std::exception& e) {
 	std::cerr << e.what() << std::endl;
 	exit(EXIT_FAILURE);
  }
}

std::vector<double> tmap::discretize(const std::function<double(int)> & f, const int & binmin , const int & binmax) {
  std::vector<double> out(binmax - binmin, 0.);
  for (int i = 0; i<(binmax-binmin); i++){
	out[i] = f(binmin + i);
  }
  return out;  
}

int tmap::bin(const double & tau) {
  std::map<double, int>::iterator upper;
  upper = mymap.upper_bound(tau);
  
#ifndef NDEBUG
  std::map<double, int>::iterator lower;
  lower = mymap.lower_bound(tau);

  assert(lower->second == upper->second);
#endif
  
  return std::prev(upper,1)->second;  
}

ArrayXd tmap::print() {
  ArrayXd output(mymap.size()-1);
  std::map<double, int>::iterator mapit = mymap.begin();
  for (int i = 0 ; i < (mymap.size()-1); i++) {
	output(i)= (std::next(mapit, i)->first + std::next(mapit, i+1)->first)/2;
  }
  return output;  
}


void tmap::print_all() {
  for (const auto & mapit : mymap) {
	std::cout << mapit.first << '\t' << mapit.second << '\n';
  }
  std::cout<< std::endl;
}