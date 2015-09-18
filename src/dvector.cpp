#include "dvector.h"

std::vector<double> vadd(const std::vector<double> & vec1, const std::vector<double> & vec2){
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

std::vector<double> vsub(const std::vector<double> & vec1, const std::vector<double> & vec2){
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

double vsq(const std::vector<double> & vec1){
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