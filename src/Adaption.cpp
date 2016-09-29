#include "Adaption.h"

namespace pt = boost::property_tree;

Adaption::Adaption(const int & oss, const double & time, const double & normal, const double & ende, const pt::ptree & config) : tls(time),nnorms(normal),nends(ende),normmin(config.get<double>("Norm_Points")),endmin(config.get<double>("End_Points")),texp(config.get<double>("RunTime")){
  ordstsz = oss;
  enddev = (nends/endmin) -1.;
  normdev = (nnorms/normmin) -1.;
  countdev = (normdev < enddev) ? normdev : enddev;
  tdev = tls/texp -1.;

}

bottleneck Adaption::whichbtlnk() {
  try{
	const std::array<double, 3> compare = {normdev, enddev, tdev};
	std::size_t max = std::min_element(compare.begin(), compare.end()) - compare.begin();
	if ( max > 2) {throw std::out_of_range("Adaption compare failure!");}
	bottleneck tmp = static_cast<bottleneck>(max);
	return tmp;
  } catch (std::exception& e){
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

int Adaption::ordstszadapt() {
  bottleneck which = whichbtlnk();
  switch(which) {
	case NORM: ordstsz = (int)((ordstsz/(1.+tdev))+0.5); break;
	case END: ordstsz = (int)((ordstsz/(1.+tdev))+0.5); break;
	case TIME: ordstsz = (int)((ordstsz*(1.+countdev))+0.5); break;
  	default: std::cerr << "Default chosen in Adaption!" << std::endl; break;
  }
  if (ordstsz == 0) {ordstsz =1;}
  return ordstsz;
}

void Adaption::printall() {
  std::cout<< ordstsz << '\n' << nnorms <<'\n'<<  nends << '\n' << tls <<"\n\n";
  
  std::cout << normmin <<'\n' << endmin <<'\n'<< texp<< "\n\n";
	
  std::cout << countdev <<'\n' << normdev <<'\n' << enddev <<'\n' << tdev << "\n\n";
  
  std::cout  << static_cast<int>(whichbtlnk()) << '\n' << std::endl;
  
}
