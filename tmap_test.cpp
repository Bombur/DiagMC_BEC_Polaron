#include "tmap.h"
#include "dvector.h"
//config
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace pt = boost::property_tree;

int main() {

pt::ptree config;
pt::read_json("DiagMC_BEC.json", config);

std::vector<int> bins = {0, 100, 200, 250, 300};
std::vector<double> taus = {0, 2.5, 13.6825, 14.8227, 15.1727};

std::vector<std::function<double(int)> > fvec = create_fvec(config, "Functions");

for (int i= 0; i<4; i++){
  std::cout << fvec[i](bins[i+1]) << '\n';
}
std::cout<<std::endl;

tmap testmap(create_fvec(config, "Functions"), bins, taus);

  
//testmap.print_all();

std::cout <<testmap.bin(14.9)<< "\n\n" << testmap.print()<< std::endl;

bins = as_vector<int>(config, "Bins");
taus = as_vector<double>(config, "Taus");
std::cout<<taus <<'\t' <<bins <<std::endl;
return 0;
}
