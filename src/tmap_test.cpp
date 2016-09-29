#include "tmap.h"
#include "dvector.h"
//config
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace pt = boost::property_tree;

int main() {

pt::ptree config;
pt::read_json("DiagMC_BEC.json", config);

tmap testmap(create_fvec(config, "Functions"), as_vector<int>(config, "Bins"), as_vector<double>(config, "Taus"));

std::cout <<testmap.bin(9.8e-4)<< "\n\n" << testmap.print()<< std::endl;

return 0;
}
