#include "Ising.h"
#include "NeighborTable.h"

namespace pt = boost::property_tree;

DiagMC::DiagMC(const pt::ptree & config) { 
  Nbr = NeighborTable(config.get<int>("Dimension")).get_Table();
  N = Nbr.rows();
  E=0;
  beta=0;
  sigma = VectorXi::Ones(N);
  m=sigma.dot(VectorXi::Ones(N));
}


Ising_markov_method::Ising_markov_method(const pt::ptree & config, double be) : Ising2d(config), Meas_its(config.get<int>("Its_per_Measure")), 
								Test_its(config.get<int>("Its_per_Test")), Write_its(config.get<int>("Its_per_Write")), RunTime(config.get<int>("RunTime")) {
  try {
	beta=be;
	
	//Seed for Random Number Generator
	unsigned int seed = 1;

	//Random Number Generator
	std::mt19937 generator (seed);
	std::uniform_int_distribution<int> uni_dist1(1,N);
	irnd = std::bind ( uni_dist1, generator);
	std::uniform_real_distribution<double> uni_dist2(0,1);
	drnd = std::bind ( uni_dist2, generator);
  
	Data = MatrixX2i::Zero(Write_its/Meas_its, 2);
  
	binl=config.get<int>("Binning_Length");
	meanbuffer.setZero();
	stats.setZero();
  
	std::stringstream convert, convert2;
	convert<<beta;
	convert2<<RunTime;
	path=config.get<std::string>("Path");
	Datafile.open(path+"data/Data_beta_"+convert.str()+"_Run_Time_"+convert2.str()+".txt");
	Statsfile.open(path+"data/Statisics_beta_"+convert.str()+"_Run_Time_"+convert2.str()+".txt");
	if (!(Datafile.is_open() && Statsfile.is_open())) {throw openwritefile();}
	Datafile_pos=Datafile.tellp();
	Statfile_pos=Statsfile.tellp();
  }
  catch (std::exception& e){
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}
