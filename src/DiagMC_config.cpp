#include "DiagMC.h"

namespace pt = boost::property_tree;

DiagMC::DiagMC(const pt::ptree & config):p(config.get<double>("Momentum")), mu(config.get<double>("Chemical_Potential")), taumax(config.get<double>("Tau_max")), taubin(config.get<double>("Tau_bin")), Meas_its(config.get<int>("Its_per_Measure")), Test_its(config.get<int>("Its_per_Test")), Write_its(config.get<int>("Its_per_Write")), RunTime(config.get<int>("RunTime")) { 
  try{
	E = pow(p,2)/2 - mu;
	G0p = (1-exp(-E*taumax))/E;
	tau = config.get<double>("Tau_start");
  
	alpha = config.get<double>("Coupling_Strength");
   
	//Seed for Random Number Generator
	unsigned int seed = 1;

	//Random Number Generator
	std::mt19937 generator (seed);
	std::uniform_real_distribution<double> uni_dist2(0,1);
	drnd = std::bind ( uni_dist2, generator);
	
	Data = VectorXd::Zero(taubin);
  
	binl=config.get<int>("Binning_Length");
	anabuffer.setZero();
	stats.setZero();
  
	std::stringstream convert, convert2;
	convert<<p;
	convert2<<RunTime;
	path=config.get<std::string>("Path");
	Datafile.open(path+"data/Data_p_"+convert.str()+"_Run_Time_"+convert2.str()+".txt");
	Statsfile.open(path+"data/Statistics_p_"+convert.str()+"_Run_Time_"+convert2.str()+".txt");
	if (!(Datafile.is_open() && Statsfile.is_open())) {throw openwritefile();}
	Datafile_pos=Datafile.tellp();
	Statfile_pos=Statsfile.tellp();
  }
  catch (std::exception& e){
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

DiagMC::~DiagMC() {
  Datafile.close();
  Statsfile.close();  
}

