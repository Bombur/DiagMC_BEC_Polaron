#include "DiagMC.h"

namespace pt = boost::property_tree;

class oor_Probs: public std::exception {
  virtual const char* what() const throw()
  {
    return "The Probabilities do not add up to 1!";
  }
};  

class openwritefile: public std::exception {
  virtual const char* what() const throw()
  {
    return "Write File not open!";
  }
};

DiagMC::DiagMC(const int & seed, const pt::ptree & config):p(config.get<double>("Momentum")), mu(config.get<double>("Chemical_Potential")), taumax(config.get<double>("Tau_max")), taubin(config.get<int>("Tau_bin")),alpha(2*sqrt(2)*M_PI*config.get<double>("Coupling_Strength")), omegap(config.get<double>("Omega_Phonon")), maxord(config.get<double>("Max_Order")),  
										ctcor(config.get<double>("Correction_tau")),  qcor(config.get<double>("Correction_dq")),  dtins(config.get<double>("Insert_at_End_DT")), dqins(config.get<double>("Insert_at_End_DQ")),
										Prem(config.get<double>("Remove_Probability")), Pins(config.get<double>("Insert_Probability")), Pct(config.get<double>("Change_tau_Probability")), Pctho(config.get<double>("Change_tau_in_HO_Probability")), Psw(config.get<double>("Swap_Probability")), Pdq(config.get<double>("DQ_Probability")), Piae(config.get<double>("IAE_Probability")), Prae(config.get<double>("RAE_Probability")),
										Meas_its(config.get<int>("Its_per_Measure")), Test_its(config.get<int>("Its_per_Test")), Write_its(config.get<int>("Its_per_Write")), RunTime(config.get<int>("RunTime")) { 
  try{
	if (fabs(1-Prem-Pins-Pct -Pctho -Piae -Prae - Pdq-Psw) > 0.0000001) {throw oor_Probs();}
	
	E = pow(p,2.)/2. - mu;
	G0p = (1.-exp(-E*taumax))/E;
  
	//Random Number Generator
	std::mt19937 generator (seed);
	std::uniform_real_distribution<double> uni_dist2(0,1);
	drnd = std::bind ( uni_dist2, generator);
	
	Data = MatrixXi::Zero(taubin, 4);
	
	diag.set(p, config.get<double>("Tau_start"), drnd);
  
	updatestat = MatrixXd::Zero(11,6);
	orderstat = VectorXi::Zero(40);
  
	std::stringstream convert, convert2, convert3, convert4, convert5, convert6, convert7; //p, mu, taumax, alpha, wp, RunTime
	convert<<p;
	convert2<<mu;
	convert3<<taumax;
	convert4<<alpha;
	convert5<<omegap;
	convert6<<RunTime;
	convert7<<seed;
	path=config.get<std::string>("Path");
	
	//Datafile.open(path+"data/Data_p_"+convert.str()+"_mu_"+convert2.str()+"_tmax_"+convert3.str()+"_a_"+convert4.str()+"_wp_"+convert5.str()+"_Run_Time_"+convert6.str()+"_core_"+convert7.str()+".txt");
	//udsfile.open(path+"data/udstat_p_"+convert.str()+"_mu_"+convert2.str()+"_tmax_"+convert3.str()+"_a_"+convert4.str()+"_wp_"+convert5.str()+"_Run_Time_"+convert6.str()+"_core_"+convert7.str()+".txt");
	//osfile.open(path+"data/ostat_p_"+convert.str()+"_mu_"+convert2.str()+"_tmax_"+convert3.str()+"_a_"+convert4.str()+"_wp_"+convert5.str()+"_Run_Time_"+convert6.str()+"_core_"+convert7.str()+".txt");
	//tsfile.open(path+"data/tstat_p_"+convert.str()+"_mu_"+convert2.str()+"_tmax_"+convert3.str()+"_a_"+convert4.str()+"_wp_"+convert5.str()+"_Run_Time_"+convert6.str()+"_core_"+convert7.str()+".txt");
	//if (!(Datafile.is_open() && udsfile.is_open() && osfile.is_open() && tsfile.is_open())) {throw openwritefile();}
	//Datafile_pos=Datafile.tellp();
	//uds_pos=udsfile.tellp();
	//os_pos=udsfile.tellp();
	//ts_pos=udsfile.tellp();
	
	global_weight=exp(-E*config.get<double>("Tau_start"));
	
	lu = "nothing";
  }
  catch (std::exception& e){
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

DiagMC::~DiagMC() {
  //Datafile.close();
  //udsfile.close();
  //osfile.close();
  //tsfile.close();
}

