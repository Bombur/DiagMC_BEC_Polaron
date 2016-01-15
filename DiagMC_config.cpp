#include "DiagMC.h"

namespace pt = boost::property_tree;

DiagMC::DiagMC(const int & seed, const pt::ptree & config):p(config.get<double>("Momentum")), qc(config.get<double>("Q_Cutoff")), mu(config.get<double>("Chemical_Potential")), alpha(config.get<double>("Alpha")), relm(config.get<double>("Impurity_Mass")), maxord(config.get<int>("Max_Order")),  
										ctcor(config.get<double>("Correction_tau")),  qcor(config.get<double>("Correction_dq")),  dtins(config.get<double>("Insert_FP_DT")), dqins(config.get<double>("Insert_FP_DQ")), fw(config.get<double>("Fake_Weight")), sigfac(config.get<double>("QSigma_Factor")),
										Prem(config.get<double>("Remove_Probability")), Pins(config.get<double>("Insert_Probability")), Pct(config.get<double>("Change_tau_Probability")), Pctho(config.get<double>("Change_tau_in_HO_Probability")), Psw(config.get<double>("Swap_Probability")), Pdq(config.get<double>("DQ_Probability")), Piae(config.get<double>("IAE_Probability")), Prae(config.get<double>("RAE_Probability")),
										Meas_its(config.get<int>("Its_per_Measure")), Test_its(config.get<int>("Its_per_Test")), Write_its(config.get<int>("Its_per_Write")), RunTime(config.get<int>("RunTime")),
										ordstsz(config.get<int>("Order_Step_Size")), normmin(config.get<int>("Norm_Points")), endmin(config.get<int>("End_Points")), TotRunTime(config.get<int>("Total_Time")), TotMaxOrd(config.get<int>("Total_Max_Order")),
										taumap(tmap(create_fvec(config, "Functions"), as_vector<int>(config, "Bins"), as_vector<double>(config, "Taus")))
	
{ 
  try{
	if (fabs(1-Prem-Pins-Pct -Pctho -Piae -Prae - Pdq-Psw) > 0.0000000001) {throw oor_Probs();}

	E = pow(p,2.)/2. - mu;
	G0p = (1.-exp(-E*taumap.taumax))/E;
	
#ifdef FP
	wp = config.get<int>("FP_Omega");
#endif
  
	//Random Number Generator
	std::mt19937 generator (seed);
	std::uniform_real_distribution<double> uni_dist2(0,1);
	drnd = std::bind ( uni_dist2, generator);
	
	Data = ArrayXXd::Zero(taumap.taubin, 5);
	
	diag.set(p, config.get<double>("Tau_start"), drnd);
	
	
	//Cumulative SE Calculation
	minord = 0;
	ordstep = 0;
#ifdef SECUMUL
	if ((TotMaxOrd < ordstsz) && (TotMaxOrd > 0)) {
	  maxord = TotMaxOrd;
	} else {
	  maxord = minord + ordstsz;
	}
	SEib= ArrayXd::Zero(taumap.taubin);
	Norms = ArrayXd::Zero(taumap.taubin);
	Ends = ArrayXd::Zero(taumap.taubin);
	nnorms.assign(1, 0.);
	nends.assign(1, 0.);
#endif
  
	updatestat = ArrayXXd::Zero(11,6);
	orderstat = ArrayXi::Zero(40);
  
	std::stringstream convert, convert2, convert3, convert4, convert5, convert6, convert7; //p, mu, taumap.taumax, alpha, wp, RunTime
	convert<<p;
	convert2<<mu;
	convert3<<taumap.taumax;
	convert4<<alpha;
	convert5<<relm;
	convert6<<RunTime;
	convert7<<seed;
	path=config.get<std::string>("Path");
	
	global_weight=exp(-E*config.get<double>("Tau_start"));
	
	lu = "nothing";
  }
  catch (std::exception& e){
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

DiagMC::~DiagMC() {

}

