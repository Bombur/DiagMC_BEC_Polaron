#include "DiagMC.h"

namespace pt = boost::property_tree;

DiagMC::DiagMC(const int & thread, const pt::ptree & config):p(config.get<double>("Momentum")), qc(config.get<double>("Q_Cutoff")), mu(config.get<double>("Chemical_Potential")), alpha(config.get<double>("Alpha")), relm(config.get<double>("Impurity_Mass")),
										wmax(config.get<double>("Omega_max")), wbin(config.get<int>("Omega_bin")), qsigma(config.get<int>("Q_Gauss")),
										maxord(config.get<int>("Max_Order")),  
										ctcor(config.get<double>("Correction_tau")),  qcor(config.get<double>("Correction_dq")),  dtins(config.get<double>("Insert_FP_DT")), dqins(config.get<double>("Insert_QRange")), fw(config.get<double>("Fake_Weight")), fwzero(config.get<double>("Fake_Weight_Zero")), fwone(config.get<double>("Fake_Weight_One")), sigfac(config.get<double>("QSigma_Factor")),
										fog0seg0(config.get<bool>("FOG0SEG0")), ct_taucut(config.get<bool>("CT_Tau_Cut")),ins_taucut(config.get<bool>("Ins_Tau2_Cut")),ins_tau2lin(config.get<bool>("Ins_Tau2_Linear")),ic_taucut(config.get<bool>("IC_Tau_Cut")),
										Prem(config.get<double>("Remove_Probability")), Pins(config.get<double>("Insert_Probability")), Pct(config.get<double>("Change_tau_Probability")), Pctho(config.get<double>("Change_tau_in_HO_Probability")), Psw(config.get<double>("Swap_Probability")), Pdq(config.get<double>("DQ_Probability")), Pic(config.get<double>("IC_Probability")), Prc(config.get<double>("RC_Probability")),
										Meas_its(config.get<int>("Its_per_Measure")), Test_its(config.get<int>("Its_per_Test")), Write_its(config.get<int>("Its_per_Write")), RunTime(config.get<int>("RunTime")),
										ordstsz(config.get<int>("Order_Step_Size")), normmin(config.get<int>("Norm_Points")), endmin(config.get<int>("End_Points")), TotRunTime(config.get<int>("Total_Time")), TotMaxOrd(config.get<int>("Total_Max_Order")),
										taumap(tmap(create_fvec(config, "Functions"), as_vector<int>(config, "Bins"), as_vector<double>(config, "Taus")))
{ 
  try{
	if (fabs(1-Prem-Pins-Pct -Pctho -Pic -Prc - Pdq-Psw) > 0.0000000001) {throw oor_Probs();}

#ifdef FP
	mass = 1.;
	wp = config.get<int>("FP_Omega");
#endif
	
#ifdef BEC
	mass = relm /sqrt(2.);
#endif
	
	E = pow(p,2.)/2./mass - mu;
	G0p = (1.-expfun(-E*taumap.taumax))/E;
  
	//Random Number Generator
	pcg32 generator(config.get<pcg32::state_type>("Seed"), thread);
	//std::mt19937 generator (thread);
	std::uniform_real_distribution<double> uni_dist2(0,1);
	drnd = std::bind ( uni_dist2, generator);
	
	Data = ArrayXXd::Zero(taumap.taubin, 5);
	
	diag.set(p, config.get<double>("Tau_start"), drnd);
	
	//Estimator
	std::vector<double> wstmp = as_vector<double>(config, "Ws_for_Epol");
	ws = Map<const Array<double, 1, Dynamic> > (wstmp.data(), wstmp.size());
	Epol = ArrayXXd::Zero(taumap.taubin, ws.size());
	SE = ArrayXXd::Zero(taumap.taubin, 3);
	G0SEiw = ArrayXXcd::Zero(wbin, 2);

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
	orderstat = ArrayXd::Zero(config.get<int>("Order_Stat_Size"));
	qstat = ArrayXXd::Zero(config.get<int>("QStat_Size"),2);
	tstat = ArrayXXd::Zero(taumap.taubin,4);
	testhisto = ArrayXXd::Zero(2*taumap.taubin+2, 2);
  
	path=config.get<std::string>("Path");
	
	global_weight=expfun(-E*config.get<double>("Tau_start"));
	
	lu = "nothing";
  }
  catch (std::exception& e){
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

DiagMC::~DiagMC() {

}

