#include "DiagMC.h"

namespace pt = boost::property_tree;

DiagMC::DiagMC(const int & thread, const pt::ptree & config):p(config.get<double>("Momentum")), qc(config.get<double>("Q_Cutoff")), mu(config.get<double>("Chemical_Potential")), alpha(config.get<double>("Alpha")), relm(config.get<double>("Impurity_Mass")),
										wmax(config.get<double>("Omega_max")), wbin(config.get<int>("Omega_bin")), 
										maxord(config.get<int>("Max_Order")),  
										ctcor(config.get<double>("Correction_tau")),  qcor(config.get<double>("Correction_dq")),  dtins(config.get<double>("Insert_FP_DT")), dqins(config.get<double>("Insert_QRange")), 
										fw(config.get<double>("Fake_Weight")),
										qsigma(config.get<int>("Q_Gauss")), explim(config.get<double>("Exp_Arg_Limit")), doitlim(config.get<int>("Do_While_Limit")),
										fo_insrem(config.get<bool>("FO_InsRem")), fog0seg0(config.get<bool>("FOG0SEG0")),
										ct_taucut(config.get<bool>("CT_Tau_Cut")),ct_lin(config.get<bool>("CT_Tau_Linear")),ins_tau2cut(config.get<bool>("Ins_Tau2_Cut")),ins_tau2lin(config.get<bool>("Ins_Tau2_Linear")),ic_taulin(config.get<bool>("IC_Tau2_Linear")),ic_taucut(config.get<bool>("IC_Tau2_Cut")),
										Prem(config.get<double>("Remove_Probability")), Pins(config.get<double>("Insert_Probability")), Pct(config.get<double>("Change_tau_Probability")), Pctho(config.get<double>("Change_tau_in_HO_Probability")), Psw(config.get<double>("Swap_Probability")), Pdq(config.get<double>("DQ_Probability")), Pic(config.get<double>("IC_Probability")), Prc(config.get<double>("RC_Probability")),
										Meas_its(config.get<int>("Its_per_Measure")), Test_its(config.get<int>("Its_per_Test")), Write_its(config.get<int>("Its_per_Write")), RunTime(config.get<int>("RunTime")),
										ordstsz(config.get<int>("Order_Step_Size")), normmin(config.get<int>("Norm_Points")/config.get<int>("NCores")), endmin(config.get<int>("End_Points")/config.get<int>("NCores")), TotRunTime(config.get<int>("Total_Time")), TotMaxOrd(config.get<int>("Total_Max_Order")),
										taumap(tmap(create_fvec(config), as_vector<int>(config, "Bins"), as_vector<double>(config, "Taus"), config.get<int>("Tminit"), config.get<int>("Tmaxit")))
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
	std::vector<double> taus = as_vector<double>(config, "Taus");
	G0p1arm = (1.-expfun(-E*(taumap.taug0max - taumap.taumin)))/E;
	G0p2arms = (G0p1arm/E) - (expfun(-E*(taumap.taug0max - taumap.taumin))/E*(taumap.taug0max - taumap.taumin));
	G0ptmima = expfun(-E*taumap.taumin)*(1.-expfun(-E*(taumap.taug0max-taumap.taumin)))/E;
  
	//Random Number Generator
	pcg32 generator(config.get<pcg32::state_type>("Seed"), thread);
	//std::mt19937 generator (thread);
	std::uniform_real_distribution<double> uni_dist2(0,1);
	drnd = std::bind ( uni_dist2, generator);
	
	Data = ArrayXXd::Zero(taumap.taubin, 5);
	
	double tstart = (config.get<double>("Tau_start")>taumap.taumin)? config.get<double>("Tau_start"): (((taumap.taumax-taumap.taumin)/10) + taumap.taumin);
	diag.set(p, tstart, drnd);
	
	//Estimator
	std::vector<double> wstmp = as_vector<double>(config, "Ws_for_Epol");
	ws = Map<const Array<double, 1, Dynamic> > (wstmp.data(), wstmp.size());
	Epol = ArrayXXd::Zero(taumap.taubin, ws.size());
	SE = ArrayXXd::Zero(taumap.taubin, 3);
	G0SEiw = ArrayXXcd::Zero(wbin, 2);

	//Fake_Weight_vector
	int start = config.get<int>("Fake_Weight_to_ord.start");
	int end = config.get<int>("Fake_Weight_to_ord.end");
	double fwtmp =config.get<double>("Fake_Weight_to_ord.fw");
	fwtab.assign(end, 1.);
	fwtab[0] = (1./config.get<double>("Fake_Weight_Zero"));
	fwtab[1] = (1./config.get<double>("Fake_Weight_One"));
	for (int i = start; i< fwtab.size(); i++){fwtab[i] *= (1./fwtmp);}
		
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
	
	fw_last.fill(0);
	fw_counts.fill(0.);
	fw_max = config.get<double>("Fake_Weight_Max");
	desi_ordrat= config.get<double>("Desired_Min_Max_Ratio");
#endif
  
	updatestat = ArrayXXd::Zero(13,6);
	overflowstat = ArrayXXd::Zero(13,6);
	orderstat = ArrayXd::Zero(config.get<int>("Order_Stat_Size"));
	qstat = ArrayXXd::Zero(config.get<int>("QStat_Size"),2);
	tstat = ArrayXXd::Zero(taumap.taubin,4);
	testhisto = ArrayXXd::Zero(taumap.taubin+2, 4);
  
	path=config.get<std::string>("Path");
	
	global_weight= logG0el(diag.get_p(0), diag.get_tfin(0), diag.get_tinit(0));

	lu = "nothing";
  } 
  catch (std::exception& e){
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

DiagMC::~DiagMC() {

}

