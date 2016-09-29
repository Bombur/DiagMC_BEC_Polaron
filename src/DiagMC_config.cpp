#include "DiagMC.h"

namespace pt = boost::property_tree;

DiagMC::DiagMC(const int & thread, const pt::ptree & config):p(config.get<double>("Momentum")), qc(config.get<double>("Q_Cutoff")), mu(config.get<double>("Chemical_Potential")), alpha(config.get<double>("Alpha")), relm(config.get<double>("Impurity_Mass")),
										wmax(config.get<double>("Omega_max")), wbin(config.get<int>("Omega_bin")), 
										maxord(config.get<int>("Max_Order")),  
										SE_eigen(config.get<bool>("SE_Eigen_Array")), Ep_bin_each_step(config.get<bool>("Ep_bin_in_each_step")), G0SEiw_meas(config.get<bool>("G0SEiw_meas")),
										ctcor(config.get<double>("Correction_tau")),  qcor(config.get<double>("Correction_dq")),  dtins(config.get<double>("Insert_FP_DT")), dqins(config.get<double>("Insert_QRange")), 
										fw(config.get<double>("Fake_Weight")),
										qsigma(config.get<int>("Q_Gauss")), explim(config.get<double>("Exp_Arg_Limit")), doitlim(config.get<int>("Do_While_Limit")),
										fo_insrem(config.get<bool>("FO_InsRem")), fog0seg0(config.get<bool>("FOG0SEG0")),
										ct_taucut(config.get<bool>("CT_Tau_Cut")),ct_lin(config.get<bool>("CT_Tau_Linear")),ins_tau2cut(config.get<bool>("Ins_Tau2_Cut")),ins_tau2lin(config.get<bool>("Ins_Tau2_Linear")),ic_taulin(config.get<bool>("IC_Tau2_Linear")),ic_taucut(config.get<bool>("IC_Tau2_Cut")),
										fst_Ep_meas(config.get<bool>("1st_Ep_meas")),
										Prem(config.get<double>("Remove_Probability")), Pins(config.get<double>("Insert_Probability")), Pct(config.get<double>("Change_tau_Probability")), Pctho(config.get<double>("Change_tau_in_HO_Probability")), Psw(config.get<double>("Swap_Probability")), Pdq(config.get<double>("DQ_Probability")), Pic(config.get<double>("IC_Probability")), Prc(config.get<double>("RC_Probability")),
										Meas_its(config.get<int>("Its_per_Measure")), Test_its(config.get<int>("Its_per_Test")), Write_its(config.get<int>("Its_per_Write")), Ep_meas_it(config.get<int>("Ep_Meas")),
										RunTime(config.get<int>("RunTime")), ThermTime(config.get<int>("Therm_Time")),
										ord_tab(as_vector<int>(config, "Order_Step_Table")), normmin(static_cast<int>(config.get<double>("Norm_Points")/config.get<double>("NCores"))), endmin(static_cast<int>(config.get<double>("End_Points")/config.get<double>("NCores"))), TotRunTime(config.get<int>("Total_Time")), TotMaxOrd(config.get<int>("Total_Max_Order")),
										which_fw_ad(config.get<int>("Which_fw_adapt()")), desi_rat(config.get<double>("Desi_fw_rat")),fw_max(config.get<double>("Fake_Weight_Max")),
										taumap(tmap(create_fvec(config), as_vector<int>(config, "Bins"), as_vector<double>(config, "Taus"), config.get<int>("Tminit"), config.get<int>("Tmaxit")))
{ 
  try{
	if (fabs(1-Prem-Pins-Pct -Pctho -Pic -Prc - Pdq-Psw) > 0.0000000001) {throw oor_Probs();}
#ifdef SECUMUL
	if (ThermTime > TotRunTime)
#else
	if (ThermTime > RunTime)
#endif
	{throw Therm();}

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
	Eptmp.assign(ws.size(), 0.);
	Epol.assign(ws.size(), 0.);
	SE = ArrayXXd::Zero(taumap.taubin,3);
	SEacc = create_empty_SE_acc("step0");
	SEacc << alps::accumulators::LogBinningAccumulator<std::vector<double>>("ord1");
	SEtmp.assign(taumap.taubin, 0.);
	G0SEiw = ArrayXXcd::Zero(wbin, 2);
	//binning
	counts = create_empty_count_acc(0);
	Epbin = create_empty_Ep_acc("step0"); //this is just possible if ws was initialzed before
	last_g0_count =0.;
	last_measured_Epol = Epol;
	Ep_intv = create_empty_Ep_acc("step0");
	last_count_g0attau0 = 0.;
	G0tau0 << alps::accumulators::FullBinningAccumulator<double>("tau0");
	G0tau0_each << alps::accumulators::LogBinningAccumulator<double>("tau0");
	ordesti = create_empty_ordesti(0);

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
	ordstsz = ord_tab[0];
	if (config.get<bool>("1st_Step_Max_Ord")){
		if ((maxord>TotMaxOrd) && (TotMaxOrd>0)) {maxord = TotMaxOrd;}
	} else {
		if ((TotMaxOrd < ordstsz) && (TotMaxOrd > 0)) {maxord = TotMaxOrd;}
	   	else {maxord = minord + ordstsz;}
	}
	SEib= ArrayXd::Zero(taumap.taubin);
	Norms = ArrayXd::Zero(taumap.taubin);
	Ends = ArrayXd::Zero(taumap.taubin);
	nnorms.assign(1, 0.);
	nends.assign(1, 0.);
	
	if (which_fw_ad==1) {
		fw_last.assign(2, 0);
		fw_counts.assign(2, 0.);
	} else if (which_fw_ad ==2){
		fw_last.assign(maxord-minord+1, 0);
		fw_counts.assign(maxord-minord+1, 0.);
	}
	fw_vec.assign(maxord-minord, 1.);
#endif
  
	updatestat = ArrayXXd::Zero(14,6);
	overflowstat = updatestat;
	updatestat.row(13) = ArrayXd::Ones(6);
	orderstat = ArrayXd::Zero(config.get<int>("Order_Stat_Size"));
	qstat = ArrayXXd::Zero(config.get<int>("QStat_Size"),2);
	tstat = ArrayXXd::Zero(taumap.taubin,4);
	testhisto = ArrayXXd::Zero(taumap.taubin+2, 4);
  
	path=config.get<std::string>("Path");
	
	global_weight= logG0el(diag.get_p(0), diag.get_tfin(0), diag.get_tinit(0));

	lu = "nothing";
	if (thread == 0){taumap.write_norm_tab();}
  } 
  catch (std::exception& e){
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}

DiagMC::~DiagMC() {

}

alps::accumulators::accumulator_set DiagMC::create_empty_Ep_acc(const std::string & name){
  accumulators_type Ep_empty;
  Ep_empty << alps::accumulators::LogBinningAccumulator<std::vector<double>>(name);
  return Ep_empty;
}

alps::accumulators::accumulator_set DiagMC::create_empty_SE_acc(const std::string & name){
  accumulators_type SE_empty;
  SE_empty << alps::accumulators::LogBinningAccumulator<std::vector<double>>(name);
  return SE_empty;
}

alps::accumulators::accumulator_set DiagMC::create_empty_count_acc(const int & order){
  accumulators_type count_empty;
  count_empty << alps::accumulators::LogBinningAccumulator<double>("norm"+std::to_string(order));
  count_empty << alps::accumulators::LogBinningAccumulator<double>("end"+std::to_string(order));
  return count_empty;
}

alps::accumulators::accumulator_set DiagMC::create_empty_ordesti(const int & order){
  accumulators_type ordesti_empty;
  ordesti_empty << alps::accumulators::LogBinningAccumulator<double>("step"+std::to_string(order));
  return ordesti_empty;
}

