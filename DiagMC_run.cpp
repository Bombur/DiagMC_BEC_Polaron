// Froehlich Polaron
// Hans Peter Guertner
#include <iostream>
#include <string>
#include <ctime>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "DiagMC.h"
#include "Adaption.h"
#include "RunException.h"
#include "tmap.h"
#include <assert.h>
#include <alps/accumulators.hpp>
#include <alps/hdf5.hpp>

#if defined(_OPENMP)
#include <omp.h>
#endif  

namespace pt = boost::property_tree;
using namespace std::chrono;
using namespace Eigen;


int main() {
	#if defined(_OPENMP)
	  const int nseeds = omp_get_max_threads();
	  //const int nseeds = 4;
	#else
	  const int nseeds = 1;
	#endif
	  
	std::cout << "# Using " << nseeds << " cores." << std::endl;
	
	pt::ptree config;
	pt::read_json("DiagMC_BEC.json", config);
	
	//constants and parameters
	double onecore_time= 0;
  	
	//container for pointer to all parallel DiagMCs to reuse parallelism
	std::vector<DiagMC *> DiagMC_pl(nseeds); 
	
	//Initialization
	#pragma omp parallel for num_threads(nseeds) 
	for (int seed = 0; seed < nseeds; seed++)
	{
		#pragma omp critical  // to have the same staring point
		{
		  DiagMC_pl[seed] = new DiagMC(seed,config);
		  std::cout << "# Initialized thread " << seed << std::endl;
		}
	}
	
	//Total Statistics and Data
	std::vector<int> bintemp = as_vector<int>(config, "Bins");
	std::vector<double> tautemp = as_vector<double>(config, "Taus");
	
	std::vector<ArrayXXd> Data_tot(nseeds);
	
	//Estimators
	std::vector<double> ws = as_vector<double>(config, "Ws_for_Epol");
	std::vector<ArrayXd> Epcum;										//save averages from Ep cumulative of each order step
	ArrayXXd Eptot= ArrayXXd::Zero(ws.size(), nseeds);				//Ep tot of each seed 
	//ArrayXXd Eptest = ArrayXXd::Zero(bintemp[bintemp.size()-1], ws.size());
	std::vector<ArrayXXd> SE(nseeds, ArrayXXd::Zero(bintemp[bintemp.size()-1], 4)); 
	std::vector<ArrayXXcd> G0SEiw(nseeds, ArrayXXcd::Zero(config.get<int>("Omega_bin"), 3)); 
	
	//count0 binning
	alps::accumulators::accumulator_set counts;
	counts = DiagMC_pl[0]->create_empty_count_acc(0);
	
	//Ep Binning
	alps::accumulators::accumulator_set Epbin;
	Epbin = DiagMC_pl[0]->create_empty_Ep_acc("step0");
	
	//Ep integral binning
	alps::accumulators::accumulator_set Ep_intv;
	Ep_intv = DiagMC_pl[0]->create_empty_Ep_acc("step0");
	
	//G0 Test Binning container
	//alps::accumulators::accumulator_set G0tau0;
	//G0tau0 = DiagMC_pl[0]->get_G0tau0();
	//alps::accumulators::result_set G0tau0_each;
	
	//Statistics
	ArrayXXd uds_tot = DiagMC_pl[0]->get_uds();
	ArrayXXd ofs_tot = DiagMC_pl[0]->get_ofs();
	ArrayXd os_tot = DiagMC_pl[0]->get_os();
	ArrayXd ts_tot = ArrayXd::Zero(10);
	ArrayXXd qs_tot = DiagMC_pl[0]->get_qs();
	ArrayXXd taus_tot = ArrayXXd::Zero(bintemp[bintemp.size()-1], 5); // Tau Histogram
	ArrayXXd testhisto = ArrayXXd::Zero(bintemp[bintemp.size()-1]+2, 4);
	
	//Total time Control
	steady_clock::time_point Cumt_be = steady_clock::now();
	steady_clock::time_point Cumdt_be = steady_clock::now();
#ifdef SECUMUL
	  double Cumt =0.;
	  double Cumdt= 0.;
	  
	  //Statistics and Data for SECUMUL
	  ArrayXXd SEib = ArrayXXd::Zero(bintemp[bintemp.size()-1],nseeds);
	  std::vector<ArrayXXd> Norms;
	  std::vector<ArrayXXd> Ends;
	  std::vector<double> nnorms;
	  std::vector<double> nends;
	  std::vector<double> mins;
	  std::vector<double> maxs;
	  std::vector<double> stptime;
	  std::vector<double> bottlenecks;
	  ArrayXXd fw_stat = ArrayXXd::Zero(1, nseeds); 

	  //temporary Containers to transfer Data
	  ArrayXXd Normstemp(bintemp[bintemp.size()-1],nseeds);
	  ArrayXXd Endstemp(bintemp[bintemp.size()-1],nseeds);
	  
	  // -------------------- Order Step Loop----------------
	  //order steps iterator
	  int ordit =0;
	  
	  double speedcumul = 0.;
	  
	  //Total Max Order
	  int TotMaxOrdcheck = -1;
	  
	  std::cout << "# Order Loop Starting! "<< std::endl;
	  do{
		std::cout <<'\n'<< "# ------ Order Step " << ordit << "------" <<'\n' << std::endl;
		Cumdt_be = steady_clock::now();
		
		//temporary Containers to transfer Data
		double nnormstemp = 0.;
	  	double nendstemp = 0.;
		double nnormstemp2 = 0.;
	  	double nendstemp2 = 0.;
		VectorXd fw_stat_tmp= VectorXd::Zero(nseeds);
		
		//variable to make all threads to do the next order step at the same time
		int norm_allthreads=0;
		int end_allthreads=0;
#endif
		
		
	  #pragma omp parallel for ordered schedule(static, 1)
	  for (int seed = 0; seed < nseeds; seed++)
	  {	
		DiagMC * bec = DiagMC_pl[seed];
  	  
		//time control
		steady_clock::time_point time_begin = steady_clock::now();  //start time
		ArrayXd timestat = ArrayXd::Zero(10);
		steady_clock::time_point dt_be = steady_clock::now();
  	
		double nseconds =0.;
		double dt= 0.;
		double speed = 0.; 
		
		//measure control
		int meas =0;
		
		//iterator
		int it= 1;
		
#ifdef SECUMUL//Control of Do Loop
		bool maxordcheck = false;
		bool normcheck = false;
		bool endcheck = false;
		bool normadd1 = false;
		bool endadd1 = false;
		double Tott_check=0.;
#endif
		
		do {
		  double action = bec->drnd();
			
		  if (action < bec->Prem) {
			steady_clock::time_point t_rem_be = steady_clock::now();
			bec->remove();
			steady_clock::time_point t_rem_end = steady_clock::now();
			timestat(3) += static_cast<double>((t_rem_end-t_rem_be).count())* steady_clock::period::num / steady_clock::period::den ;
							
		  }
		  else if ((action-bec->Prem)<bec->Pins) {
			steady_clock::time_point t_ins_be= steady_clock::now();
			bec->insert();
			steady_clock::time_point t_ins_end = steady_clock::now();
			timestat(2) += static_cast<double>((t_ins_end-t_ins_be).count())* steady_clock::period::num / steady_clock::period::den ;
		  }
		  else if ((action-bec->Prem-bec->Pins) < bec->Pct) {
			
			steady_clock::time_point t_ct_be= steady_clock::now();
			bec->change_tau();
			steady_clock::time_point t_ct_end = steady_clock::now();
			timestat(0) += static_cast<double>((t_ct_end-t_ct_be).count())* steady_clock::period::num / steady_clock::period::den ;
		  }
			
		  else if ((action-bec->Prem-bec->Pins-bec->Pct) < bec->Pdq) {
			steady_clock::time_point t_dq_be= steady_clock::now();
			bec->dq();
			steady_clock::time_point t_dq_end = steady_clock::now();
			timestat(5) += static_cast<double>((t_dq_end-t_dq_be).count())* steady_clock::period::num / steady_clock::period::den ;
				
		  }
			
		  else if ((action- bec->Prem- bec->Pins - bec->Pct - bec->Pdq) < bec->Psw) {
				  
			steady_clock::time_point t_sw_be= steady_clock::now();
			bec->swap();
			steady_clock::time_point t_sw_end = steady_clock::now();
			timestat(4) += static_cast<double>((t_sw_end-t_sw_be).count())* steady_clock::period::num / steady_clock::period::den ;
				 
		  }
				
		  else if ((action- bec->Prem- bec->Pins - bec->Pct - bec->Pdq- bec->Psw) < bec->Pctho) {
				  
			steady_clock::time_point t_ctho_be= steady_clock::now();
			bec->ct();
			steady_clock::time_point t_ctho_end = steady_clock::now();
			timestat(1) += static_cast<double>((t_ctho_end-t_ctho_be).count())* steady_clock::period::num / steady_clock::period::den ;
				 
		  }

		  else if ((action- bec->Prem- bec->Pins - bec->Pct - bec->Pdq - bec->Psw - bec->Pctho) < bec->Pic) {
				  
			steady_clock::time_point t_ic_be= steady_clock::now();
			bec->inscrossed();
			steady_clock::time_point t_ic_end = steady_clock::now();
			timestat(6) += static_cast<double>((t_ic_end-t_ic_be).count())* steady_clock::period::num / steady_clock::period::den ;
				 
		  }
				
		  else if ((action- bec->Prem- bec->Pins - bec->Pct - bec->Pdq - bec->Psw - bec->Pctho - bec->Pic) < bec->Prc) {
				  
			steady_clock::time_point t_rc_be = steady_clock::now();
			bec->remcrossed();
			steady_clock::time_point t_rc_end = steady_clock::now();
			timestat(7) += static_cast<double>((t_rc_end-t_rc_be).count())* steady_clock::period::num / steady_clock::period::den ;
		  }

		  else { 
			assert(0);
		  }
			
		  if ((it%bec->Meas_its) ==0) {
			steady_clock::time_point t_meas_be= steady_clock::now();
			double thtime =((t_meas_be-Cumdt_be).count()* steady_clock::period::num / steady_clock::period::den);
			if (thtime > bec->ThermTime){
				//bec->meas_histo();
				meas += bec->measure();
				bec->meas_ordstat();
				if ((meas % bec->Ep_meas_it)==0){
					//bec->meas_G0tau0();
					bec->meas_Ep_intv();
				}
			} else if (thtime<(bec->ThermTime*0.8)) {
				bec->meas_ordstat();
			} else {
				bec->set_os_to_zero();
			}	
			steady_clock::time_point t_meas_end = steady_clock::now();
			timestat(8) += static_cast<double>((t_meas_end-t_meas_be).count())* steady_clock::period::num / steady_clock::period::den;
		  }
			
		  if ( (it%bec->Test_its) ==0) {
			steady_clock::time_point t_test_be= steady_clock::now();
			bec->meas_histo(); 
			try{
			  bec->test();
			} catch (std::exception& e){
			  std::cerr << e.what() << std::endl;
			  exit(EXIT_FAILURE);
			}
			std::cout << "# Test ok thread " << seed << std::endl;
			steady_clock::time_point t_test_end = steady_clock::now();
			//exit(EXIT_FAILURE);
			timestat(9) += static_cast<double>((t_test_end-t_test_be).count())* steady_clock::period::num / steady_clock::period::den ;
		  }
		  
#ifndef SECUMUL
		  if ( (it%bec->Write_its) ==0) {
			//time 
			steady_clock::time_point time_end = steady_clock::now();
			nseconds = static_cast<double>(duration_cast<seconds>(time_end-time_begin).count());
		  	steady_clock::time_point dt_end = steady_clock::now();
			dt = static_cast<double>(duration_cast<seconds>(dt_end-dt_be).count());
			
			speed = static_cast<double>(it)/dt;
			std::cout << "# " << seed<<  ": Time on core  " << nseconds << '\t' << "dt  " << dt << '\t' << "speed [its/s]: " << speed << '\n' << std::endl;
			it = 0;
			dt_be = steady_clock::now();
		  }	
		  it++;

		} while (nseconds < bec->RunTime - dt);
		
		#pragma omp critical  // to have the same staring point
		{	
		  std::cout << "# Combine thread " << seed << std::endl;
		  Data_tot[seed] = bec->get_Data();
		  Eptot.col(seed) = bec->get_Ep();
		  //Eptest += bec->get_Eptest();
		  SE[seed] = bec->get_SE();
		  G0SEiw[seed] = bec->get_G0SEiw();
		  
		  //binning
		  counts.merge(bec->get_counts());
		  if (bec->Ep_bin_each_step){Epbin.merge(bec->get_Epacc());}
		  Ep_intv.merge(bec->get_Ep_intv());
		  //G0tau0 = bec->get_G0tau0();
		  //G0tau0_each = bec->get_G0tau0_each();
		  
		  uds_tot += bec->get_uds();
		  ofs_tot += bec->get_ofs();
		  os_tot += bec->get_os();
		  ts_tot += timestat;
		  qs_tot += bec->get_qs();
		  taus_tot += bec->get_taus();
		  testhisto += bec->get_testhisto();
		  bec->get_ordesti(seed);
	
		  onecore_time += nseconds;
		}
	  }//end of omp for loop
	  
	  std::cout << "Writing binning results" << std::endl;
	  //count0 binning analysis
	  alps::accumulators::result_set counts_res(counts);
	  counts_res["norm0"] *= counts_res["norm0"].count();
	  
	  alps::accumulators::result_set Ep_intv_res(Ep_intv);

	  //alps::accumulators::result_set G0tau0_res(G0tau0);
	  //std::cout << G0tau0_res["tau0"] << std::endl;
	  //std::cout << G0tau0_each["tau0"] << std::endl;
	  
	  alps::accumulators::result_set Epbin_res(Epbin);
	  if (DiagMC_pl[0]->Ep_bin_each_step){
		 Epbin_res["step0"] *= Epbin_res["step0"].count()*DiagMC_pl[0]->get_Ep_norm()/counts_res["norm0"].mean<double>();
		}
	  }
	  
	  // Note the file is opened with write permission.
	  std::string Ep_filename("data/Ep/Epbin.h5");
	  alps::hdf5::archive oar(Ep_filename, "a");
	  oar["Ep"] << Epbin_res;
	  oar["counts"]<<counts_res;
	  oar["Ep_intv"] << Ep_intv_res;
	  //oar["G0tau0"]<< G0tau0_res;
	  //oar["G0tau0_each"]<< G0tau0_each;
	  oar.close();
	
	  
#else	//-------------2nd Part of Order Loop--------------------------

		  if ( (it%bec->Write_its) ==0) {
			//fw_adapt()
			double thtime =static_cast<double>((steady_clock::now()-Cumdt_be).count())* steady_clock::period::num / steady_clock::period::den;
			if (thtime < (0.8*bec->ThermTime)){ 
				if (bec->get_which_fw_adapt()==1) {bec->fw_adapt();}
				else if (bec->get_which_fw_adapt()==2) {bec->fw_vec_adapt();}
			} else if (thtime > bec->ThermTime){
				//Check if Norm and End numbers are reached
				if ((bec->get_normostat() > bec->normmin-1) && !normadd1) {
					normadd1 = true;
					#pragma omp atomic update
						norm_allthreads++;
				}
				if ((bec->get_endostat() > bec->endmin-1) && !endadd1) {
					endadd1 = true;
					#pragma omp atomic update
						end_allthreads++;
				}
			
				int ntmp;
				int etmp;
				#pragma omp atomic read
					ntmp = norm_allthreads;
				#pragma omp atomic read
					etmp = end_allthreads;
			  
				if (ntmp == nseeds) {normcheck = true;} // checking if we reached the norm minimum
				if (etmp == nseeds) {endcheck = true;} // checking if we reached the end minimum
				std::cout <<'\n' << bec->get_normostat() << '\t' << bec->normmin-1 << '\t' << bec->get_endostat() << '\t'<< bec->endmin-1 << std::endl;
				std::cout << normadd1 << '\t' << ntmp << '\t'<< endadd1 << '\t'<< etmp << '\n' << std::endl;
			}

			//time 
			steady_clock::time_point time_end = steady_clock::now();
			nseconds = static_cast<double>(duration_cast<seconds>(time_end-time_begin).count());
		  	steady_clock::time_point dt_end = steady_clock::now();
			dt = static_cast<double>(duration_cast<seconds>(dt_end-dt_be).count());
		  
			std::cout << "# " << seed<<  ": Time in Order Step  " << nseconds << '\t' << "dt  " << dt << '\n' << std::endl;
			dt_be = steady_clock::now();
			
			//To check Total Run Time
			steady_clock::time_point Tott_check_end = steady_clock::now();
			Tott_check = static_cast<double>(duration_cast<seconds>(Tott_check_end-Cumt_be).count());
			
			speedcumul+=static_cast<double>(it);
			it=0;
		  }	
		  it++;
			
		  //to go from one order step to the next 
		  //the diagram has to be in the maximum order of before
		  if (nseconds > bec->RunTime - dt && normcheck && endcheck) {
			if(bec->get_order()== bec->get_max_order()) {maxordcheck = true;}
		  }
		  
		} while ((nseconds < bec->RunTime - dt || !(normcheck) || !(endcheck) ||!(maxordcheck)) && (Tott_check < bec->TotRunTime - dt));
		
		
		#pragma omp critical  // to have the same staring point
		{
		  //first we need to calculate the normalization
		  nnormstemp += bec->normcalc(); 
		  nendstemp += bec->endcalc();
		  nnormstemp2 += bec->get_normostat(); 
		  nendstemp2 += bec->get_endostat();
		  
		  Eptot.col(seed) += bec->get_Ep();
		  if (ordit == 0) {Data_tot[seed] = bec->get_Data();}
		  ArrayXXd tmp = bec->get_SE();
		  SE[seed].rightCols(tmp.cols()-1) += tmp.rightCols(tmp.cols()-1);
		  SE[seed].col(0) = tmp.col(0);
		  ArrayXXcd tmp2 = bec->get_G0SEiw(); 
		  G0SEiw[seed].rightCols(tmp2.cols()-1) += tmp2.rightCols(tmp2.cols()-1);
		  G0SEiw[seed].col(0) = tmp2.col(0); 
		  
		  //binning
		  if (bec->Ep_bin_each_step){Epbin.merge(bec->get_Epacc());}
		  Ep_intv.merge(bec->get_Ep_intv());
		  counts.merge(bec->get_counts());
		  
		  //stats
		  uds_tot += bec->get_uds();
		  ofs_tot += bec->get_ofs();
		  os_tot += bec->get_os();
		  ts_tot += timestat;
		  qs_tot += bec->get_qs();
		  tmp = bec->get_taus();
		  taus_tot.rightCols(tmp.cols()-1) += tmp.rightCols(tmp.cols()-1);
		  taus_tot.col(0) = tmp.col(0);
		  fw_stat_tmp(seed) = bec->get_fw();
		  bec->get_ordesti(seed);
	
		  onecore_time += nseconds;
		  
		  //transfering Data
		  SEib.col(seed) += bec->get_SEib();
		  Normstemp.col(seed) = bec->get_NormDiag();
		  Endstemp.col(seed) = bec->get_EndDiag();
		  std::cout << "# Transfered Data thread " << seed << std::endl;
		}
		} //end of omp for loop
		
		//Not Parallel anymore
		//Combining Data
		{
		  //Estimator calc
		  Epcum.push_back(Eptot.rowwise().mean());
	  	  		  
	  	  //To be able to always print a Norm End comparison 
		  if (ordit == 0) {
			Norms.push_back(Normstemp);
			Norms.push_back(Normstemp);
			Ends.push_back(Endstemp);
			Ends.push_back(Endstemp);
		  } else if (ordit == 1) {
			Norms[1] = Normstemp;
			Ends[1] = Endstemp;
		  } else if (ordit == 2 || ordit==3){
			Norms.push_back(Normstemp);
			Ends.push_back(Endstemp);
		  } else {
			Norms[2]=Norms[3];
			Norms[3]=Normstemp;
			Ends[2]=Ends[3];
			Ends[3]=Endstemp;
		  }	
		  nnormstemp /= static_cast<double>(nseeds);
		  nendstemp /= static_cast<double>(nseeds);
		  nnorms.push_back(nnormstemp2);
		  nends.push_back(nendstemp2);
		  mins.push_back((DiagMC_pl[0]->get_minmax())[0]);
		  maxs.push_back((DiagMC_pl[0]->get_minmax())[1]);
		  eigen_push_back(fw_stat, fw_stat_tmp);
		  std::cout << "# Combined Data!" << std::endl;
		}
			  
		//count binning analysis
		alps::accumulators::result_set counts_res(counts);
		counts_res["norm"+std::to_string(ordit)] *= counts_res["norm"+std::to_string(ordit)].count();
		counts_res["end"+std::to_string(ordit)] *= counts_res["end"+std::to_string(ordit)].count();
		
		std::cout << "# Counts comparison: " << nnormstemp*static_cast<double>(nseeds) << '\t' << counts_res["norm"+std::to_string(ordit)].mean<double>()<<std::endl;
		
		alps::accumulators::result_set Ep_intv_res(Ep_intv);

		alps::accumulators::result_set Epbin_res(Epbin);
		if (DiagMC_pl[0]->Ep_bin_each_step){
			Epbin_res["step"+std::to_string(ordit)] *= Epbin_res["step"+std::to_string(ordit)].count()*DiagMC_pl[0]->get_Ep_norm()/counts_res["norm"+std::to_string(ordit)].mean<double>();
		}
	
		// Note the file is opened with write permission.
		std::cout << "# Writing binning results"<< std::endl;
		std::string filename("data/Ep/Epbin.h5");
		alps::hdf5::archive oar(filename, "a");
		oar["Ep"] << Epbin_res;
		oar["counts"] << counts_res;
		oar["Ep_intv"] << Ep_intv_res;
		oar.close();
		
		//setting binning container to zero again
		Epbin = DiagMC_pl[0]->create_empty_Ep_acc("step"+std::to_string(ordit+1));
		Ep_intv = DiagMC_pl[0]->create_empty_Ep_acc("step"+std::to_string(ordit+1));
		counts = DiagMC_pl[0]->create_empty_count_acc(ordit+1);
		
		// time per  order step
		steady_clock::time_point Cumdt_end = steady_clock::now();
		Cumdt = static_cast<double>(duration_cast<seconds>(Cumdt_end-Cumdt_be).count());
		speedcumul /= static_cast<double>(nseeds);
		speedcumul /= Cumdt;
		
		//Save Step Time for Order Step
		stptime.push_back(Cumdt);
		
		//Order Step Adaption
		Adaption OSA(DiagMC_pl[0]->get_ordstsz(), Cumdt, nnormstemp, nendstemp, config);
		bottlenecks.push_back(static_cast<double>(OSA.whichbtlnk()));
		if (config.get<bool>("ADAPTION")){
		  int tempordstsz = OSA.ordstszadapt();
		  std::cout<< "#Adapted Order Step Size from " << DiagMC_pl[0]->get_ordstsz() << " to " << tempordstsz << '\n';
		  for (auto i : DiagMC_pl) {
			i->set_ordstsz(tempordstsz);
		  }
		}
				
	  //Check if we reached Total Maximum oder	
	  	if ((DiagMC_pl[0]->TotMaxOrd == -1) || (DiagMC_pl[0]->TotMaxOrd > DiagMC_pl[0]->get_order() + DiagMC_pl[0]->get_ordstsz())) {
		  #pragma omp parallel for ordered schedule(static, 1)
		  for (int seed = 0; seed < nseeds; seed++)
		  {
			if (config.get<bool>("Averge_Norm_Stpw")){
			  DiagMC_pl[seed]->set_av_nei(nnormstemp, nendstemp, ordit);
			}
			DiagMC_pl[seed]->ord_step();
		  }
		} else if (DiagMC_pl[0]->TotMaxOrd == DiagMC_pl[0]->get_order()) { // we reached TotMaxOrd
		  TotMaxOrdcheck = ordit + 1;
		} else { // we have to change the last order step to reach Tot Max Ord
		  const int laststsz = DiagMC_pl[0]->TotMaxOrd - DiagMC_pl[0]->get_order();
		  std::cout<<"#Last Step! Step Size changed to " << laststsz << std::endl;
		  #pragma omp parallel for ordered schedule(static, 1)
		  for (int seed = 0; seed < nseeds; seed++)
		  {	
			if (config.get<bool>("Averge_Norm_Stpw")){
			  DiagMC_pl[seed]->set_av_nei(nnormstemp, nendstemp, ordit);
			}
			DiagMC_pl[seed]->set_ordstsz(laststsz);
			DiagMC_pl[seed]->ord_step();
		  }
		  TotMaxOrdcheck = ordit + 2;
		}
		  
		
		ordit++;
		
		steady_clock::time_point Cumt_end = steady_clock::now();
		Cumt = static_cast<double>(duration_cast<seconds>(Cumt_end-Cumt_be).count());
		
		std::cout << "\n#Total Time " << Cumt << '\t' << "dt  " << Cumdt << '\t' << "speed [its/s]: " << speedcumul <<'\n' << std::endl;
		
		} while ((Cumt < DiagMC_pl[0]->TotRunTime - Cumdt) && (TotMaxOrdcheck != ordit) );
#endif
	  
	  //Deleting the vector of DiagMC Pointers
	  #pragma omp for ordered schedule(static, 1)
	  for (int seed = 0; seed < nseeds; seed++)
	  {	
		delete DiagMC_pl[seed];
	  }
	
	//Eptest /= static_cast<double>(nseeds);
	uds_tot /= static_cast<double>(nseeds);
	ofs_tot /= static_cast<double>(nseeds);
	uds_tot.col(4)=uds_tot.col(3)/uds_tot.col(1);
	uds_tot.col(5)=uds_tot.col(3)/uds_tot.col(0);
	ofs_tot.col(4)=ofs_tot.col(2)/uds_tot.col(2);
	ofs_tot.col(5)=ofs_tot.col(3)/uds_tot.col(3);
	os_tot/=  static_cast<double>(nseeds);
	ts_tot /= onecore_time;
	qs_tot /= static_cast<double>(nseeds);
	taus_tot /= static_cast<double>(nseeds);
	testhisto /= static_cast<double>(nseeds);
	
	//averages	
	//all
#ifdef SECUMUL
	ArrayXXd all_orders(Data_tot[0].rows(), 3);  // 0:tau, 1:average, 2:error
	all_orders << Data_tot[0].col(0), ArrayXXd::Zero(Data_tot[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  all_orders.col(1) += SEib.col(i);
	all_orders.col(1) /= static_cast<double>(nseeds);
	
	//Compare first end and Norm Diagram
	ArrayXXd Norm_first(Data_tot[0].rows(), 3);  // 0:tau, 1:average, 2:error
	Norm_first << Data_tot[0].col(0), ArrayXXd::Zero(Data_tot[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  Norm_first.col(1) += Norms[1].col(i);
	Norm_first.col(1) /= static_cast<double>(nseeds);
	
	ArrayXXd End_first(Data_tot[0].rows(), 3);  // 0:tau, 1:average, 2:error
	End_first << Data_tot[0].col(0), ArrayXXd::Zero(Data_tot[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  End_first.col(1) += Ends[0].col(i);
	End_first.col(1) /= static_cast<double>(nseeds);
	
	//Compare last end and Norm Diagram
	ArrayXXd Norm_last(Data_tot[0].rows(), 3);  // 0:tau, 1:average, 2:error
	Norm_last << Data_tot[0].col(0), ArrayXXd::Zero(Data_tot[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  Norm_last.col(1) += Norms[Norms.size()-1].col(i);
	Norm_last.col(1) /= static_cast<double>(nseeds);
	
	ArrayXXd End_last(Data_tot[0].rows(), 3);  // 0:tau, 1:average, 2:error
	End_last << Data_tot[0].col(0), ArrayXXd::Zero(Data_tot[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  End_last.col(1) += Ends[Ends.size()-2].col(i);
	End_last.col(1) /= static_cast<double>(nseeds);
#else
	ArrayXXd all_orders(Data_tot[0].rows(), 3);  // 0:tau, 1:average, 2:error
	all_orders << Data_tot[0].col(0), ArrayXXd::Zero(Data_tot[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  all_orders.col(1) += Data_tot[i].col(1);
	all_orders.col(1) /= static_cast<double>(nseeds);
#endif
	
	//zero order
	ArrayXXd zero_order(Data_tot[0].rows(), 3);  // 0:tau, 1:average, 2:error
	zero_order << Data_tot[0].col(0) , ArrayXXd::Zero(Data_tot[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  zero_order.col(1) += Data_tot[i].col(2);
	zero_order.col(1) /=  static_cast<double>(nseeds);
	
	//first order
	ArrayXXd first_order(Data_tot[0].rows(), 3);  // 0:tau, 1:average, 2:error
	first_order << Data_tot[0].col(0) , ArrayXXd::Zero(Data_tot[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  first_order.col(1) += Data_tot[i].col(3);
	first_order.col(1) /=  static_cast<double>(nseeds);
	
	//second order
	ArrayXXd second_ordera(Data_tot[0].rows(), 3);  // 0:tau, 1:average, 2:error
	second_ordera << Data_tot[0].col(0) , ArrayXXd::Zero(Data_tot[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  second_ordera.col(1) += Data_tot[i].col(4);
	second_ordera.col(1) /=  static_cast<double>(nseeds);
	
	ArrayXXd second_orderb(Data_tot[0].rows(), 3);  // 0:tau, 1:average, 2:error
	second_orderb << Data_tot[0].col(0) , ArrayXXd::Zero(Data_tot[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  second_orderb.col(1) += Data_tot[i].col(5);
	second_orderb.col(1) /=  static_cast<double>(nseeds);
	
	//errors
	if (nseeds > 3) {
	  for (int i=0; i<nseeds ; i++) {
#ifdef SECUMUL
		all_orders.col(2) += (SEib.col(i)-all_orders.col(1)).pow(2); 
		Norm_first.col(2) += (Norms[1].col(i)-Norm_first.col(1)).pow(2);
		End_first.col(2) += (Ends[0].col(i)-End_first.col(1)).pow(2);
		Norm_last.col(2) += (Norms[Norms.size() -1].col(i)-Norm_last.col(1)).pow(2);
		End_last.col(2) += (Ends[Ends.size() -2].col(i)-End_last.col(1)).pow(2);
#else
		all_orders.col(2) += (Data_tot[i].col(1)-all_orders.col(1)).pow(2); 
#endif
		zero_order.col(2) += (Data_tot[i].col(2)-zero_order.col(1)).pow(2);
		first_order.col(2) += (Data_tot[i].col(3)-first_order.col(1)).pow(2); 
		second_ordera.col(2) += (Data_tot[i].col(4)-second_ordera.col(1)).pow(2); 
		second_orderb.col(2) += (Data_tot[i].col(5)-second_orderb.col(1)).pow(2); 
	  }
	  
	  double norm = 1/static_cast<double>(nseeds*(nseeds-1));
#ifdef SECUMUL
	  Norm_first.col(2) *= norm;
	  Norm_first.col(2) = Norm_first.col(2).sqrt();
	  End_first.col(2) *= norm;
	  End_first.col(2) = End_first.col(2).sqrt();
	  Norm_last.col(2) *= norm;
	  Norm_last.col(2) = Norm_last.col(2).sqrt();
	  End_last.col(2) *= norm;
	  End_last.col(2) = End_last.col(2).sqrt();
#endif
	  all_orders.col(2) *= norm;
	  all_orders.col(2) = all_orders.col(2).sqrt();
	  zero_order.col(2) *= norm;
	  zero_order.col(2) = zero_order.col(2).sqrt();
	  first_order.col(2) *= norm;
	  first_order.col(2) = first_order.col(2).sqrt();
	  second_ordera.col(2) *= norm;
	  second_ordera.col(2) = second_ordera.col(2).sqrt();
	  second_orderb.col(2) *= norm;
	  second_orderb.col(2) = second_orderb.col(2).sqrt();
	  }
	  
	std::ofstream all("data/all_orders");  
	all << all_orders << '\n';
	all.close();
	
	std::ofstream zero("data/zero_order");
	zero << zero_order << '\n';
	zero.close();
	
	std::ofstream first("data/first_order");
	first << first_order << '\n';
	first.close();
	
	std::ofstream seconda("data/second_ordera");
	seconda << second_ordera << '\n';
	seconda.close();
	
	std::ofstream secondb("data/second_orderb");
	secondb << second_orderb << '\n';
	secondb.close();

#ifdef SECUMUL
	std::ofstream Nf("data/secumul/Norm_first");  
	Nf << Norm_first << '\n';
	Nf.close();
	
	std::ofstream Ef("data/secumul/End_first");  
	Ef << End_first << '\n';
	Ef.close();
	
	std::ofstream Nl("data/secumul/Norm_last");  
	Nl << Norm_last << '\n';
	Nl.close();
	
	std::ofstream El("data/secumul/End_last");  
	El << End_last << '\n';
	El.close();
#endif

	
//---------------------------------------------------ESTIMATORS
	ArrayXXd Eptot_out(ws.size(),nseeds + 1); 
	Eptot_out.col(0) =  Map<const Array<double, 1, Dynamic> > (ws.data(), ws.size());   //Max Order of Step
	Eptot_out.rightCols(nseeds) = Eptot;
	
	std::ofstream Eptot_outfile("data/Ep/Epvsws");
	Eptot_outfile << Eptot_out << '\n';
	Eptot_outfile.close();
	
	/*
	ArrayXXd Eptest_out(Eptest.rows(), Eptest.cols()+1); 
	Eptest_out.col(0) = Data_tot[0].col(0);   //tau
	Eptest_out.rightCols(Eptest.cols()) = Eptest;
	
	std::ofstream Eptest_outfile("data/Ep/Eptest");
	Eptest_outfile << Eptest_out << '\n';
	Eptest_outfile.close();
	*/

#ifdef SELFENERGY
#ifdef SECUMUL
	ArrayXXd Ep_out(Epcum.size(), ws.size()+1);
	Ep_out.col(0) =  Map<const Array<double, 1, Dynamic> > (maxs.data(), mins.size());   //Max Order of Step
	for (int i =0; i < Epcum.size() ; i++)
	  Ep_out.rightCols(ws.size()).row(i) = Epcum[i];
	
	std::ofstream Ep_outfile("data/Ep/Epvsord");
	Ep_outfile << Ep_out << '\n';
	Ep_outfile.close();
#endif
#endif
	
	//SE 1st
	ArrayXXd SE_1st(SE[0].rows(), 3);  // 0:tau, 1:average, 2:error
	SE_1st << SE[0].col(0) , ArrayXXd::Zero(SE[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  SE_1st.col(1) += SE[i].col(2);
	SE_1st.col(1) /=  static_cast<double>(nseeds);
	
	//SE >1
	ArrayXXd SE_big1(SE[0].rows(), 3);  // 0:tau, 1:average, 2:error
	SE_big1 << SE[0].col(0) , ArrayXXd::Zero(SE[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  SE_big1.col(1) += SE[i].col(3);
	SE_big1.col(1) /=  static_cast<double>(nseeds);
	
	//SE all
	ArrayXXd SE_all(SE[0].rows(), 3);  // 0:tau, 1:average, 2:error
	SE_all << SE[0].col(0) , ArrayXXd::Zero(SE[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  SE_all.col(1) += SE[i].col(1);
	SE_all.col(1) /=  static_cast<double>(nseeds);
	
	if (nseeds >3) {
	  for (int i=0; i<nseeds; i++){
		SE_1st.col(2) += (SE[i].col(2)-SE_1st.col(1)).pow(2);
		SE_big1.col(2) += (SE[i].col(3)-SE_big1.col(1)).pow(2);
		SE_all.col(2) += (SE[i].col(1)-SE_all.col(1)).pow(2);
	  }
	  double norm = 1/static_cast<double>(nseeds*(nseeds-1));
	  SE_1st.col(2) *= norm;
	  SE_1st.col(2) = SE_1st.col(2).sqrt();
	  SE_big1.col(2) *= norm;
	  SE_big1.col(2) = SE_big1.col(2).sqrt();
	  SE_all.col(2) *= norm;
	  SE_all.col(2) = SE_all.col(2).sqrt();
	}
	
	std::ofstream SE_1st_out("data/SE/SE_1");  
	SE_1st_out<< SE_1st << '\n';
	SE_1st_out.close();

	std::ofstream SE_big1_out("data/SE/SE_>1");  
	SE_big1_out<< SE_big1 << '\n';
	SE_big1_out.close();
	
	std::ofstream SE_all_out("data/SE/SE_all");  
	SE_all_out<< SE_all << '\n';
	SE_all_out.close();
	
	//G0SEiw 1st
	ArrayXXcd G0SEiw_1st(G0SEiw[0].rows(), 3);  // 0:omega, 1:average, 2:error
	G0SEiw_1st << G0SEiw[0].col(0) , ArrayXXcd::Zero(G0SEiw[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  G0SEiw_1st.col(1) += G0SEiw[i].col(2);
	G0SEiw_1st.col(1) /=  static_cast<double>(nseeds);
	
	//SE all
	ArrayXXcd G0SEiw_all(G0SEiw[0].rows(), 3);  // 0:tau, 1:average, 2:error
	G0SEiw_all << G0SEiw[0].col(0) , ArrayXXcd::Zero(G0SEiw[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  G0SEiw_all.col(1) += G0SEiw[i].col(1);
	G0SEiw_all.col(1) /=  static_cast<double>(nseeds);

	if (nseeds >3) {
	  for (int i=0; i<nseeds; i++){
		G0SEiw_1st.col(2) += (G0SEiw[i].col(2)-G0SEiw_1st.col(1)).pow(2);
		G0SEiw_all.col(2) += (G0SEiw[i].col(1)-G0SEiw_all.col(1)).pow(2);
	  }
	  double norm = 1/static_cast<double>(nseeds*(nseeds-1));
	  G0SEiw_1st.col(2) *= norm;
	  G0SEiw_1st.col(2) = G0SEiw_1st.col(2).sqrt();
	  G0SEiw_all.col(2) *= norm;
	  G0SEiw_all.col(2) = G0SEiw_all.col(2).sqrt();
	}

	std::ofstream G0SEiw_1st_out("data/G0SEiw/G0SEiw_1");  
	G0SEiw_1st_out<< G0SEiw_1st << '\n';
	G0SEiw_1st_out.close();

	std::ofstream G0SEiw_all_out("data/G0SEiw/G0SEiw_all");  
	G0SEiw_all_out << G0SEiw_all << '\n';
	G0SEiw_all_out.close();
	  
	
// ------------------------------------------------------STATISTICS
	std::ofstream uds_outfile("data/stats/uds_tot"); 
	std::ofstream ofs_outfile("data/stats/ofs_tot"); 
	std::ofstream os_outfile("data/stats/os_tot");
	std::ofstream ts_outfile("data/stats/ts_tot");
	std::ofstream qs_outfile("data/stats/qs_tot");
	std::ofstream taus_outfile("data/stats/taus_tot");
	
	std::string udcols[] = {"ATTEMPTED", "POSSIBLE", "REJECTED", "ACCEPTED", "ACCEPTANCE RATIO POSSIBLE", "ACCEPTANCE RATIO TOTAL"};
	std::string udrows[] = {"CHANGE TAU", "INSERT", "REMOVE", "SWAP", "SWAPoc", "SWAPco", "SWAPoocc", "CT in HO", "DQ", "IC", "RC", "FOIns", "FORem"};

	
	uds_outfile << "Update Statistics" << '\n';
    uds_outfile << "*************************************" << '\n';
    uds_outfile << "COLUMNS" << '\n';
	int it =0;
	for (auto name : udcols){
	  uds_outfile << it<<":"<< name << '\n';
	  it++;
	}
    uds_outfile << "*************************************" << '\n';
	it =0;
	for (auto name : udrows){
	  uds_outfile << name << ": " <<'\t'<< uds_tot.block(it, 0, 1, 6) << '\n';
	  it++;
	}
    uds_outfile << '\n';
    
    std::string ofcols[] = {"OVERFLOW", "UNDERFLOW", "REJECTED BY FLOWERROR", "ACCEPTED BY FLOWERROR", "RELATIVE REJECTION" , "RELATIVE ACCEPTANCE"};
	ofs_outfile << "Overflow Statistics" << '\n';
    ofs_outfile << "*************************************" << '\n';
    ofs_outfile << "COLUMNS" << '\n';
	it =0;
	for (auto name : ofcols){
	  ofs_outfile << it<<":"<< name << '\n';
	  it++;
	}
    ofs_outfile << "*************************************" << '\n';
	it =0;
	for (auto name : udrows){
	  ofs_outfile << name << ": " <<'\t'<< ofs_tot.block(it, 0, 1, 6) << '\n';
	  it++;
	}
	ofs_outfile << "Ep binning: " <<'\t'<< ofs_tot.block(13, 0, 1, 6) << '\n';
    ofs_outfile << '\n';
	
	MatrixXd orderprint(os_tot.size(),2);
	orderprint << ArrayXd::LinSpaced(config.get<int>("Order_Stat_Size"), 0, config.get<double>("Order_Stat_Size")-1), os_tot;
	
	os_outfile << "Order Statistics" << '\n';
	os_outfile << orderprint  << '\n';
	
	std::string tsrows[] = {"CHANGE TAU", "CT in HO", "INSERT", "REMOVE", "SWAP", "DQ", "IC", "RC", "MEAS", "TEST"};
	ts_outfile << "Time Statistics [%]" << std:: endl;
	it =0;
	for (auto name :tsrows) {
	  ts_outfile << name << ": " << '\t' << ts_tot(it) << '\n';
	  it +=1;
	}
	
	ArrayXXd qprint(qs_tot.rows(),3);
	qprint.col(0) = ArrayXd::LinSpaced(config.get<int>("QStat_Size"), 0, config.get<double>("Q_Cutoff"));
	qprint.rightCols(2) << qs_tot;
	qs_outfile << "Phonon Momentum Statistics" << '\n';
	qs_outfile << qprint  << '\n';
	
	taus_outfile << "Tau Histogram" << '\n';
	taus_outfile << taus_tot  << '\n';
	
	uds_outfile.close();
	ofs_outfile.close();
	os_outfile.close();
	ts_outfile.close();
	qs_outfile.close();
	taus_outfile.close();
	
	std::ofstream testhisto_of("data/stats/testhisto");
	testhisto_of << testhisto << '\n';
	testhisto_of.close();
	
#ifdef SECUMUL
	ArrayXXd minmax_out(mins.size(), 7);
	minmax_out.col(0) = Map<const Array<double, 1, Dynamic> > (mins.data(), mins.size());
	minmax_out.col(1) = Map<const Array<double, 1, Dynamic> > (maxs.data(), maxs.size());
	minmax_out.col(2) = Map<const Array<double, 1, Dynamic> > (nnorms.data(), nnorms.size());
	minmax_out.col(3) = Map<const Array<double, 1, Dynamic> > (nends.data(), nends.size());
	minmax_out.col(4) = Map<const Array<double, 1, Dynamic> > (stptime.data(), stptime.size());
	minmax_out.col(5) = Map<const Array<double, 1, Dynamic> > (bottlenecks.data(), bottlenecks.size());
	minmax_out.col(6) = minmax_out.col(3) / minmax_out.col(2);
	std::ofstream minmax_outfile("data/stats/minmax");
	minmax_outfile << "Minimum Maximum Order Statistics" << '\n';
    minmax_outfile << "*************************************" << '\n';
	minmax_outfile << "Min \t Max  \t nMin \t nMax \t time [s] \t Bottleneck \t nRatio\n" ;
	minmax_outfile << minmax_out;
	minmax_outfile.close();
	
	std::ofstream fw_outfile("data/stats/fw_stat");
	fw_outfile << fw_stat;
	fw_outfile.close();
#endif	
	
//	std::cout << "# Analysing data, result = " << system(config.get<std::string>("Ana_File").c_str()) << std::endl;
  
	return 0;
}



  
