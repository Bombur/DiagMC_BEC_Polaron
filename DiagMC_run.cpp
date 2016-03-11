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
		DiagMC_pl[seed] = new DiagMC(seed,config);
		
		#pragma omp critical  // to have the same staring point
		{
		  std::cout << "# Initialized thread " << seed << std::endl;
		}
	}
	
	//Total Statistics and Data
	std::vector<int> bintemp = as_vector<int>(config, "Bins");
	
	std::vector<ArrayXXd> Data_tot(nseeds);
	
	std::vector<double> ws = as_vector<double>(config, "Ws_for_Epol");
	std::vector<ArrayXd> Ep_tot;										//save averages from Ep of each order step
	ArrayXXd Epcum = ArrayXXd::Zero(ws.size(), nseeds);					//Ep cumulative of each seed 
	//ArrayXXd Eptest = ArrayXXd::Zero(bintemp[bintemp.size()-1], ws.size());

	ArrayXXd uds_tot = ArrayXXd::Zero(11,6);
	ArrayXd os_tot = ArrayXd::Zero(config.get<int>("Order_Stat_Size"));
	ArrayXd ts_tot = ArrayXd::Zero(10);
	ArrayXd qs_tot = ArrayXd::Zero(config.get<int>("QStat_Size"));
	
	  
#ifdef SECUMUL
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

	  	//Total time Control
	  steady_clock::time_point Cumt_be = steady_clock::now();  //start time
	  steady_clock::time_point Cumdt_be;
	  double Cumt =0.;
	  double Cumdt= 0.;
	  
	  //temporary Containers to transfer Data
	  ArrayXXd SEibtemp(bintemp[bintemp.size()-1],nseeds);
	  ArrayXXd Normstemp(bintemp[bintemp.size()-1],nseeds);
	  ArrayXXd Endstemp(bintemp[bintemp.size()-1],nseeds);
	  
	  // -------------------- Order Step Loop----------------
	  //order steps iterator
	  int ordit =0;
	  
	  //Total Max Order
	  int TotMaxOrdcheck = -1;
	  
	  std::cout << "# Order Loop Starting! "<< std::endl;
	  do{
		std::cout <<'\n'<< "# ------ Order Step " << ordit << "------" <<'\n' << std::endl;
		steady_clock::time_point Cumdt_be = steady_clock::now();
		
		 //temporary Containers to transfer Data
		double nnormstemp = 0.;
	  	double nendstemp = 0.;
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
		
		//measure control
		int meas =0;
		
		//iterator
		int it= 1;
		
#ifdef SECUMUL//Control of Do Loop
		bool maxordcheck = false;
		bool normcheck = false;
		bool endcheck = false;
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

#ifdef FP
		  else if ((action- bec->Prem- bec->Pins - bec->Pct - bec->Pdq - bec->Psw - bec->Pctho) < bec->Piae) {
				  
			steady_clock::time_point t_iae_be= steady_clock::now();
			//bec->insatend();
			steady_clock::time_point t_iae_end = steady_clock::now();
			timestat(6) += static_cast<double>((t_iae_end-t_iae_be).count())* steady_clock::period::num / steady_clock::period::den ;
				 
		  }
				
		  else if ((action- bec->Prem- bec->Pins - bec->Pct - bec->Pdq - bec->Psw - bec->Pctho - bec->Piae) < bec->Prae) {
				  
			steady_clock::time_point t_rae_be = steady_clock::now();
			//bec->rematend();
			steady_clock::time_point t_rae_end = steady_clock::now();
			timestat(7) += static_cast<double>((t_rae_end-t_rae_be).count())* steady_clock::period::num / steady_clock::period::den ;
		  }
#endif
		  else { 
			assert(0);
		  }
			
		  if ((it%bec->Meas_its) ==0) {
			steady_clock::time_point t_meas_be= steady_clock::now();
#ifdef SECUMUL
			meas += bec->measure(ordit);
#else
			meas += bec->measure(0);
#endif
			steady_clock::time_point t_meas_end = steady_clock::now();
			timestat(8) += static_cast<double>((t_meas_end-t_meas_be).count())* steady_clock::period::num / steady_clock::period::den ;
		  }
			
		  if ( (it%bec->Test_its) ==0) {
			steady_clock::time_point t_test_be= steady_clock::now();
			
			try{
			  bec->test();
#ifdef SELFENERGY
			  if (meas < (bec->Test_its/bec->Meas_its) -1) {
				std::cerr << meas << '\t'<< bec->Test_its/bec->Meas_its << std::endl;
				throw meascheck();}
#endif
#ifdef SECUMUL
#endif
			} catch (std::exception& e){
			  std::cerr << e.what() << std::endl;
			  exit(EXIT_FAILURE);
			}
			std::cout << "# Test ok thread " << seed << std::endl;
			
			steady_clock::time_point t_test_end = steady_clock::now();
			timestat(9) += static_cast<double>((t_test_end-t_test_be).count())* steady_clock::period::num / steady_clock::period::den ;
		  }
		  
#ifndef SECUMUL
		  if ( (it%bec->Write_its) ==0) {
			//time 
			steady_clock::time_point time_end = steady_clock::now();
			nseconds = static_cast<double>(duration_cast<seconds>(time_end-time_begin).count());
		  	steady_clock::time_point dt_end = steady_clock::now();
			dt = static_cast<double>(duration_cast<seconds>(dt_end-dt_be).count());
		  
			std::cout << "# " << seed<<  ": Time on core  " << nseconds << '\t' << "dt  " << dt << '\n' << std::endl;
			dt_be = steady_clock::now();

		  }	
		  it++;


		
		} while (nseconds < bec->RunTime - dt);
		
		#pragma omp critical  // to have the same staring point
		{	
		  std::cout << "# Combine thread " << seed << std::endl;
		  Data_tot[seed] = bec->get_Data();
		  Epcum.col(seed) = bec->get_Ep();
		  //Eptest += bec->get_Eptest();
		  uds_tot += bec->get_uds();
		  os_tot += bec->get_os();
		  ts_tot += timestat;
		  qs_tot += bec->get_qs();
	
		  onecore_time += nseconds;
		}
	  }//end of omp for loop
	  
	  
#else	//-------------2nd Part of Order Loop--------------------------

		  if ( (it%bec->Write_its) ==0) {
			//time 
			steady_clock::time_point time_end = steady_clock::now();
			nseconds = static_cast<double>(duration_cast<seconds>(time_end-time_begin).count());
		  	steady_clock::time_point dt_end = steady_clock::now();
			dt = static_cast<double>(duration_cast<seconds>(dt_end-dt_be).count());
		  
			std::cout << "# " << seed<<  ": Time in Order Step  " << nseconds << '\t' << "dt  " << dt << '\n' << std::endl;
			dt_be = steady_clock::now();
			if (bec->normcalc() > bec->normmin-1) {normcheck = true;} // checking if we reached the norm minimum
			if (bec->endcalc() > bec->endmin-1) {endcheck = true;} // checking if we reached the end minimum

			//To check Total Run Time
			steady_clock::time_point Tott_check_end = steady_clock::now();
			Tott_check = static_cast<double>(duration_cast<seconds>(Tott_check_end-Cumt_be).count());
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
		  if (ordit == 0) {Data_tot[seed] = bec->get_Data();}
		  Epcum.col(seed) += bec->get_Ep();
		  uds_tot += bec->get_uds();
		  os_tot += bec->get_os();
		  ts_tot += timestat;
		  qs_tot += bec->get_qs();
	
		  onecore_time += nseconds;
		  
		  //transfering Data
		  SEibtemp.col(seed) = bec->get_SEib();
		  Normstemp.col(seed) = bec->get_NormDiag();
		  Endstemp.col(seed) = bec->get_EndDiag();
		  nnormstemp += bec->normcalc();
		  nendstemp += bec->endcalc();
		  std::cout << "# Transfered Data thread " << seed << std::endl;
		}
		} //end of omp for loop
		
		//Not Parallel anymore
		//Combining Data
		{
		  //Estimator calc
		  Ep_tot.push_back(Epcum.rowwise().mean());
	  	  		  
		  SEib += SEibtemp;
		  
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
		  nnorms.push_back(nnormstemp);
		  nends.push_back(nendstemp);
		  mins.push_back((DiagMC_pl[0]->get_minmax())[0]);
		  maxs.push_back((DiagMC_pl[0]->get_minmax())[1]);
		  std::cout << "# Combined Data!" << std::endl;
		}
			
		// time per  order step
		steady_clock::time_point Cumdt_end = steady_clock::now();
		Cumdt = static_cast<double>(duration_cast<seconds>(Cumdt_end-Cumdt_be).count());
		
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
		
		std::cout << "\n#Total Time " << Cumt << '\t' << "dt  " << Cumdt << '\n' << std::endl;
		
		} while ((Cumt < DiagMC_pl[0]->TotRunTime - Cumdt) && (TotMaxOrdcheck != ordit) );
#endif
	   
	  
	  //Deleting the vector of DiagMC Pointers
	  #pragma omp for ordered schedule(static, 1)
	  for (int seed = 0; seed < nseeds; seed++)
	  {	
		delete DiagMC_pl[seed];
	  }
	
	//Eptest /= static_cast<double>(nseeds);
	uds_tot.col(4)=uds_tot.col(3)/uds_tot.col(1);
	uds_tot.col(5)=uds_tot.col(3)/uds_tot.col(0);
	os_tot/=  static_cast<double>(nseeds);
	ts_tot /= onecore_time;
	qs_tot /= static_cast<double>(nseeds);
	
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
	
	//Estimators	
	ArrayXXd Ep_out(Ep_tot.size(), ws.size()+1);
	Ep_out.col(0) =  Map<const Array<double, 1, Dynamic> > (maxs.data(), mins.size());   //Max Order of Step
	for (int i =0; i < Ep_tot.size() ; i++)
	  Ep_out.rightCols(ws.size()).row(i) = Ep_tot[i];
	
	std::ofstream Ep_outfile("data/Ep/Epvsord");
	Ep_outfile << Ep_out << '\n';
	Ep_outfile.close();
#endif
	
	ArrayXXd Epcum_out(ws.size(),nseeds + 1); 
	Epcum_out.col(0) =  Map<const Array<double, 1, Dynamic> > (ws.data(), ws.size());   //Max Order of Step
	Epcum_out.rightCols(nseeds) = Epcum;
	
	std::ofstream Epcum_outfile("data/Ep/Epvsws");
	Epcum_outfile << Epcum_out << '\n';
	Epcum_outfile.close();
	
	/*
	ArrayXXd Eptest_out(Eptest.rows(), Eptest.cols()+1); 
	Eptest_out.col(0) = Data_tot[0].col(0);   //tau
	Eptest_out.rightCols(Eptest.cols()) = Eptest;
	
	std::ofstream Eptest_outfile("data/Ep/Eptest");
	Eptest_outfile << Eptest_out << '\n';
	Eptest_outfile.close();
	*/
	
	//statistics
	std::ofstream uds_outfile("data/stats/uds_tot"); 
	std::ofstream os_outfile("data/stats/os_tot");
	std::ofstream ts_outfile("data/stats/ts_tot");
	std::ofstream qs_outfile("data/stats/qs_tot");
	
	uds_outfile << "Update Statistics" << '\n';
    uds_outfile << "*************************************" << '\n';
    uds_outfile << "COLUMNS" << '\n' << "1:ATTEMPTED" << '\n' << "2:POSSIBLE" << '\n' << "3:REJECTED" << '\n' << "4:ACCEPTED" << '\n' << "5:ACCEPTANCE RATIO POSSIBLE" << '\n' << "6:ACCEPTANCE RATIO TOTAL" << '\n';
    uds_outfile << "*************************************" << '\n';
    uds_outfile << "CHANGE TAU: " << '\t' << uds_tot.topRows(1) << '\n' << "CT IN HO: " << '\t' << uds_tot.block(7, 0, 1, 6) << '\n' << "INSERT: " << '\t' << uds_tot.block(1, 0, 1, 6) << '\n' << "REMOVE: " << '\t' << uds_tot.block(2, 0, 1, 6) << '\n' << "SWAP:   " << '\t' << uds_tot.block(3, 0, 1, 6) << '\n'<< "SWAPoocc: " << '\t' << uds_tot.block(4, 0, 1, 6) << '\n'<< "SWAPoc: " << '\t' << uds_tot.block(5, 0, 1, 6) << '\n'<< "SWAPco: " << '\t' << uds_tot.block(6, 0, 1, 6) << '\n' << "DQ:      " << '\t' << uds_tot.block(8, 0, 1, 6)<< '\n' << "IAE:      " << '\t' << uds_tot.block(9, 0, 1, 6) << '\n' << "RAE :         " << '\t' << uds_tot.block(10, 0, 1, 6) << '\n' << '\n';
	
	MatrixXd orderprint(os_tot.size(),2);
	orderprint << ArrayXd::LinSpaced(config.get<int>("Order_Stat_Size"), 0, config.get<double>("Order_Stat_Size")-1), os_tot;
	
	os_outfile << "Order Statistics" << '\n';
	os_outfile << orderprint  << '\n';
	
	ts_outfile << "Time Statistics [%]" << std:: endl;
	ts_outfile << "CHANGE TAU:" << '\t' << ts_tot(0) << '\n' << "CT HO:     " << '\t' <<ts_tot(1)<< '\n' <<  "INSERT: " << '\t' << ts_tot(2) << '\n' << "REMOVE: " << '\t' << ts_tot(3) << '\n' <<  "SWAP:     " << '\t' << ts_tot(4) << '\n' << "DQ:      " << '\t' << ts_tot(5) << '\n' << "IAE:      " << '\t' << ts_tot(6) << '\n' << "RAE :       " << '\t' << ts_tot(7) << '\n' << "MEAS:       " << '\t' << ts_tot(8) << '\n' << "TEST:       " << '\t' << ts_tot(9) <<'\n';		
	
	MatrixXd qprint(qs_tot.size(),2);
	qprint << ArrayXd::LinSpaced(config.get<int>("QStat_Size"), 0, config.get<double>("Q_Cutoff")), qs_tot;
	qs_outfile << "Phonon Momentum Statistics" << '\n';
	qs_outfile << qprint  << '\n';
	
	uds_outfile.close();
	os_outfile.close();
	ts_outfile.close();
	qs_outfile.close();
	
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
#endif	
	
	 std::cout << "# Analysing data, result = " << system(config.get<std::string>("Ana_File").c_str()) << std::endl;
  
	return 0;
}



  
