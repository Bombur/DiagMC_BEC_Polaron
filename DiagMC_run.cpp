// Froehlich Polaron
// Hans Peter Guertner
#include <iostream>
#include <string>
#include <ctime>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "DiagMC.h"
#include "RunException.h"

#if defined(_OPENMP)
#include <omp.h>
#endif  

namespace pt = boost::property_tree;
using namespace std::chrono;
using namespace Eigen;


int main() {
	#if defined(_OPENMP)
	  const int nseeds = omp_get_max_threads();
	  //const int nseeds = 1;
	#else
	  const int nseeds = 1;
	#endif
	  
	std::cout << "# Using " << nseeds << " cores." << std::endl;
	
	pt::ptree config;
	pt::read_json("DiagMC_BEC.json", config);
	
	//Total Statistics and Data
	std::vector<ArrayXXd> Data_tot(nseeds);
	  	
	ArrayXXd uds_tot = ArrayXXd::Zero(11,6);
	ArrayXi os_tot = ArrayXi::Zero(40);
	ArrayXd ts_tot = ArrayXd::Zero(8);
	
	//Statistics and Data for SECUMUL
#ifdef SECUMUL
	ArrayXXd SEib = ArrayXXd::Zero(config.get<int>("Tau_bin"),nseeds);
	std::vector<ArrayXXd> Norms;
	std::vector<ArrayXXd> Ends;
	std::vector<int> nnorms;
	std::vector<int> nends;
	std::vector<std::vector<int>> minmax;
	
	//temporary
	ArrayXXd SEibtemp(config.get<int>("Tau_bin"),nseeds);
	ArrayXXd Normstemp(config.get<int>("Tau_bin"),nseeds);
	ArrayXXd Endstemp(config.get<int>("Tau_bin"),nseeds);
	int nnormstemp;
	int nendstemp;
		
	double ordstsz_calc=0.;
#endif
	
	//constants and parameters
	double onecore_time= 0;
		
	#pragma omp parallel num_threads(nseeds) 
	{

	  #pragma omp for ordered schedule(static, 1)
	  for (int seed = 0; seed < nseeds; seed++)
	  {
		// initialization
		DiagMC fp(seed, config);
		
		#pragma omp critical  // to have the same staring point
		{
		  std::cout << "# Initialize thread " << seed << std::endl;
		}
		
#ifdef SECUMUL
		//Total time Control
		steady_clock::time_point Cumt_be = steady_clock::now();  //start time
		steady_clock::time_point Cumdt_be = steady_clock::now();
  		double Cumt =0.;
		double Cumdt= 0.;
		
		//order steps iterator
		int ordit =0;
		
		do{
		bool maxordcheck = true;
		bool normcheck = true;
#endif
		//time control
		steady_clock::time_point time_begin = steady_clock::now();  //start time
		ArrayXd timestat = ArrayXd::Zero(8);
		steady_clock::time_point dt_be = steady_clock::now();
  	
		double nseconds =0.;
		double dt= 0.;
		
		//measure control
		int meas =0;
		
		//iterator
		int it= 1;
		do {
		  double action = fp.drnd();
			
		  if (action < fp.Prem) {
			steady_clock::time_point t_rem_be = steady_clock::now();
			fp.remove();
			steady_clock::time_point t_rem_end = steady_clock::now();
			timestat(3) += static_cast<double>((t_rem_end-t_rem_be).count())* steady_clock::period::num / steady_clock::period::den ;
							
		  }
		  else if ((action-fp.Prem)<fp.Pins) {
			steady_clock::time_point t_ins_be= steady_clock::now();
			fp.insert();
			steady_clock::time_point t_ins_end = steady_clock::now();
			timestat(2) += static_cast<double>((t_ins_end-t_ins_be).count())* steady_clock::period::num / steady_clock::period::den ;
		  }
		  else if ((action-fp.Prem-fp.Pins) < fp.Pct) {
			
			steady_clock::time_point t_ct_be= steady_clock::now();
			fp.change_tau();
			steady_clock::time_point t_ct_end = steady_clock::now();
			timestat(0) += static_cast<double>((t_ct_end-t_ct_be).count())* steady_clock::period::num / steady_clock::period::den ;
		  }
			
		  else if ((action-fp.Prem-fp.Pins-fp.Pct) < fp.Pdq) {
			steady_clock::time_point t_dq_be= steady_clock::now();
			fp.dq();
			steady_clock::time_point t_dq_end = steady_clock::now();
			timestat(5) += static_cast<double>((t_dq_end-t_dq_be).count())* steady_clock::period::num / steady_clock::period::den ;
				
		  }
			
		  else if ((action- fp.Prem- fp.Pins - fp.Pct - fp.Pdq) < fp.Psw) {
				  
			steady_clock::time_point t_sw_be= steady_clock::now();
			fp.swap();
			steady_clock::time_point t_sw_end = steady_clock::now();
			timestat(4) += static_cast<double>((t_sw_end-t_sw_be).count())* steady_clock::period::num / steady_clock::period::den ;
				 
		  }
				
		  else if ((action- fp.Prem- fp.Pins - fp.Pct - fp.Pdq- fp.Psw) < fp.Pctho) {
				  
			steady_clock::time_point t_ctho_be= steady_clock::now();
			fp.ct();
			steady_clock::time_point t_ctho_end = steady_clock::now();
			timestat(1) += static_cast<double>((t_ctho_end-t_ctho_be).count())* steady_clock::period::num / steady_clock::period::den ;
				 
		  }

#ifdef FP
		  else if ((action- fp.Prem- fp.Pins - fp.Pct - fp.Pdq - fp.Psw - fp.Pctho) < fp.Piae) {
				  
			steady_clock::time_point t_iae_be= steady_clock::now();
			//fp.insatend();
			steady_clock::time_point t_iae_end = steady_clock::now();
			timestat(6) += static_cast<double>((t_iae_end-t_iae_be).count())* steady_clock::period::num / steady_clock::period::den ;
				 
		  }
				
		  else if ((action- fp.Prem- fp.Pins - fp.Pct - fp.Pdq - fp.Psw - fp.Pctho - fp.Piae) < fp.Prae) {
				  
			steady_clock::time_point t_rae_be = steady_clock::now();
			//fp.rematend();
			steady_clock::time_point t_rae_end = steady_clock::now();
			timestat(7) += static_cast<double>((t_rae_end-t_rae_be).count())* steady_clock::period::num / steady_clock::period::den ;
		  }
#endif
		  else { 
			assert(0);
		  }
			
		  if ((it%fp.Meas_its) ==0) {
#ifdef SECUMUL
			meas += fp.measure(ordit);
#else
			meas += fp.measure(0);
#endif
		  }
			
		  if ( (it%fp.Test_its) ==0) {
			try{
			  fp.test();
#ifdef SELFENERGY
			  if (meas < fp.Meas_its) {throw meascheck();}
#endif
#ifdef SECUMUL
#endif
			} catch (std::exception& e){
			  std::cerr << e.what() << std::endl;
			  exit(EXIT_FAILURE);
			}
			std::cout << "# Test ok thread " << seed << std::endl;
		  }
		  
#ifndef SECUMUL
		  if ( (it%fp.Write_its) ==0) {
			//fp.write();
			//fp.timestats(timestat/nseconds);
		  	//fp.Stattofile(timestat/nseconds);
			
			//time 
			steady_clock::time_point time_end = steady_clock::now();
			nseconds = duration_cast<seconds>(time_end-time_begin).count();
		  	steady_clock::time_point dt_end = steady_clock::now();
			dt = duration_cast<seconds>(dt_end-dt_be).count();
		  
			std::cout << "# " << seed<<  ": Time on core  " << nseconds << '\t' << "dt  " << dt << '\n' << std::endl;
			dt_be = steady_clock::now();

		  }	
		  it++;


		
		} while (nseconds < fp.RunTime - dt);
		fp.printdiag();
		
		#pragma omp critical  // to have the same staring point
		{	
		  std::cout << "# Combine thread " << seed << std::endl;
		  Data_tot[seed] = fp.get_Data();
		  uds_tot += fp.get_uds();
		  os_tot += fp.get_os();
		  ts_tot += timestat;
	
		  onecore_time += nseconds;
		}
		
#else		//rest of write()
		  if ( (it%fp.Write_its) ==0) {
			//fp.write();
			//fp.timestats(timestat/nseconds);
		  	//fp.Stattofile(timestat/nseconds);
			
			//time 
			steady_clock::time_point time_end = steady_clock::now();
			nseconds = duration_cast<seconds>(time_end-time_begin).count();
		  	steady_clock::time_point dt_end = steady_clock::now();
			dt = duration_cast<seconds>(dt_end-dt_be).count();
		  
			std::cout << "# " << seed<<  ": Time on core  " << nseconds << '\t' << "dt  " << dt << '\n' << std::endl;
			dt_be = steady_clock::now();
			if (fp.normcalc() > fp.normmin-1) {normcheck = false;} // checking if we reached the norm minimum
		  }	
		  it++;
			
		  //to go from one order step to the next 
		  //the diagram has to be in the maximum order of before
		  if (nseconds > fp.RunTime - dt && !(normcheck)) {
			if(fp.get_order()== fp.get_max_order()) {maxordcheck = false;}
		  }
		} while (nseconds < fp.RunTime - dt || normcheck || maxordcheck);
		
		
		#pragma omp ordered  // to have the same staring point
		{	
		  std::cout << "# Combine thread " << seed << std::endl;
		  if (ordit == 0) {Data_tot[seed] = fp.get_Data();}
		  uds_tot += fp.get_uds();
		  os_tot += fp.get_os();
		  ts_tot += timestat;
	
		  onecore_time += nseconds;
		  
		  SEibtemp.col(seed) = fp.get_SEib();
		  Normstemp.col(seed) = fp.get_NormDiag();
		  Endstemp.col(seed) = fp.get_EndDiag();
		  nnormstemp += fp.normcalc();
		  nendstemp += fp.endcalc();		  
		}
		
		//transfering Data (needs to be done just once => single
		if (seed == nseeds - 1) {
		  SEib += SEibtemp;
		  Norms.push_back(Normstemp);
		  Ends.push_back(Endstemp);
		  nnorms.push_back(nnormstemp);
		  nends.push_back(nendstemp);
		  minmax.push_back(fp.get_minmax());
		}
			
		// time per  order step
		steady_clock::time_point Cumdt_end = steady_clock::now();
		Cumdt = duration_cast<seconds>(Cumdt_end-Cumdt_be).count();
		/*
		#pragma omp ordered  
		{
		  ordstsz_calc += fp.check_ordstsz(Cumdt)/nseeds;
		}
		
		#pragma omp barrier
		fp.set_ordstsz(ordstsz_calc);
		*/
		fp.ord_step();
		ordit++;
		std::cout << "# Order Step thread " << seed << std::endl;
		
		/*
		#pragma omp barrier
		#pragma omp single nowait
		{ordstsz_calc=0.;}
		*/
		//total time 
		steady_clock::time_point Cumt_end = steady_clock::now();
		Cumt = duration_cast<seconds>(Cumt_end-Cumt_be).count();
		
		} while (Cumt < fp.TotRunTime - Cumdt);
#endif
		
		
	  }
	}
	
	uds_tot.col(4)=uds_tot.col(3)/uds_tot.col(1);
	uds_tot.col(5)=uds_tot.col(3)/uds_tot.col(0);
	ts_tot /= onecore_time;
	
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
		Norm_last.col(2) += (Norms[Norms.size()-1].col(i)-Norm_first.col(1)).pow(2);
		End_last.col(2) += (Ends[Ends.size()-2].col(i)-End_first.col(1)).pow(2);
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
	std::ofstream Nf("data/Norm_first");  
	Nf << Norm_first << '\n';
	Nf.close();
	
	std::ofstream Ef("data/End_first");  
	Ef << End_first << '\n';
	Ef.close();
	
	std::ofstream Nl("data/Norm_last");  
	Nl << Norm_last << '\n';
	Nl.close();
	
	std::ofstream El("data/End_last");  
	El << End_last << '\n';
	El.close();
#endif	
	
	
	//statistics
	std::ofstream uds_outfile("data/uds_tot");
	std::ofstream os_outfile("data/os_tot");
	std::ofstream ts_outfile("data/ts_tot");
	
	uds_outfile << "Update Statistics" << '\n';
    uds_outfile << "*************************************" << '\n';
    uds_outfile << "COLUMNS" << '\n' << "1:ATTEMPTED" << '\n' << "2:POSSIBLE" << '\n' << "3:REJECTED" << '\n' << "4:ACCEPTED" << '\n' << "5:ACCEPTANCE RATIO POSSIBLE" << '\n' << "6:ACCEPTANCE RATIO TOTAL" << '\n';
    uds_outfile << "*************************************" << '\n';
    uds_outfile << "CHANGE TAU: " << '\t' << uds_tot.topRows(1) << '\n' << "CT IN HO: " << '\t' << uds_tot.block(7, 0, 1, 6) << '\n' << "INSERT: " << '\t' << uds_tot.block(1, 0, 1, 6) << '\n' << "REMOVE: " << '\t' << uds_tot.block(2, 0, 1, 6) << '\n' << "SWAP:   " << '\t' << uds_tot.block(3, 0, 1, 6) << '\n'<< "SWAPoocc: " << '\t' << uds_tot.block(4, 0, 1, 6) << '\n'<< "SWAPoc: " << '\t' << uds_tot.block(5, 0, 1, 6) << '\n'<< "SWAPco: " << '\t' << uds_tot.block(6, 0, 1, 6) << '\n' << "DQ:      " << '\t' << uds_tot.block(8, 0, 1, 6)<< '\n' << "IAE:      " << '\t' << uds_tot.block(9, 0, 1, 6) << '\n' << "RAE :         " << '\t' << uds_tot.block(10, 0, 1, 6) << '\n' << '\n';
	
	MatrixXi orderprint(os_tot.size(),2);
	orderprint << VectorXi::LinSpaced(os_tot.size(), 0, os_tot.size()-1), os_tot;
	
	os_outfile << "Order Statistics" << '\n';
	os_outfile << orderprint  << '\n';
	
	ts_outfile << "Time Statistics [%]" << std:: endl;
	ts_outfile << "CHANGE TAU:" << '\t' << ts_tot(0) << '\n' << "CT HO:     " << '\t' <<ts_tot(1)<< '\n' <<  "INSERT: " << '\t' << ts_tot(2) << '\n' << "REMOVE: " << '\t' << ts_tot(3) << '\n' <<  "SWAP:     " << '\t' << ts_tot(4) << '\n' << "DQ:      " << '\t' << ts_tot(5) << '\n' << "IAE:      " << '\t' << ts_tot(6) << '\n' << "RAE :       " << '\t' << ts_tot(7) << '\n';		

	uds_outfile.close();
	os_outfile.close();
	ts_outfile.close();
	
	int rem=std::system("rm -f ./data/*_core_*");  //removing the data from each core
	
	const std::string anacommand = "python3 data_ana.py " + config.get<std::string>("Peter_Path");
	std::cout << "# Analysing data, result = " << system(anacommand.c_str()) << std::endl;
	return 0;
}



	
