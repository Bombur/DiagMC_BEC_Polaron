// Froehlich Polaron
// Hans Peter Guertner
#include <iostream>
#include <string>
#include <ctime>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "DiagMC.h"

#if defined(_OPENMP)
#include <omp.h>
#endif 

namespace pt = boost::property_tree;
using namespace std::chrono;
using namespace Eigen;

class oor_Probs: public std::exception {
  virtual const char* what() const throw()
  {
    return "No possible action chosen";
  }
};

int main() {
  try{
	#if defined(_OPENMP)
	  const int nseeds = 4;
	#else
	  const int nseeds = 1;
	#endif
	  
	std::cout << "# Using " << nseeds << " cores." << std::endl;
	
	//Total Statistics and Data
	std::vector<MatrixXd> Data_tot(nseeds);
	  	
	MatrixXd uds_tot = MatrixXd::Zero(9,6);
	VectorXi os_tot = VectorXi::Zero(40);
	VectorXd ts_tot = VectorXd::Zero(6);
	
	//constants and parameters
	double onecore_time= 0;
	
	pt::ptree config;
	pt::read_json("DiagMC_FP.json", config);
	
	#pragma omp parallel num_threads(nseeds) 
	{

	  #pragma omp for
	  for (int seed = 0; seed < nseeds; seed++)
	  {
		// initialization
		DiagMC fp(seed, config);
		
		#pragma omp critical  // to have the same staring point
		{
		  std::cout << "# Initialize thread " << seed << std::endl;
		}
  
		steady_clock::time_point time_begin = steady_clock::now();  //start time
		VectorXd timestat = VectorXd::Zero(6);
  	
		double nseconds;
		double dt;
		do {
		  std::time_t dt_be;
		  std::time(& dt_be);
		  for (int i=0; i<fp.Write_its; i++) {
			for (int j=0; j<fp.Test_its; j++) {
			  for (int k=0; k<fp.Meas_its; k++) {
				double action = fp.drnd();
			
				if (action < fp.Prem) {
				  /*
				  steady_clock::time_point t_rem_be = steady_clock::now();
				  fp.remove();
				  steady_clock::time_point t_rem_end = steady_clock::now();
				  timestat(3) += static_cast<double>((t_rem_end - t_rem_be).count())* steady_clock::period::num / steady_clock::period::den ;
				  */
				  
				}
				else if ((action-fp.Prem)<fp.Pins) {
				  /*
				  steady_clock::time_point t_ins_be= steady_clock::now();
				  fp.insert();
				  steady_clock::time_point t_ins_end= steady_clock::now();
				  timestat(2) += static_cast<double>((t_ins_end - t_ins_be).count())* steady_clock::period::num / steady_clock::period::den ;
				  */
				  
				}
				else if ((action-fp.Prem-fp.Pins) < fp.Pct) {
				  
				  steady_clock::time_point t_ct_be= steady_clock::now();
				  fp.change_tau();
				  steady_clock::time_point t_ct_end= steady_clock::now();
				  timestat(0) += static_cast<double>((t_ct_end - t_ct_be).count())* steady_clock::period::num / steady_clock::period::den ;
			/*
				  steady_clock::time_point t_ctho_be= steady_clock::now();
				  fp.ct();		//change taus in higher order
				  steady_clock::time_point t_ctho_end= steady_clock::now();
				  timestat(1) += static_cast<double>((t_ctho_end - t_ctho_be).count())* steady_clock::period::num / steady_clock::period::den ;
			  */
				}
			
				else if ((action-fp.Prem-fp.Pins-fp.Pct) < fp.Pdq) {
				  /*
				  steady_clock::time_point t_dq_be= steady_clock::now();
				  fp.dq();
				  steady_clock::time_point t_dq_end= steady_clock::now();
				  timestat(5) += static_cast<double>((t_dq_end - t_dq_be).count())* steady_clock::period::num / steady_clock::period::den ;
				
				  */
				}
			
				else if ((action- fp.Prem- fp.Pins - fp.Pct - fp.Pdq) < fp.Psw) {
				  /*
				  steady_clock::time_point t_sw_be= steady_clock::now();
				  fp.swap();
				  steady_clock::time_point t_sw_end= steady_clock::now();
				  timestat(4) += static_cast<double>((t_sw_end - t_sw_be).count())* steady_clock::period::num / steady_clock::period::den ;
				*/
				  
				}
				
				else { 
				  throw oor_Probs();
				}
			
			
			  }
			  fp.measure(j);
			}
			fp.test();
		  }
		  		  
		  //fp.status();
		  std::cout << "# Test ok thread " << seed << std::endl;
	
		  fp.write();
	
		  steady_clock::time_point time_end = steady_clock::now();
		  steady_clock::duration time_span = time_end-time_begin;
		  nseconds = static_cast<double>(time_span.count()) * steady_clock::period::num / steady_clock::period::den;
		  
		  //fp.timestats(timestat/nseconds);
		  
		  fp.Stattofile(timestat/nseconds);
		  std::time_t dt_end;
		  std::time(& dt_end);
		  
		  dt = std::difftime(dt_end, dt_be);
		  
		  
		  std::cout << "# Time on core  " << nseconds << '\n' << std::endl;
		  
		} while (nseconds < fp.RunTime - dt);
		
		
		
		#pragma omp critical  // to have the same staring point
		{	
		  std::cout << "# Combine thread " << seed << std::endl;
		  Data_tot[seed] = fp.get_Data();
		  uds_tot += fp.get_uds();
		  os_tot += fp.get_os();
		  ts_tot += timestat;
	
		  
		  onecore_time += nseconds;
		}
		
		
	  }
	}
	
	uds_tot.col(4)=uds_tot.col(3).array()/uds_tot.col(1).array();
	uds_tot.col(5)=uds_tot.col(3).array()/uds_tot.col(0).array();
	ts_tot /= onecore_time;
	
	//averages	
	//all
	ArrayXXd all_orders(Data_tot[0].rows(), 3);  // 0:tau, 1:average, 2:error
	all_orders << Data_tot[0].col(0).array() , ArrayXXd::Zero(Data_tot[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  all_orders.col(1) += Data_tot[i].col(1).array();
	all_orders.col(1) /= nseeds;
	
	//zero order
	ArrayXXd zero_order(Data_tot[0].rows(), 3);  // 0:tau, 1:average, 2:error
	zero_order << Data_tot[0].col(0) , ArrayXXd::Zero(Data_tot[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  zero_order.col(1) += Data_tot[i].col(2).array();
	zero_order.col(1) /= nseeds;
	
	//first order
	ArrayXXd first_order(Data_tot[0].rows(), 3);  // 0:tau, 1:average, 2:error
	first_order << Data_tot[0].col(0) , ArrayXXd::Zero(Data_tot[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  first_order.col(1) += Data_tot[i].col(3).array();
	first_order.col(1) /= nseeds;
	
	//second order
	ArrayXXd second_order(Data_tot[0].rows(), 3);  // 0:tau, 1:average, 2:error
	second_order << Data_tot[0].col(0) , ArrayXXd::Zero(Data_tot[0].rows(), 2);   //tau
	for (int i =0; i < nseeds ; i++)
	  second_order.col(1) += Data_tot[i].col(4).array();
	second_order.col(1) /= nseeds;
	
	
	//errors
	if (nseeds > 3) {
	  for (int i=0; i<nseeds ; i++) {
		all_orders.col(2) += (Data_tot[i].col(1).array()-all_orders.col(1)).pow(2); 
		zero_order.col(2) += (Data_tot[i].col(2).array()-zero_order.col(1)).pow(2); 
		first_order.col(2) += (Data_tot[i].col(3).array()-first_order.col(1)).pow(2); 
		second_order.col(2) += (Data_tot[i].col(4).array()-second_order.col(1)).pow(2); 
	  }
	  double norm = 1/(nseeds*(nseeds-1));
	  all_orders.col(2) *= norm;
	  zero_order.col(2) *= norm;
	  first_order.col(2) *= norm;
	  second_order.col(2) *= norm;
	  all_orders.col(2).cwiseSqrt();
	  zero_order.col(2).cwiseSqrt();
	  first_order.col(2).cwiseSqrt();
	  second_order.col(2).cwiseSqrt();
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
	
	std::ofstream second("data/second_order");
	second << second_order << '\n';
	second.close();
	
	
	//statistics
	std::ofstream uds_outfile("data/uds_tot");
	std::ofstream os_outfile("data/os_tot");
	std::ofstream ts_outfile("data/ts_tot");
	
	uds_outfile << "Update Statistics" << '\n';
    uds_outfile << "*************************************" << '\n';
    uds_outfile << "COLUMNS" << '\n' << "1:ATTEMPTED" << '\n' << "2:POSSIBLE" << '\n' << "3:REJECTED" << '\n' << "4:ACCEPTED" << '\n' << "5:ACCEPTANCE RATIO POSSIBLE" << '\n' << "6:ACCEPTANCE RATIO TOTAL" << '\n';
    uds_outfile << "*************************************" << '\n';
    uds_outfile << "CHANGE TAU: " << '\t' << uds_tot.topRows(1) << '\n' << "CT IN HO: " << '\t' << uds_tot.block(7, 0, 1, 6) << '\n' << "INSERT: " << '\t' << uds_tot.block(1, 0, 1, 6) << '\n' << "REMOVE: " << '\t' << uds_tot.block(2, 0, 1, 6) << '\n' << "SWAP:   " << '\t' << uds_tot.block(3, 0, 1, 6) << '\n'<< "SWAPoocc: " << '\t' << uds_tot.block(4, 0, 1, 6) << '\n'<< "SWAPoc: " << '\t' << uds_tot.block(5, 0, 1, 6) << '\n'<< "SWAPco: " << '\t' << uds_tot.block(6, 0, 1, 6) << '\n' << "DQ:      " << '\t' << uds_tot.block(8, 0, 1, 6)<< '\n' << '\n';
	
	MatrixXi orderprint(os_tot.size(),2);
	orderprint << VectorXi::LinSpaced(os_tot.size(), 0, os_tot.size()-1), os_tot;
	
	os_outfile << "Order Statistics" << '\n';
	os_outfile << orderprint  << '\n';
	
	ts_outfile << "Time Statistics [%]" << std:: endl;
	ts_outfile << "CHANGE TAU:" << '\t' << ts_tot(0) << '\n' << "CT HO:     " << '\t' <<ts_tot(1)<< '\n' <<  "INSERT: " << '\t' << ts_tot(2) << '\n' << "REMOVE: " << '\t' << ts_tot(3) << '\n' <<  "SWAP:     " << '\t' << ts_tot(4) << '\n' << "DQ:      " << '\t' << ts_tot(5) << '\n';		

	uds_outfile.close();
	os_outfile.close();
	ts_outfile.close();
	
	std::system("rm -f ./data/*_core*");  //removing the data from each core
	
	
	std::cout << "# Analysing data, result = " << system("python3 data_ana.py") << std::endl;
	return 0;
  } catch (std::exception& e){
	std::cerr << e.what() << std::endl;
	exit(EXIT_FAILURE);
  }
}



	
