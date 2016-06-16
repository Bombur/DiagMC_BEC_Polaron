//config
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
 
//system
#include <cstdlib>

//exception
#include "DiagMCException.h"
#include <assert.h>

//io 
#include <iostream>
#include <fstream>
#include <string>

//random
#include <random>
#include "/project/theorie/h/H.Guertner/lib/pcg/include/pcg_random.hpp"
#include <functional>

//time
#include <chrono>
#include <ctime>
#include <ratio>

//math and container
#include "/project/theorie/h/H.Guertner/lib/Eigen/Eigen/Dense"
#include <cmath>
#include <vector> 
#include <array> 
#include <complex>
#include "Diagram.h"
#include "dvector.h"
#include "mystructs.h"
#include "tmap.h"
#include <alps/accumulators.hpp>
 
#ifndef __DIAGMC_H_INCLUDED__
#define __DIAGMC_H_INCLUDED__
 
namespace pt = boost::property_tree;
using namespace Eigen; 
 
class DiagMC { 
  typedef alps::accumulators::accumulator_set accumulators_type;
  
  private: 
	const double p, mu;			//mu and initial p 	
	const double qc;	//Q cut off
	tmap taumap;		//taubin and taumax are in here
	const double alpha; 					//coupling strength
	const double relm; 		// mI/mB relative mass
	double mass; 	//Mass of Impurity 
	const double wmax;  // omega max for G0SEiw
	const int wbin; 
	double E;
	double G0p1arm;				
	double G0p2arms;
	double G0ptmima;				
#ifndef SECUMUL
	const int maxord;				//maximum order
#else
	int maxord;
#endif 

#ifdef FP
	double wp; //Froehlich Polaron Omega
#endif
 
	//Calculation variables
	ArrayXXd Data;					//0:G(p,tau_i), 1:G0(p,tau_i), 2:G1(p,tau_i), 3:G2(p,tau_i)
	Diagram diag;

	//Estimator
	ArrayXd ws;		    	//different  omegas
	ArrayXd Eptmp;			//temporary container for binning
	ArrayXd Epol; 			//Polaron Energy Estimator for ws
	ArrayXXd SE;			//Estimator for Self Energy 0: All_Order, 1:1st, 2: >1
	ArrayXXcd G0SEiw;			//Estimator for Self Energy 0: All_Order, 1:endl1st, 2: >1
	//binning
	accumulators_type counts; //accumulator for norm counts
	accumulators_type Epbin;  //measured each time
	//Measure Ep always after a number of steps
	double last_g0_count;
	ArrayXd last_measured_Epol;
	accumulators_type Ep_intv; //measured every 500 times including count number
	accumulators_type ordesti; //order estimator

	//To Test Error Convergence
	double last_count_g0attau0;
	accumulators_type G0tau0; //Test container for G_0(Tau)
	accumulators_type G0tau0_each;
	
	//Parameter
	const double ctcor;		//maximum correction in ctho
	const double qcor;		//maximum correction in dq
	const double dtins;			// for insert at end
	const double dqins;			// Insert Q Range
	double fw;  			// fake weight (correction for G0) 
	std::vector<double> fwtab;
	const int qsigma; //sigma for Q selection
	const double explim;  // limit for exponential function Overflow and Underflow
	const int doitlim; //limit for iterations in do while loops

	//Cumulative SE Calculation
	int minord; // minimum Order (Norm Order)
	int ordstsz; // order step Size
	int ordstep; // current order step
#ifdef SECUMUL 
	ArrayXd SEib; //Self Energy inbetween minord and maxord
	ArrayXd Norms; // Norm Diagrams for each order step
	ArrayXd Ends; // End Diagrams for each step
	std::vector<double> nnorms; //counts of being in norm Diagrams
	std::vector<double> nends; // counts of being in end Diagrams
	
	//for fw_adapt
	//0:minord, 1:maxord
	const int which_fw_ad; // desired order ratio of maxord/minord 
	std::vector<int> fw_last;
	std::vector<double> fw_counts;	//container to keep track of max min ratio
	const double fw_max;
	const double desi_rat;
	std::vector<double> fw_vec;//vec of size ordstepsize for adaption
#endif
			
	//rows of stats: 0:change tau, 1:insert, 2:remove 3:swap, 4:swapoc, 5:swapco, 6:swapoocc, 7:ct_ho, 8:dq , 9:inscrossed, 10: remcrossed, 11:fo_insert, 12:fo_remove
	ArrayXXd  updatestat; 		//0:attempted, 1:possible, 2:rejected, 3:accepted, 4:acceptance ratio possible, 5:acceptance ratio total  
	//addition row for accumulators 13:accumulator
	ArrayXXd overflowstat;  //0:Overflow, 1:Underflow, 2:rejected, 3:accepted, 4:rejected ratio rejected, 5:acceptance ratio accepted
	ArrayXd orderstat;
	ArrayXXd qstat;		//0:all Orders, 1: 1st_Order
	ArrayXXd tstat;		//Length of Propergator 0:all Orders, 1: 1st reducible, 2:Rainbow, 3:crossed
	ArrayXXd testhisto;

 	//io variables 
	std::string path;
	
	//test variables
	double global_weight; 
	
	//Bools instead of Pragmas
	const bool fo_insrem;
	const bool fog0seg0;
	const bool ct_taucut;
	const bool ct_lin;
	const bool ins_tau2cut;
	const bool ins_tau2lin;
	const bool ic_taulin;
	const bool ic_taucut;
	const bool fst_Ep_meas;

	std::string lu;				//last update for error message
  public: 
	const double Prem, Pins, Pct, Pctho, Psw, Pdq, Pic, Prc;		//probabilities to choose remove or insert branch
	
	//Random Function
	std::function<double()> drnd;
	
	//Execution variables
	const int Meas_its, Test_its, Write_its;
	const int Ep_meas_it;
	const bool Ep_bin_each_step;
	const int RunTime;			//in seconds (in SECUMUL) Time for one order Step
	const int ThermTime;
  
	//Cumulative SE Calculation
	const int normmin; // minimum number of points to do the next step
	const int endmin; // minimum number of points to do the next step
	const int TotRunTime; // Total Time
	const int TotMaxOrd; // Total Maximum Order
#ifdef SECUMUL
	void set_ordstsz(const int & tmp) {ordstsz = tmp;};
	int get_ordstsz() {return ordstsz;};
	void ord_step();  // sets the minimum and maximum order for the next step
	int get_order() {return diag.get_order();}
	int get_max_order() {return maxord;}			//return current maximum order
	double normcalc();  //Calculates how often we have been in the norm Diagram
	double endcalc(); //Calculates how often we have been in the end Diagram
	double get_normostat() {return orderstat(minord);} ; //returns how often we actually have been in the different orders from orderstat
	double get_endostat() {return orderstat(maxord);};
  
	void set_av_nei(const double & av_normi, const double & av_endi, const int & ordit); //Sets the averaged number of all threads to the Norm/End of step ordit
	std::array<double, 2> get_minmax();
	
	//returning Data
	double pref_calc();
	ArrayXd get_SEib();
	ArrayXd get_NormDiag();
	ArrayXd get_EndDiag();
	
	//Fake Weight Adaption
	int get_which_fw_adapt() {return which_fw_ad;}
	int fw_adapt();
	int fw_vec_adapt();
	double get_fw() {return fw;};
	std::vector<double> get_fw_vec() {return fw_vec;}
#endif
	
	//DiagMC_config.cpp
	DiagMC(const int &, const pt::ptree &);
	~DiagMC();
	accumulators_type create_empty_Ep_acc();
	accumulators_type create_empty_count_acc(const int &);
	accumulators_type create_empty_ordesti(const int &);
	
	//DiagMC_updates.cpp
	double G0el(const std::array< double, 3 > & p, const double & tfin, const double & tinit);
	double Dph(const std::array< double, 3 > & q, const double & tfin, const double & tinit);
	double Vq2(const std::array< double, 3 > & q);
	double logG0el(const std::array< double, 3 > & p, const double & tfin, const double & tinit);
	double logDph(const std::array< double, 3 > & q, const double & tfin, const double & tinit);
	double logVq2(const std::array< double, 3 > & q);
	double disprel(const std::array< double, 3 > & q); //dispersion relation
	double irwsinV(const std::array< double, 3 > & p, const std::array< double, 3 > & q, const double & tfin, const double & tinit); //Insert Remove weigth without V(q)^2
	double logirwsinV(const std::array< double, 3 > & p, const std::array< double, 3 > & q, const double & tfin, const double & tinit);
	double tau2pref(const std::array< double, 3 > & p, const std::array< double, 3 > & q); //prefactor of tau2
	int change_tau();   //only 0 order
	int ct();
	int insert();
	int remove();
	int swap();
	int dq();
	int inscrossed();
	int remcrossed();
	int fo_insert();  //first Order Insert
	int fo_remove();  //first Order Insert
	
	//DiagMC_measure.cpp
	double fw_for_meas(const int & ord, const bool & fw_ad=false);
	int measure();
	void meas_histo();
	void meas_ordstat();
	void set_os_to_zero();
	
	//DiagMC_estimator.cpp
	void fill_Epcont(const double & cor); 
	void meas_Epol();
	void meas_SE(const double & cor); 
	void meas_G0SEiw(const double & cor); 
	//binning
	void bin_Epol();
	void meas_Ep_intv(); // measure just evere Ep_meas_it steps
	void meas_G0tau0();
	void meas_G0tau0_each();
	void meas_ordesti();
	
	//DiagMC_io.cpp
	ArrayXXd get_Data();
	ArrayXXd get_Eptest(); 
	ArrayXd get_Ep();		//return Polaron Energy depending omegas (per Step for SECUMUL)
	ArrayXXd get_SE();
	ArrayXXcd get_G0SEiw();
	//binning
	accumulators_type get_counts(); //return accumulator for norm counts
	accumulators_type get_Epacc(); //return accumulator measuring each time
	double get_counts_for_Ep_norm(); //to calculate Norm later
	double get_Ep_norm(); 
	accumulators_type get_Ep_intv(); //return accumulator measuring sum of Ep_meas_it
	accumulators_type get_G0tau0();
	alps::accumulators::result_set get_G0tau0_each();
	void get_ordesti(const int &);//Write each core seperat

	ArrayXXd get_uds() {return updatestat;}
	ArrayXXd get_ofs();
	ArrayXd get_os() {return orderstat;}
	ArrayXXd get_qs() {return qstat;}
	ArrayXXd get_taus();
	ArrayXXd get_testhisto();
  
	
	//DiagMC_test.cpp
	void test(); 
	double weight_calc();
	void printall();
	void printdiag();
};

   
#endif
