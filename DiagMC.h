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
#include "Diagram.h"
#include "dvector.h"
#include "mystructs.h"
#include "tmap.h"

#ifndef __DIAGMC_H_INCLUDED__
#define __DIAGMC_H_INCLUDED__
 
namespace pt = boost::property_tree;
using namespace Eigen; 

class DiagMC { 
  private: 
	const double p, mu;			//mu and initial p 	
	const double qc;	//Q cut off
	tmap taumap;		//taubin and taumax are in here
	const double alpha; 					//coupling strength
	const double relm; 		// mI/mB relative mass
	double E;
	double G0p;				//Green's Function G0(p)
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
	ArrayXXd Epol; 			//Polaron Energy

	//Parameter
	const double ctcor;		//maximum correction in ctho
	const double qcor;		//maximum correction in dq
	const double dtins;			// for insert at end
	const double dqins;			// for insert at endend
	const double fw;  			// fake weight (correction for G0) 
	const double fwzero; 		//fake weight for zero order
	const double fwone;		//fake weight for order <=1	
	const double sigfac;		//factor for sigma of Gauss distribution for choosing q

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
#endif
			
	//rows of stats: 0:change tau, 1:insert, 2:remove 3:swap, 4:swapoocc, 5:swapoc, 6:swapco, 7:ct_ho, 8:dq , 9:insatend, 10: rematend
	ArrayXXd  updatestat; 		//0:attempted, 1:possible, 2:rejected, 3:accepted, 4:acceptance ratio possible, 5:acceptance ratio total
	ArrayXd orderstat;
	ArrayXd qstat;
  
	//io variables 
	std::string path;
	
	//test variables
	double global_weight; 
	std::string lu;				//last update for error message
	
  public: 
	const double Prem, Pins, Pct, Pctho, Psw, Pdq, Piae, Prae;		//probabilities to choose remove or insert branch
	
	//Random Function
	std::function<double()> drnd;
	
	//Execution variables
	const int Meas_its, Test_its, Write_its;
	const int RunTime;			//in seconds (in SECUMUL) Time for one order Step
	
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
	void set_av_nei(const double & av_normi, const double & av_endi, const int & ordit); //Sets the averaged number of all threads to the Norm/End of step ordit
	std::array<double, 2> get_minmax();
	
	//returning Data
	double pref_calc();
	ArrayXd get_SEib();
	ArrayXd get_NormDiag();
	ArrayXd get_EndDiag();
#endif
	
	//DiagMC_config.cpp
	DiagMC(const int &, const pt::ptree &);
	~DiagMC();
	
	//DiagMC_updates.cpp
	int change_tau();   //only 0 order
	int ct();
	int insert();
	int remove();
	int swap();
	int dq();
	int insatend();
	int rematend();
	int measure(const int & ordstp);
	
	//DiagMC_estimator.cpp
	double G0el(const std::array< double, 3 > & p, const double & tfin, const double & tinit);
	double Dph(const std::array< double, 3 > & q, const double & tfin, const double & tinit);
	double Vq2(const std::array< double, 3 > & q);
	void meas_Epol(const double & cor);
	ArrayXXd get_Eptest(); 
	ArrayXd get_Ep();		//return Polaron Energy depending omegas (per Step for SECUMUL)
	
	//DiagMC_io.cpp
	void write();
	void Stattofile(const ArrayXd &);
	
	//for multible cores
	ArrayXXd get_Data();
	ArrayXXd get_uds() {return updatestat;}
	ArrayXd get_os() {return orderstat;}
	ArrayXd get_qs() {return qstat;}
	
	 
	//DiagMC_test.cpp
	void test(); 
	double weight_calc();
	void printall();
	void printdiag();
};



   
#endif
