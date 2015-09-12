//config
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

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

//exceptions
#include <exception>
#include <stdexcept>
 
//math and container
#include "/project/theorie/h/H.Guertner/lib/Eigen/Eigen/Dense"
#include <cmath>
#include "Diagram.h"
#include "dvector.h"
 

#ifndef __DIAGMC_H_INCLUDED__
#define __DIAGMC_H_INCLUDED__

namespace pt = boost::property_tree;
using namespace Eigen;

class DiagMC {
  private:
	const double p, mu;			//mu and initial p 		
	const double taumax;			//imaginary time cut off
	const int taubin;			//bins für 0 bis tau
	const double alpha; 					//coupling strength
	const double omegap;
	double E;
	double G0p;				//Green's Function G0(p)
	
	
	//Calculation variables
	MatrixXi Data;					//0:G(p,tau_i), 1:G0(p,tau_i), 2:G1(p,tau_i), 3:G2(p,tau_i)
	Diagram diag;
	const double ctcor;		//maximum correction in ctho
	const double qcor;		//maximum correction in ctho
	
			
	//rows of stats: 0:change tau, 1:insert, 2:remove 3:swap, 4:swapoocc, 5:swapoc, 6:swapco, 7:ct_ho, 8:dq
	MatrixXd  updatestat; 		//0:attempted, 1:possible, 2:rejected, 3:accepted, 4:acceptance ratio possible, 5:acceptance ratio total
	VectorXi orderstat;
	
	//io variables
	std::string path;
	std::ofstream Datafile;
	std::ofstream udsfile;			//update statistcs
	std::ofstream osfile;			//order statistics
	std::ofstream tsfile;			//time statistics
	long long Datafile_pos;
	long long uds_pos;
	long long os_pos;
	long long ts_pos;
	
  public:
	const double Prem, Pins, Pct, Psw, Pdq;		//probabilities to choose remove or insert branch
	
	//Random Function
	std::function<double()> drnd;
	
	//Execution variables
	const int Meas_its, Test_its, Write_its;
	const int RunTime;			//in seconds
	
	//DiagMC_config.cpp
	DiagMC(const pt::ptree &);
	~DiagMC();

	//DiagMC.cpp
	int change_tau();   //only 0 order
	int ct();
	int insert();
	int remove();
	int swap();
	int dq();
	
	//DiagMC_updates.cpp
	int get_order() {return diag.get_order();};
	void measure(const int & whichmeas);
	
	//DiagMC_estimator.cpp
	double G0el(const std::vector< double > & p, const double & tfin, const double & tinit);
	double Dph(const double & tfin, const double & tinit);
	
	//DiagMC_io.cpp
	void write();
	void Stattofile(const VectorXd &);
	
		
	void status();
	void updatestat();
	void orderstats();
	void timestats(const VectorXd &);
	
	//for multible cores
	const MatrixXd & get_Data();
	const MatrixXd & get_uds() {return updatestat;}
	const MatrixXi & get_os() {return orderstat;}
	
	
	//DiagMC_test.cpp
	void test(); 				
	
};



   
#endif