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
	double tau;
	
	
	//Calculation variables
	MatrixXi Data;					//0:G(p,tau_i), 1:G0(p,tau_i), 2:G1(p,tau_i), 3:G2(p,tau_i)
	Diagram diag;
			
	//Analysis variables
	int binl;						//binning length
	Vector3d anabuffer;				//0:mean, 1:binning error, 2:integrated correlation time
	
	//rows of stats: 0:change tau, 1:insert, 2:remove
	Matrix<double,4,6>  stats; 		//0:attempted, 1:possible, 2:rejected, 3:accepted, 4:acceptance ratio possible, 5:acceptance ratio total
	
	//io variables
	std::string path;
	std::ofstream Datafile;
	std::ofstream Statsfile;
	long long Datafile_pos;
	long long Statfile_pos;
	
  public:
	const double Prem, Pins;		//probabilities to choose remove or insert branch
	
	//Random Function
	std::function<double()> drnd;
	
	//Execution variables
	const int Meas_its, Test_its, Write_its;
	const int RunTime;			//in seconds
	
	//DiagMC_config.cpp
	DiagMC(const pt::ptree &);
	~DiagMC();

	//DiagMC_run.cpp
	void change_tau();
	int insert();
	int remove();
	
	//DiagMC_updates.cpp
	int get_order() {return diag.get_order();};
	void measure(const int & whichmeas);
	
	//DiagMC_estimator.cpp
	double mean(const VectorXd & measdata);
	void bin_ana(const VectorXd & measdata);
	
	//DiagMC_io.cpp
	void write();
	void Stattofile();
	void status();
	void printmean();
	void printstats();
	
	//DiagMC_test.cpp
	void test(); 				
	
};



   
#endif