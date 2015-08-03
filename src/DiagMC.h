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
#include "myexceptions.h"

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
	const double taubin;			//bins für 0 bis tau
	const double Prem, Pins;		//probabilities to choose remove or insert branch
	double E;
	double G0p;				//Green's Function G0(p)
	double tau;
	double alpha; 					//coupling strength

	
	//Calculation variables
	std::function<double()> drnd;
	VectorXd Data;					//0:tau_i, 1:G(p,tau_i)
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
	//Execution variables
	const int Meas_its, Test_its, Write_its;
	const int RunTime;					//in seconds
	
	//DiagMC_config.cpp
	DiagMC(const pt::ptree &);
	~DiagMC();

	//DiagMC_run.cpp
	void initialize();							
	void change_tau();
	void insert();
	void remove();
	
	//DiagMC_updates.cpp
	void measure(const int & whichmeas);
	
	//DiagMC_estimator.cpp
	void status();
	double mean(const VectorXd & measdata);
	void bin_ana(const VectorXd & measdata);
	
	//DiagMC_io.cpp
	void write();
	void Stattofile();
	//void meantofile(std::ofstream & file, const long long writing_pos);
	void printmean();
	void printstats();
	
	//Test functions
	void test(const int &); 				//Ising.cpp

	void printallvariables();	//Ising_io.cpp
	
};



   
#endif