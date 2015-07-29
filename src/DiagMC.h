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

#ifndef __DIAGMC_H_INCLUDED__
#define __DIAGMC_H_INCLUDED__

namespace pt = boost::property_tree;
using namespace Eigen;

class DiagMC {
  private:
	const double p, mu;				
	const double taumax;			//imaginary time cut off
	const double taubin;
	double E;
	double G0p;
	double tau;
	double alpha; 					//coupling strength
	
	//Calculation variables
	std::function<double()> drnd;
	VectorXd Data;					//0:tau_i, 1:G(p,tau_i)
		
	//Analysis variables
	int binl;						//binning length
	Vector3d anabuffer;				//0:mean, 1:binning error, 2:integrated correlation time
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

	//DiagMC.cpp
	void initialize();							
	void change_tau();		
	//void beta_sweep(double begin, double end, double step, int binning_length);
	
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