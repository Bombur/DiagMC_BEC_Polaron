//exceptions
#include <exception>
#include <stdexcept>

//io
#include <iostream>
#include <string>

//math and container
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

#ifndef __DVECTOR_H_INCLUDED__
#define __DVECTOR_H_INCLUDED__

class dnf_vecsize: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Vector sizes do not fit!";
  }
};

class no_propose: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Nothing proposed!";
  }
};

std::vector<double> vadd(const std::vector<double> & vec1, const std::vector<double> & vec2);

std::vector<double> vsub(const std::vector<double> & vec1, const std::vector<double> & vec2);

double vsq(const std::vector<double> & vec1);

#endif

