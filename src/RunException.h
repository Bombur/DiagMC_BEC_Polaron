//Exceptions for DiagMC class

#include <exception>
#include <stdexcept>

#ifndef __RUNEXCEPTION_H_INCLUDED__
#define __RUNEXCEPTION_H_INCLUDED__


#ifdef SELFENERGY
class meascheck: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Not every Diagram is a SE Diagram! Measurement Check failed";
  }
};
#endif


#endif