//Exceptions for DiagMC class

#include <exception>
#include <stdexcept>

#ifndef __DIAGMCEXCEPTION_H_INCLUDED__
#define __DIAGMCEXCEPTION_H_INCLUDED__


//DiagMC_test.cpp
#ifdef SECUMUL
class ossoor: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Order Step Size out of range!";
  }
};

#endif

class oor_Gp: public std::exception
{
  virtual const char* what() const throw()
  {
    return "G(p) is out of range!";
  }
};

class oor_tau: public std::exception
{
  virtual const char* what() const throw()
  {
    return "tau is out of range!";
  }
};

class data_empty: public std::exception
{
  virtual const char* what() const throw()
  {
    return "G(p, tau_i) is empty!";
  }
};

class greenerr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Error in the Data of the Green's function!";
  }
};

class weight_check: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Weight Check failed!";
  }
};

class fake_check: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Fake Check failed!";
  }
};

class staterr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Statistics do not match!";
  }
};

class ins_rem: public std::exception
{
  virtual const char* what() const throw()
  {
    return "More Removes than Inserts!";
  }
};

class maxoerr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Maximum Order does not hold!";
  }
};

//DiagMC_config.cpp 
class oor_Probs: public std::exception {
  virtual const char* what() const throw()
  {
    return "The Probabilities do not add up to 1!";
  }
};  

//DiagMC_io.cpp 
class openwritefile: public std::exception {
  virtual const char* what() const throw()
  {
    return "Write File not open!";
  }
};

class filepos: public std::exception {
  virtual const char* what() const throw()
  {
    return "Positions in Write File do not fit!";
  }
};


class Therm: public std::exception {
  virtual const char* what() const throw()
  {
    return "Thermalization time is bigger than RunTime!";
  }
};


class Red: public std::exception {
  virtual const char* what() const throw()
  {
    return "The diagramm is reducible!";
  }
};

class oor_g0sum: public std::exception {
  virtual const char* what() const throw()
  {
    return "Data.sum() and Norms.sum() do not fit!";
  }
};



#endif
