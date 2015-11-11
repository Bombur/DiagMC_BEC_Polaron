//Exceptions for DiagMC class

#include <exception>
#include <stdexcept>

#ifndef __DIAGRAMEXCEPTION_H_INCLUDED__
#define __DIAGRAMEXCEPTION_H_INCLUDED__


//Diagram_test.cpp
class dnf_diagvec: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Diagram vector sizes do not fit!";
  }
};

class timeserr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Error in Diagram::times!";
  }
};

class phproperr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Error in Diagram::phprop!";
  }
};

class elproperr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Error in Diagram::elprop!";
  }
};

class momerr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Momentum is not conserved!";
  }
};

class qerr: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Q is not conserved!";
  }
};

class propopen: public std::exception
{
  virtual const char* what() const throw()
  {
    return "A propergator is not closed";
  }
};


class swpos: public std::exception
{
  virtual const char* what() const throw()
  {
    return "Number of possible vertices for swap do not fit number of zero loops!";
  }
};

class qoor: public std::exception
{
  virtual const char* what() const throw()
  {
    return "q out of range!";
  }
};












#endif