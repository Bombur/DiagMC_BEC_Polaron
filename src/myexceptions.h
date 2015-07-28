#include <exception>
#include <stdexcept>

#ifndef __MYEXCEPTIONS_H_INCLUDED__
#define __MYEXCEPTIONS_H_INCLUDED__

using namespace std;

//ios exceptions
class openwritefile: public exception
{
  virtual const char* what() const throw()
  {
    return "Write File not open!";
  }
};


class filepos: public exception
{
  virtual const char* what() const throw()
  {
    return "Positions in Write File do not fit!";
  }
};




//tests
class out_of_range_E: public exception
{
  virtual const char* what() const throw()
  {
    return "The mean energy is out of range";
  }
};

class : public exception
{
  virtual const char* what() const throw()
  {
    return "The mean magnetization is less than 0!";
  }
};

class Weighterror: public exception
{
  virtual const char* what() const throw()
  {
    return "The weight and acceptance ratio do not fit!";
  }
};






#endif