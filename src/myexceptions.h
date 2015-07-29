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










#endif