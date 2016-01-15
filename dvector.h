//exceptions 
#include <assert.h>

//io
#include <iostream>
#include <fstream>
#include <string>

//math and container
#include <vector>
#include <array>
#include <cmath>

//template
#include <type_traits>

#ifndef __DVECTOR_H_INCLUDED__
#define __DVECTOR_H_INCLUDED__

//is_container()

template <typename T1, std::size_t N = 0> struct is_container : std::false_type {};
template <typename T1> struct is_container<std::vector<T1> > : std::true_type { };
template <typename T1, std::size_t N> struct is_container<std::array<T1, N> > : std::true_type { };

//is_number
template<class T> struct is_number : std::integral_constant<bool, std::is_integral<T>::value ||std::is_floating_point<T>::value>::type {};

template <typename vecT, typename std::enable_if<is_container<vecT>::value, int>::type = 0>
vecT operator+=(vecT & lhs, const vecT & rhs){
  assert(lhs.size()== rhs.size());
  int it=0;
  for (auto & x : lhs ){
 	  x += rhs[it];
	  it++;
  }
  return lhs;
}

template <typename vecT, typename std::enable_if<is_container<vecT>::value, int>::type = 0>
vecT operator+(vecT lhs, const vecT & rhs){
  return lhs += rhs;
}

template <typename vecT, typename std::enable_if<is_container<vecT>::value, int>::type = 0>
vecT operator-=(vecT & lhs, const vecT & rhs){
  assert(lhs.size()== rhs.size());
  int it=0;
  for (auto & x : lhs ){
 	  x -= rhs[it];
	  it++;
  }
  return lhs;
}
 
template <typename vecT, typename std::enable_if<is_container<vecT>::value, int>::type = 0>
vecT operator-(vecT lhs, const vecT & rhs){
  return lhs -= rhs;
}


//product with scalar
 
template <typename vecT, typename valT, typename std::enable_if<is_container<vecT>::value && is_number<valT>::value, int>::type = 0>
vecT operator*=(vecT & lhs, const valT & rhs){
  for (auto &x : lhs ){
 	  x *= rhs;
  }
  return lhs;
}

template <typename vecT, typename valT, typename std::enable_if<is_container<vecT>::value && is_number<valT>::value, int>::type = 0>
vecT operator*(vecT lhs, const valT & rhs){
  return lhs *= rhs;
}

template <typename vecT, typename valT, typename std::enable_if<is_container<vecT>::value && is_number<valT>::value, int>::type = 0>
vecT operator*(const valT & lhs, vecT rhs){
  return rhs *= lhs;
}


template <typename vecT, typename valT, typename std::enable_if<is_container<vecT>::value && is_number<valT>::value, int>::type = 0>
vecT operator/=(vecT & lhs, const valT & rhs){
  return lhs *= (1/rhs);
}

template <typename vecT, typename valT, typename std::enable_if<is_container<vecT>::value && is_number<valT>::value, int>::type = 0>
vecT operator/(vecT lhs, const valT & rhs){
  return lhs *= (1/rhs);
}



//square


template <typename vecT, typename std::enable_if<is_container<vecT>::value, int>::type = 0>
double vsq(const vecT & vec){
  double out = 0;
  for (auto x : vec) {
	out += x*x;
  }
  return out;
}


template <typename T1, std::size_t N = 0> struct is_container_of_num : std::false_type {};
template<> struct is_container_of_num<std::vector< int >> : std::true_type { };
template<> struct is_container_of_num<std::vector< double >> : std::true_type { };
template <std::size_t N> struct is_container_of_num<std::array<int, N> > : std::true_type { };
template <std::size_t N> struct is_container_of_num<std::array<double, N> > : std::true_type { };

//print
template <typename vecT, typename std::enable_if<is_container_of_num<vecT>::value, int>::type = 0>
std::ostream& operator<<(std::ostream& os, const vecT & vec) {
  os<< "(";
  auto it= vec.begin();
  for (auto x : vec) {
	if (it != vec.end()-1){ 
	  os << x << ", ";
	  it++;
	} else {
	  os << x;
	}
  }
  os << ")";
  return os;
}

template <typename T1, std::size_t N = 0> struct is_container_in_container : std::false_type {};
template <typename T1> struct is_container_in_container<std::vector<std::vector<T1>> > : std::true_type { };

template <typename vecvecT, typename std::enable_if<is_container_in_container<vecvecT>::value, int>::type = 0>
std::ostream& operator<<(std::ostream& os, const vecvecT & vec) {
  for (std::size_t i =0; i< vec[0].size(); i++){ 
	auto it= vec.begin();
	for (auto x : vec) {
	  if (it != vec.end()-1){ 
		os << x[i] << '\t';
		it++;
	  } else {
	  os << x[i] << '\n';
	  }
	}
  }
  return os;
}



#endif

