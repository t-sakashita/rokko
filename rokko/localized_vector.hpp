/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_LOCALIZED_VECTOR_HPP
#define ROKKO_LOCALIZED_VECTOR_HPP

#include <Eigen/Dense>
#include <iostream>

namespace rokko {

template<typename T, int ROWS = Eigen::Dynamic>
class localized_vector : public Eigen::Matrix<T, ROWS, 1> {
public:
  typedef T value_type;
  typedef Eigen::Matrix<value_type, ROWS, 1> super_type;
  typedef localized_vector<value_type> vector_type;

  localized_vector() : super_type() {}
  localized_vector(int size) : super_type(size) {}

  template <typename U>
  localized_vector(U const & other) : super_type(other) {}
  template <typename U>
  localized_vector& operator=(T const& other) { super_type::operator=(other); return *this; }

  void print() const { std::cout << *this << std::endl; }
};

typedef localized_vector<float> flvector;
typedef localized_vector<double> dlvector;
typedef localized_vector<std::complex<float> > clvector;
typedef localized_vector<std::complex<double> > zlvector;

} // namespace rokko

#endif // ROKKO_LOCALIZED_VECTOR_HPP
