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

#ifndef ROKKO_LOCALIZED_ROWVECTOR_HPP
#define ROKKO_LOCALIZED_ROWVECTOR_HPP

#include <Eigen/Dense>

namespace rokko {

template<typename T, int COLS = Eigen::Dynamic>
class localized_rowvector : public Eigen::Matrix<T, 1, COLS> {
public:
  typedef T value_type;
  typedef Eigen::Matrix<value_type, 1, COLS> super_type;
  typedef localized_rowvector<value_type> vector_type;

  localized_rowvector() : super_type() {}
  localized_rowvector(int size) : super_type(size) {}

  template <typename U>
  localized_rowvector(U const & other) : super_type(other) {}
  template <typename U>
  localized_rowvector& operator=(T const& other) { super_type::operator=(other); return *this; }
};

typedef localized_rowvector<float> flvector;
typedef localized_rowvector<double> dlvector;
typedef localized_rowvector<std::complex<float> > clvector;
typedef localized_rowvector<std::complex<double> > zlvector;

} // namespace rokko

#endif // ROKKO_LOCALIZED_ROWVECTOR_HPP
