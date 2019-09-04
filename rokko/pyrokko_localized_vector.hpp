/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef PYROKKO_LOCALIZED_VECTOR_HPP
#define PYROKKO_LOCALIZED_VECTOR_HPP

#include <rokko/localized_vector.hpp>

namespace rokko {

class wrap_localized_vector {
public:

  wrap_localized_vector() {
    _ptr = new localized_vector<double>();
  }

  wrap_localized_vector(int size) {
    _ptr = new localized_vector<double>(size);
  }

  localized_vector<double>& obj() {
    return *_ptr;
  }

  localized_vector<double> const & obj() const {
    return *_ptr;
  }

  void print() const {
    obj().print();
  }

private:
  localized_vector<double>* _ptr;
};

} // end namespace rokko

#endif // PYROKKO_LOCALIZED_VECTOR_HPP
