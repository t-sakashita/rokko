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


namespace rokko {

class wrap_localized_vector {
public:

  wrap_localized_vector() {
    _ptr = new Eigen::VectorXd();
  }

  wrap_localized_vector(int size) {
    _ptr = new Eigen::VectorXd(size);
  }

  Eigen::VectorXd& obj() {
    return *_ptr;
  }

  Eigen::VectorXd const & obj() const {
    return *_ptr;
  }

  void print() const {
    std::cout << obj() << std::endl;
  }

private:
  Eigen::VectorXd* _ptr;
};

} // end namespace rokko

#endif // PYROKKO_LOCALIZED_VECTOR_HPP
