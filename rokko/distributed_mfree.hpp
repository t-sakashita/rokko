/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_DISTRIBUTED_MFREE_HPP
#define ROKKO_DISTRIBUTED_MFREE_HPP

namespace rokko {

class distributed_mfree {
public:
  distributed_mfree() {}
  ~distributed_mfree() {}

  virtual void multiply(const double* x, double* y) const = 0;
  virtual int get_dim() const = 0;
  virtual int get_local_offset() const = 0;
  virtual int get_num_local_rows() const = 0;
};

class distributed_mfree_slepc : public rokko::distributed_mfree {
public:
  virtual void diagonal(double* x) const = 0;
};

} // end namespace rokko

#endif // ROKKO_DISTRIBUTED_MFREE_HPP
