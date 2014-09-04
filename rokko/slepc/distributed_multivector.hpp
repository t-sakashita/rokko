/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SLEPC_DISTRIBUTED_MULTIVECTOR_H
#define ROKKO_SLEPC_DISTRIBUTED_MULTIVECTOR_H



namespace rokko {

class distributed_multivector_slepc {
public:
  distributed_multivector_slepc() {}
  distributed_multivector_slepc(int blocksize) : map_(map) {
  }
  ~distributed_multivector_slepc() {}
  void init_random() {
  }
private:
  mapping_1d map_;
};

} // namespace rokko

#endif // ROKKO_SLEPC_DISTRIBUTED_MULTIVECTOR_H
