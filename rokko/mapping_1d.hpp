/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014-2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MAPPING_1D_HPP
#define ROKKO_MAPPING_1D_HPP

#include <rokko/grid_1d.hpp>

#include <vector>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>

namespace rokko {

class mapping_1d {
public:
  mapping_1d() : dim_(0), map_(dim_, 0, Epetra_MpiComm(MPI_COMM_WORLD)) {}
  explicit mapping_1d(int dim, grid_1d const& g) :
    dim_(dim), map_(dim_, 0, Epetra_MpiComm(g.get_comm())) {
    rows_.resize(map_.NumMyElements());
    map_.MyGlobalElements(&rows_[0]);
  }
  int dimension() const { return dim_; }
  int num_rows() const { return rows_.size(); }
  std::vector<int> const& rows() const { return rows_; }
  Epetra_Map const& get_epetra_map() const { return map_; }
private:
  int dim_;
  Epetra_Map map_;
  std::vector<int> rows_;
};

} // namespace rokko

#endif // ROKKO_MAPPING_1D_HPP
