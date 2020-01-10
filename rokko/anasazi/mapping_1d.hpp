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

#ifndef ROKKO_ANASAZI_MAPPING_1D_HPP
#define ROKKO_ANASAZI_MAPPING_1D_HPP

#include <rokko/mpi_communicator.hpp>

#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Teuchos_RCPDecl.hpp>

namespace rokko {

namespace anasazi {

class mapping_1d {
public:
  explicit mapping_1d() : mapping_1d(0) {}
  explicit mapping_1d(int dim) : mapping_1d(dim, mpi_comm(MPI_COMM_WORLD)) {}
  explicit mapping_1d(int dim, mpi_comm const& mpi_comm_in) :
    dim_(dim), mpi_comm_(mpi_comm_in), ep_comm_(Epetra_MpiComm(mpi_comm_in.get_comm())), map_(dim_, 0, ep_comm_), num_local_rows_(map_.NumMyElements()) {
  }
  int get_dim() const { return dim_; }
  int get_num_local_rows() const { return num_local_rows_; }
  Epetra_MpiComm const& get_epetra_comm() const { return ep_comm_; }
  Epetra_Map const& get_epetra_map() const { return map_; }
  mpi_comm const& get_mpi_comm() const { return mpi_comm_; }

private:
  int dim_;
  mpi_comm mpi_comm_;
  Epetra_MpiComm ep_comm_;
  Epetra_Map map_;
  int num_local_rows_;
};

} // namespace anasazi

} // namespace rokko

#endif // ROKKO_ANASAZI_MAPPING_1D_HPP
