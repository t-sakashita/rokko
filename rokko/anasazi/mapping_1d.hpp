/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_ANASAZI_MAPPING_1D_HPP
#define ROKKO_ANASAZI_MAPPING_1D_HPP

#include <rokko/mpi_communicator.hpp>
#include <rokko/mapping_1d.hpp>

#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>

namespace rokko {

namespace anasazi {

class mapping_1d : public detail::ps_mapping_1d_base {
public:
  explicit mapping_1d() : mapping_1d(0) {}
  explicit mapping_1d(int dim) : mapping_1d(dim, mpi_comm{MPI_COMM_WORLD}) {}
  explicit mapping_1d(int dim, mpi_comm const& mpi_comm_in)
    : detail::ps_mapping_1d_base(dim, mpi_comm_in), ep_comm_(new Epetra_MpiComm(mpi_comm_in.get_comm())), map_(new Epetra_Map(get_dim(), 0, *ep_comm_)) {
  }
  void init(int dim, mpi_comm const& mpi_comm_in) {
    set_dim(dim);
    set_mpi_comm(mpi_comm_in);
    ep_comm_ = new Epetra_MpiComm(mpi_comm_in.get_comm());
    map_ = new Epetra_Map(dim, 0, *ep_comm_);
  }
  int get_num_local_rows() const { return map_->NumMyElements();; }
  int start_row() const {
    return get_epetra_map().MinMyGID();
  }
  int end_row() const {
    return get_epetra_map().MaxMyGID() + 1; // to follow C++ convention
  }

  Epetra_MpiComm const& get_epetra_comm() const { return *ep_comm_; }
  Epetra_Map const& get_epetra_map() const { return *map_; }
  const ps_mapping_1d_base* get_impl() const { return this; }

private:
  const Epetra_MpiComm* ep_comm_;
  const Epetra_Map* map_;
};

} // namespace anasazi

} // namespace rokko

#endif // ROKKO_ANASAZI_MAPPING_1D_HPP
