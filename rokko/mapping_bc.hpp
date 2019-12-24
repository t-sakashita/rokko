/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MAPPING_BC_HPP
#define ROKKO_MAPPING_BC_HPP

#include <array>
#include <rokko/grid.hpp>
#include <rokko/matrix_major.hpp>
#include <rokko/mapping_global2local.hpp>
#include <rokko/mapping_local2array.hpp>
#include <rokko/scalapack.hpp>

namespace rokko {

template <typename MATRIX_MAJOR>
class mapping_bc : public mapping_global2local, public mapping_local2array<MATRIX_MAJOR> {
public:
  using matrix_major = MATRIX_MAJOR;
  explicit mapping_bc() {}
  explicit mapping_bc(int global_dim, int block_size)
    : mapping_global2local(global_dim, block_size, grid()),
      mapping_local2array<MATRIX_MAJOR>() {
    set_blacs_descriptor();
  }
  // default_mapping func of solver interface for elpa, elemental, scalapack
  explicit mapping_bc(int global_dim, int block_size, MATRIX_MAJOR matrix_major, grid const& g_in)
    : mapping_global2local(global_dim, block_size, g_in),
      mapping_local2array<MATRIX_MAJOR>() {
    set_blacs_descriptor();
  }
  
  // default_mapping func of solver interface for eigen_exa
  explicit mapping_bc(int global_dim, int block_size, std::array<int,2> const& padded_size,
                      grid const& g_in)
    : mapping_global2local(global_dim, block_size, g_in),
      mapping_local2array<MATRIX_MAJOR>(padded_size) {
    set_blacs_descriptor();
  }
  
  explicit mapping_bc(int global_dim, int block_size, int lld, grid const& g_in)
    : mapping_global2local(global_dim, block_size, g_in),
      mapping_local2array<MATRIX_MAJOR>(lld) {
    set_blacs_descriptor();
  }
  
  explicit mapping_bc(int global_dim, int block_size, grid const& g_in)
    : mapping_global2local(global_dim, block_size, g_in),
      mapping_local2array<MATRIX_MAJOR>() {
    set_blacs_descriptor();
  }
  
  explicit mapping_bc(std::array<int,2> const& global_size_in, std::array<int,2> const& block_size_in, grid const& g_in)
    : mapping_global2local(global_size_in, block_size_in, g_in),
      mapping_local2array<MATRIX_MAJOR>() {
    set_blacs_descriptor();
  }

  void set_blacs_descriptor() {
    int m = mapping_global2local::get_m_global();
    int n = mapping_global2local::get_n_global();
    int mb = mapping_global2local::get_mb();
    int nb = mapping_global2local::get_nb();
    int blacs_context = mapping_global2local::get_grid().get_blacs_context();
    int lld = mapping_local2array<matrix_major>::get_lld();
    if (lld == 0)  lld = 1;  // to avoid segmentation fault in descinit
    int info = scalapack::descinit(blacs_descriptor, m, n, mb, nb, 0, 0, blacs_context, lld);
    if (info) {
      std::cerr << "error info=" << info << " at descinit function of descA " << "mA="
              << m << "  nA=" << n << "  lld=" << lld << "." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, info);
    }
  }
  const std::array<int,9>& get_blacs_descriptor() const { return blacs_descriptor; }

private:
  std::array<int,9> blacs_descriptor;
};

} // namespace rokko

#endif // ROKKO_MAPPING_BC_HPP
