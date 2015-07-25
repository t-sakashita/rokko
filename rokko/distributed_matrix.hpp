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

#ifndef ROKKO_DISTRIBUTED_MATRIX_HPP
#define ROKKO_DISTRIBUTED_MATRIX_HPP

#include <rokko/grid.hpp>
#include <rokko/matrix_major.hpp>
#include <rokko/mapping_bc.hpp>
#include <rokko/blacs/blacs_wrap.h>
#include <rokko/pblas.h>

#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <boost/type_traits/is_same.hpp>
#include <boost/throw_exception.hpp>

namespace rokko {

class parallel_dense_solver;

template<typename T, typename MATRIX_MAJOR = rokko::matrix_col_major>
class distributed_matrix {
public:
  typedef T value_type;
  distributed_matrix(mapping_bc<MATRIX_MAJOR> const& map_in) {
    bool is_col_major = boost::is_same<MATRIX_MAJOR, matrix_col_major>::value;
    if (is_col_major != map.is_col_major()) {
      std::cerr << "error: matrix major of tamplate parameter and one of given mapping are different." << std::endl;
      throw;
    }
    map = map_in;
    allocate_array();
  }
  template<typename SOLVER>
  distributed_matrix(int m_global_in, const grid& g_in, SOLVER const& solver_in)// need fix
    : map(solver_in.optimized_mapping(m_global_in, g_in)) {
    allocate_array();
  }
  template<typename SOLVER>
  distributed_matrix(int m_global_in, int n_global_in, const grid& g_in, SOLVER const& solver_in) {
    bool is_col_major = boost::is_same<MATRIX_MAJOR, matrix_col_major>::value;
    if (is_col_major != map.is_col_major()) {
      std::cerr << "error: matrix major of tamplate parameter and one of given mapping are different." << std::endl;
      throw;
    }
    map = solver_in.optimized_mapping(m_global_in, g_in);
    allocate_array();
  }
  ~distributed_matrix() {
    delete[] array;
    array = 0;
  }

  void allocate_array() {
    array = new value_type[map.get_length_array()];
    if (array == 0) {
      std::cerr << "failed to allocate array." << std::endl;
      MPI_Abort(map.get_comm(), 3);
    }
  }

  value_type* get_array_pointer() { return array; }
  const value_type* get_array_pointer() const { return array; }

  void set_mapping(mapping_bc<MATRIX_MAJOR> const& map_in) { map = map_in; }

  mapping_bc<MATRIX_MAJOR>& get_mapping() { return map; }

  void set_local(int local_i, int local_j, value_type value) {
    array[map.get_array_index(local_i, local_j)] = value;
  }

  void update_local(int local_i, int local_j, value_type value) {
    array[map.get_array_index(local_i, local_j)] += value;
  }

  value_type get_local(int local_i, int local_j) const {
    return array[map.get_array_index(local_i, local_j)];
  }

  void set_global(int global_i, int global_j, value_type value) {
    if ((map.is_gindex(global_i, global_j)))
      set_local(map.translate_g2l_row(global_i), map.translate_g2l_col(global_j), value);
  }

  void update_global(int global_i, int global_j, value_type value) {
    if ((map.is_gindex(global_i, global_j)))
      update_local(map.translate_g2l_row(global_i), map.translate_g2l_col(global_j), value);
  }

  value_type get_global(int global_i, int global_j) const {
    return get_local(map.translate_g2l_row(global_i), map.translate_g2l_col(global_j));
  }

  value_type get_global_checked(int global_i, int global_j) const {
    if ((map.is_gindex(global_i, global_j))) {
      return get_local(map.translate_g2l_row(global_i), map.translate_g2l_col(global_j));
    } else {
      throw std::out_of_range("element not on this process.");
    }
  }

  void set_zeros() {
    for (int local_i=0; local_i<map.get_m_local(); ++local_i) {
      for (int local_j=0; local_j<map.get_n_local(); ++local_j) {
        set_local(local_i, local_j, 0);
      }
    }
  }

  template<class FUNC>
  void generate(FUNC func) {
    for(int local_i = 0; local_i < map.get_m_local(); ++local_i) {
      for(int local_j = 0; local_j < map.get_n_local(); ++local_j) {
        int global_i = map.translate_l2g_row(local_i);
        int global_j = map.translate_l2g_col(local_j);
        set_local(local_i, local_j, func(global_i, global_j));
      }
    }
  }

  void print(std::ostream& os = std::cout) const;

  // map member function
  int get_mb() const { return map.get_mb(); }
  int get_nb() const { return map.get_nb(); }

  int get_m_global() const { return map.get_m_global(); }
  int get_n_global() const { return map.get_n_global(); }
  int get_m_local() const { return map.get_m_local(); }
  int get_n_local() const { return map.get_n_local(); }
  int translate_l2g_row(const int& local_i) const { return map.translate_l2g_row(local_i); }
  int translate_l2g_col(const int& local_j) const { return map.translate_l2g_col(local_j); }
  int translate_g2l_row(const int& global_i) const { return map.translate_g2l_row(global_i); }
  int translate_g2l_col(const int& global_j) const { return map.translate_g2l_col(global_j); }
  bool is_gindex_myrow(const int& global_i) const { return map.is_gindex_myrow(global_i); }
  bool is_gindex_mycol(const int& global_j) const { return map.is_gindex_mycol(global_j); }
  bool is_gindex(const int& global_i, const int& global_j) const { return map.is_gindex(global_i, global_j); }

  int get_length_array() const { return map.get_length_array(); }
  int get_lld() const { return map.get_lld(); };
  int get_default_length_array() const { return map.get_default_length_array(); }
  int get_array_index(int local_i, int local_j) const { return map.get_array_index(local_i, local_j); }
  bool is_row_major() const { return map.is_row_major(); }
  bool is_col_major() const { return map.is_col_major(); }
  grid const& get_grid() const { return map.get_grid(); }
  int get_nprow() const { return map.get_nprow(); }
  int get_npcol() const { return map.get_npcol(); }
  int get_nprocs() const { return map.get_nprocs(); }
  int get_myrank() const { return map.get_myrank(); }
  int get_myrow() const { return map.get_myrow(); }
  int get_mycol() const { return map.get_mycol(); }

  // temporary
  void set_default_lld() { map.set_default_lld(); }
  void set_default_length_array() { map.set_default_length_array(); }

private:
  mapping_bc<MATRIX_MAJOR> map;
  value_type* array;
};

template<typename T, typename MATRIX_MAJOR>
void distributed_matrix<T, MATRIX_MAJOR>::print(std::ostream& os) const {
  for (int proc = 0; proc < map.get_nprocs(); ++proc) {
    if (proc == map.get_myrank()) {
      os << "Rank = " << map.get_myrank() << ", myrow = " << map.get_myrow() << ", mycol = " << map.get_mycol() << std::endl;
      for (int local_i = 0; local_i < map.get_m_local(); ++local_i) {
        for (int local_j = 0; local_j < map.get_n_local(); ++local_j)
          os << "  " << get_local(local_i, local_j);
        os << std::endl;
      }
      os.flush();
    }
    MPI_Barrier(map.get_comm());
  }
}

template<typename T, typename MATRIX_MAJOR>
void print_matrix(const rokko::distributed_matrix<T, MATRIX_MAJOR>& mat) { mat.print(); }

template<typename T, typename MATRIX_MAJOR>
std::ostream& operator<<(std::ostream& os, rokko::distributed_matrix<T, MATRIX_MAJOR> const& mat) {
  mat.print(os);
  return os;
}

// C = alpha A * B + beta C
template<typename T, typename MATRIX_MAJOR>
void product(typename distributed_matrix<T, MATRIX_MAJOR>::value_type alpha,
             const distributed_matrix<T, MATRIX_MAJOR>& matA, bool transA,
             const distributed_matrix<T, MATRIX_MAJOR>& matB, bool transB,
             typename distributed_matrix<T, MATRIX_MAJOR>::value_type beta,
             distributed_matrix<T, MATRIX_MAJOR>& matC) {
  int ictxt = ROKKO_blacs_get(-1, 0);
  char char_grid_major = (matA.get_grid().is_row_major() ? 'R' : 'C');
  ROKKO_blacs_gridinit(&ictxt, char_grid_major, matA.get_nprow(), matA.get_npcol());

  int descA[9], descB[9], descC[9];
  int info = ROKKO_descinit(descA, matA.get_m_global(), matA.get_n_global(), matA.get_mb(),
                            matA.get_nb(), 0, 0, ictxt, matA.get_lld());
  info = ROKKO_descinit(descB, matB.get_m_global(), matB.get_n_global(), matB.get_mb(),
                        matB.get_nb(), 0, 0, ictxt, matB.get_lld());
  info = ROKKO_descinit(descC, matC.get_m_global(), matC.get_n_global(), matC.get_mb(),
                        matC.get_nb(), 0, 0, ictxt, matC.get_lld());

  char char_transA = (transA ? 'T' : 'N');
  char char_transB = (transB ? 'T' : 'N');
  PBLASE_pgemm(char_transA, char_transB, matA.get_m_global(), matB.get_n_global(),
               matA.get_n_global(), alpha, matA.get_array_pointer(), 1, 1, descA,
               matB.get_array_pointer(), 1, 1, descB, beta,
               matC.get_array_pointer(), 1, 1, descC);
  ROKKO_blacs_gridexit(&ictxt);
}

// Y = alpha A * X + beta Y
template<typename T, typename MATRIX_MAJOR>
void product_v(typename distributed_matrix<T, MATRIX_MAJOR>::value_type alpha,
               const distributed_matrix<T, MATRIX_MAJOR>& matA, bool transA,
               const distributed_matrix<T, MATRIX_MAJOR>& vecX, bool transX, int xindex,
               typename distributed_matrix<T, MATRIX_MAJOR>::value_type beta,
               distributed_matrix<T, MATRIX_MAJOR>& vecY, bool transY, int yindex) {
  int ictxt = ROKKO_blacs_get(-1, 0);
  char char_grid_major = (matA.get_grid().is_row_major() ? 'R' : 'C');
  ROKKO_blacs_gridinit(&ictxt, char_grid_major, matA.get_nprow(), matA.get_npcol());

  int descA[9], descX[9], descY[9];
  int info = ROKKO_descinit(descA, matA.get_m_global(), matA.get_n_global(), matA.get_mb(),
                            matA.get_nb(), 0, 0, ictxt, matA.get_lld());
  info = ROKKO_descinit(descX, vecX.get_m_global(), vecX.get_n_global(), vecX.get_mb(),
                        vecX.get_nb(), 0, 0, ictxt, vecX.get_lld());
  info = ROKKO_descinit(descY, vecY.get_m_global(), vecY.get_n_global(), vecY.get_mb(),
                        vecY.get_nb(), 0, 0, ictxt, vecY.get_lld());

  char char_transA = (transA ? 'T' : 'N');
  int ix = (transX ? xindex + 1 : 1);
  int jx = (transX ? 1 : xindex + 1);
  int incx = (transX ? vecX.get_m_global() : 1);
  int iy = (transY ? yindex + 1: 1);
  int jy = (transY ? 1 : yindex + 1);
  int incy = (transY ? vecY.get_m_global() : 1);
  PBLASE_pgemv(char_transA, matA.get_m_global(), matA.get_n_global(), alpha,
               matA.get_array_pointer(), 1, 1, descA,
               vecX.get_array_pointer(), ix, jx, descX, incx, beta,
               vecY.get_array_pointer(), iy, jy, descY, incy);
  ROKKO_blacs_gridexit(&ictxt);
}

// dot = X * Y
template<typename T, typename MATRIX_MAJOR> 
T dot_product(const distributed_matrix<T, MATRIX_MAJOR>& vecX, bool transX, int xindex,
              const distributed_matrix<T, MATRIX_MAJOR>& vecY, bool transY, int yindex) {
  int ictxt = ROKKO_blacs_get(-1, 0);
  char char_grid_major = (vecX.get_grid().is_row_major() ? 'R' : 'C');
  ROKKO_blacs_gridinit(&ictxt, char_grid_major, vecX.get_nprow(), vecX.get_npcol());

  int descX[9], descY[9];
  int info = ROKKO_descinit(descX, vecX.get_m_global(), vecX.get_n_global(), vecX.get_mb(),
                            vecX.get_nb(), 0, 0, ictxt, vecX.get_lld());
  info = ROKKO_descinit(descY, vecY.get_m_global(), vecY.get_n_global(), vecY.get_mb(),
                        vecY.get_nb(), 0, 0, ictxt, vecY.get_lld());

  int n = (transX ? vecX.get_n_global() : vecX.get_m_global());
  int ix = (transX ? xindex + 1 : 1);
  int jx = (transX ? 1 : xindex + 1);
  int incx = (transX ? vecX.get_m_global() : 1);
  int iy = (transY ? yindex + 1: 1);
  int jy = (transY ? 1 : yindex + 1);
  int incy = (transY ? vecY.get_m_global() : 1);
  T dot = PBLASE_pdotc(n, vecX.get_array_pointer(), ix, jx, descX, incx,
                       vecY.get_array_pointer(), iy, jy, descY, incy);
  ROKKO_blacs_gridexit(&ictxt);
  return dot;
}

} // namespace rokko

#endif // ROKKO_DISTRIBUTED_MATRIX_HPP
