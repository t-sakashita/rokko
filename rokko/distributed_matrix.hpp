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

template<typename T, typename MATRIX_MAJOR = rokko::matrix_row_major>
class distributed_matrix {
public:
  typedef T value_type;
  
  template<typename SOLVER>
  distributed_matrix(int m_global_in, int n_global_in, const grid& g_in, SOLVER const& solver_in) :
    g(g_in) {
    initialize(m_global_in, n_global_in, g_in, solver_in);
  }
  template<typename SOLVER>
  distributed_matrix(mapping_bc const& map, SOLVER const& solver_in) {
    initialize(map, solver_in);
  }
  ~distributed_matrix() {
    delete[] array;
    array = 0;
  }

  void initialize(mapping_bc const& map) {
    set_mapping(map);
    array = new value_type[length_array];
    if (array == 0) {
      std::cerr << "failed to allocate array." << std::endl;
      MPI_Abort(g.get_comm(), 3);
    }
  }
  template<typename SOLVER>
  void initialize(int m_global, int n_global, const grid& g_in, SOLVER const& solver_in) {
    initialize(m_global, g_in, solver_in);
  }
  template<typename SOLVER>
  void initialize(int dim, const grid& g_in, SOLVER const& solver_in) {
    initialize(solver_in.optimized_mapping(g_in, dim), solver_in);
  }
  template<typename SOLVER>
  void initialize(mapping_bc const& map, SOLVER const& solver_in) {
    set_mapping(map);
    set_local_size(calculate_row_size(), calculate_col_size());
    solver_in.optimized_matrix_size(*this);
    array = new value_type[length_array];
    if (array == 0) {
      std::cerr << "failed to allocate array." << std::endl;
      MPI_Abort(g.get_comm(), 3);
    }
  }

  value_type* get_array_pointer() { return array; }
  const value_type* get_array_pointer() const { return array; }

  int get_length_array() const { return length_array; }
  const mapping_bc& get_mapping() { return map_; }
  void set_mapping(mapping_bc const& map) {
    map_ = map;
    // substitute sizes from map_
    g = map_.get_grid();
    myrank = g.get_myrank(); nprocs = g.get_nprocs();
    myrow = g.get_myrow(); mycol = g.get_mycol();
    nprow = g.get_nprow(); npcol = g.get_npcol();
    m_global = map_.get_dim();  n_global = m_global;
    mb = map_.get_block_size();  nb = mb;
    // set strides
    stride_myrow = myrow * mb;
    stride_nprow = mb * (nprow - 1);
    stride_mycol = mycol * nb;
    stride_npcol = nb * (npcol - 1);
  }
  void set_length_array(int value) { length_array = value; }
  void set_block_size(int mb_in, int nb_in) {
    mb = mb_in;
    nb = nb_in;
  }

  const grid& get_grid() const { return g; }

  int get_mb() const { return mb; }
  int get_nb() const { return nb; }

  int get_nprow() const { return nprow; }
  int get_npcol() const { return npcol; }
  int get_nprocs() const { return nprocs; }
  int get_myrank() const { return myrank; }

  int get_m_global() const { return m_global; }
  int get_n_global() const { return n_global; }

  int get_m_local() const { return m_local; }
  int get_n_local() const { return n_local; }

  int get_myrow() const { return myrow; }
  int get_mycol() const { return mycol; }


  void set_local_size(int m_local_in, int n_local_in) {
    m_local = m_local_in;
    n_local = n_local_in;
  }

  void set_default_local_size() {
    set_local_size(calculate_row_size(), calculate_col_size());
  }

  int calculate_row_size(int proc_row) const {
    int tmp = m_global / mb;
    int local_num_block_rows = (tmp - proc_row -1) / nprow + 1;
    int rest_block_row = tmp % nprow; // size of a residue block (< mb)
    int local_rest_block_rows;
    if (proc_row == rest_block_row)
      local_rest_block_rows = m_global % mb;
    else
      local_rest_block_rows = 0;

    return  local_num_block_rows * mb + local_rest_block_rows;
  }

  int calculate_row_size() const {
    return calculate_row_size(myrow);
  }

  int calculate_col_size(int proc_col) const {
    int tmp = n_global / nb;
    int local_num_block_cols = (tmp - proc_col -1) / npcol + 1;
    int rest_block_col = tmp % npcol; // size of a residue block (< nb)
    int local_rest_block_cols;
    if (proc_col == rest_block_col) {
      local_rest_block_cols = n_global % nb;
    } else {
      local_rest_block_cols = 0;
    }
    return local_num_block_cols * nb + local_rest_block_cols;
  }

  int calculate_col_size() const {
    return calculate_col_size(mycol);
  }

  int get_lld() const { return lld; };

  void set_lld(int value) { lld = value; };
  void set_default_lld() { set_lld(get_default_lld()); }

  int get_default_lld() const {
    return boost::is_same<MATRIX_MAJOR, matrix_row_major>::value ? n_local : m_local;
  }
  int get_default_length_array() const {
    return boost::is_same<MATRIX_MAJOR, matrix_row_major>::value ?
      (m_local * lld) : (lld * n_local);
  }
  void set_default_length_array() { set_length_array(get_default_length_array()); }

  int get_array_index(int local_i, int local_j) const {
    return boost::is_same<MATRIX_MAJOR, matrix_row_major>::value ?
      (local_i * lld + local_j) : (local_i + local_j * lld);
  }
  
  int translate_l2g_row(const int& local_i) const {
    return stride_myrow + local_i + (local_i / mb) * stride_nprow;
  }

  int translate_l2g_col(const int& local_j) const {
    return stride_mycol + local_j + (local_j / nb) * stride_npcol;
  }

  int translate_g2l_row(const int& global_i) const {
    int local_offset_block = global_i / mb;
    return (local_offset_block - myrow) / nprow * mb + global_i % mb;
  }

  int translate_g2l_col(const int& global_j) const {
    const int local_offset_block = global_j / nb;
    return (local_offset_block - mycol) / npcol * nb + global_j % nb;
  }

  bool is_gindex_myrow(const int& global_i) const {
    int local_offset_block = global_i / mb;
    return (local_offset_block % nprow) == myrow;
  }

  bool is_gindex_mycol(const int& global_j) const {
    int local_offset_block = global_j / nb;
    return (local_offset_block % npcol) == mycol;
  }

  bool is_gindex(const int& global_i, const int& global_j) const {
    return is_gindex_myrow(global_i) && is_gindex_mycol(global_j);
  }

  void set_local(int local_i, int local_j, value_type value) {
    array[get_array_index(local_i, local_j)] = value;
  }

  void update_local(int local_i, int local_j, value_type value) {
    array[get_array_index(local_i, local_j)] += value;
  }

  value_type get_local(int local_i, int local_j) const {
    return array[get_array_index(local_i, local_j)];
  }

  void set_global(int global_i, int global_j, value_type value) {
    if ((is_gindex(global_i, global_j)))
      set_local(translate_g2l_row(global_i), translate_g2l_col(global_j), value);
  }

  void update_global(int global_i, int global_j, value_type value) {
    if ((is_gindex(global_i, global_j)))
      update_local(translate_g2l_row(global_i), translate_g2l_col(global_j), value);
  }

  value_type get_global(int global_i, int global_j) const {
    return get_local(translate_g2l_row(global_i), translate_g2l_col(global_j));
  }

  value_type get_global_checked(int global_i, int global_j) const {
    if ((is_gindex(global_i, global_j))) {
      return get_local(translate_g2l_row(global_i), translate_g2l_col(global_j));
    } else {
      throw std::out_of_range("element not on this process.");
    }
  }

  void set_zeros() {
    for (int local_i=0; local_i<m_local; ++local_i) {
      for (int local_j=0; local_j<n_local; ++local_j) {
        set_local(local_i, local_j, 0);
      }
    }
  }

  bool is_row_major() const {
    return boost::is_same<MATRIX_MAJOR, matrix_row_major>::value;
  }
  bool is_col_major() const {
    return boost::is_same<MATRIX_MAJOR, matrix_col_major>::value;
  }

  template<class FUNC>
  void generate(FUNC func) {
    for(int local_i = 0; local_i < m_local; ++local_i) {
      for(int local_j = 0; local_j < n_local; ++local_j) {
        int global_i = translate_l2g_row(local_i);
        int global_j = translate_l2g_col(local_j);
        set_local(local_i, local_j, func(global_i, global_j));
      }
    }
  }

  void print(std::ostream& os = std::cout) const;

private:
  mapping_bc map_;
  int m_global, n_global;
  value_type* array;
  int mb, nb;
  int m_local, n_local;
  int lld;
  int length_array;
  int stride_myrow, stride_nprow, stride_mycol, stride_npcol;
  grid g;
  // variables of class grid
  int myrank, nprocs;
  int myrow, mycol;
  int nprow, npcol;
};

template<typename T, typename MATRIX_MAJOR>
void distributed_matrix<T, MATRIX_MAJOR>::print(std::ostream& os) const {
  for (int proc = 0; proc < nprocs; ++proc) {
    if (proc == myrank) {
      os << "Rank = " << myrank << ", myrow = " << myrow << ", mycol = " << mycol << std::endl;
      for (int local_i = 0; local_i < m_local; ++local_i) {
        for (int local_j = 0; local_j < n_local; ++local_j)
          os << "  " << get_local(local_i, local_j);
        os << std::endl;
      }
      os.flush();
    }
    MPI_Barrier(g.get_comm());
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
  T dot;
  PBLASE_pdot(n, &dot, vecX.get_array_pointer(), ix, jx, descX, incx,
              vecY.get_array_pointer(), iy, jy, descY, incy);
  ROKKO_blacs_gridexit(&ictxt);
  return dot;
}

} // namespace rokko

#endif // ROKKO_DISTRIBUTED_MATRIX_HPP
