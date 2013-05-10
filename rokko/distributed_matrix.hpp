#ifndef ROKKO_DISTRIBUTED_H
#define ROKKO_DISTRIBUTED_H


#include <cstdlib>
//#include <mpi.h>
#include <rokko/grid.hpp>

namespace rokko {

class solver;

struct matrix_row_major {};

struct matrix_col_major {};

template<typename MATRIX_MAJOR = rokko::matrix_row_major>
class distributed_matrix {
public:
  template<typename GRID_MAJOR>
  distributed_matrix(int m_global, int n_global, const grid<GRID_MAJOR>& g_in)
    : m_global(m_global), n_global(n_global), myrank(g_in.myrank), nprocs(g_in.nprocs), myrow(g_in.myrow), mycol(g_in.mycol), nprow(g_in.nprow), npcol(g_in.npcol), g(g_in) {
    //g = new grid<GRID_MAJOR>(g_in);

    // ローカル行列の形状を指定
    mb = m_global / nprow;
    if (mb == 0) mb = 1;
    //mb = 10;
    nb = n_global / npcol;
    if (nb == 0) nb = 1;
    //nb = 10;
    // mbとnbを最小値にそろえる．（注意：pdsyevではmb=nbでなければならない．）
    //mb = min(mb, nb);
    mb = nb = 1;

    m_local = get_row_size();
    n_local = get_col_size();
    lld = get_lld();

#ifndef NDEBUG
    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == myrank) {
	std::cout << "proc=" << proc << std::endl;
	std::cout << "  mb=" << mb << "  nb=" << nb << std::endl;
	std::cout << "  mA=" << m_local << "  nprow=" << nprow << std::endl;
	std::cout << "  nA=" << n_local << "  npcol=" << npcol << std::endl;
	std::cout << " m_local=" << m_local << " n_local=" << n_local << std::endl;
        std::cout << " myrow=" << myrow << " mycol=" << mycol << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    array = new double[m_local * n_local];
    if (array == 0) {
      std::cerr << "failed to allocate array." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 3);
    }
  }

  template<typename GRID_MAJOR>
  distributed_matrix(int m_global, int n_global, const grid<GRID_MAJOR>& g_in, rokko::solver& solver_in)
    : m_global(m_global), n_global(n_global), myrank(g_in.myrank), nprocs(g_in.nprocs), myrow(g_in.myrow), mycol(g_in.mycol), nprow(g_in.nprow), npcol(g_in.npcol), g(g_in) {
    // determine mb, nb, lld, larray
    int larray;
    solver_in.optimized_matrix_size(m_global, nprow, npcol, mb, nb, lld, larray);

    // determine m_local, n_local from m_global, n_global, mb, nb
    m_local = get_row_size();
    n_local = get_col_size();

#ifndef NDEBUG
    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == myrank) {
	std::cout << "proc=" << proc << std::endl;
	std::cout << "  mb=" << mb << "  nb=" << nb << std::endl;
	std::cout << "  mA=" << m_local << "  nprow=" << nprow << std::endl;
	std::cout << "  nA=" << n_local << "  npcol=" << npcol << std::endl;
	std::cout << " m_local=" << m_local << " n_local=" << n_local << std::endl;
        std::cout << " myrow=" << myrow << " mycol=" << mycol << std::endl;
        std::cout << " lld=" << lld << std::endl;
        std::cout << " larray=" << larray << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    array = new double[larray];
    if (array == 0 {
      std::cerr << "failed to allocate array." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 3);
    }
  }

  ~distributed_matrix() {
    delete[] array;
    array = 0;
  }

  double* get_array() { return array; }

  int get_row_size() const {
    int tmp = m_global / mb;
    int local_num_block_rows = (tmp + nprow-1 - myrow) / nprow;
    int rest_block_row = ((m_global / mb) % nprow) % nprow; // 最後のブロックを持つプロセスの次のプロセス
    int local_rest_block_rows;
    std::cout << "local_num_block_rows=" << local_num_block_rows << std::endl;
    if (myrow == rest_block_row)
      local_rest_block_rows = m_global % mb;
    else
      local_rest_block_rows = 0;

    return  local_num_block_rows * mb + local_rest_block_rows;

    //const int local_offset_block = m_global / mb;
    //cout << "local_offset_block=" << local_offset_block;
    //return (local_offset_block - myrow) / nprow * mb + m_global % mb;
  }

  int get_col_size() const {
    int tmp = n_global / nb;
    int local_num_block_cols = (tmp + npcol-1 - mycol) / npcol;
    int rest_block_col = ((n_global / nb) % npcol) % npcol; // 最後のブロックを持つプロセスの次のプロセス
    int local_rest_block_cols;
    std::cout << "local_num_block_cols=" << local_num_block_cols << std::endl;
    if (myrow == rest_block_col)
      local_rest_block_cols = n_global % nb;
    else
      local_rest_block_cols = 0;

    return  local_num_block_cols * nb + local_rest_block_cols;

    //const int local_offset_block = n_global / nb;
    //return (local_offset_block - mycol) / npcol * nb + n_global % nb;
  }

  int get_lld() const;
  int get_array_index(int local_i, int local_j) const;

  int translate_l2g_row(const int& local_i) const {
    return (myrow * mb) + local_i + (local_i / mb) * mb * (nprow - 1);
  }

  int translate_l2g_col(const int& local_j) const {
    return (mycol * nb) + local_j + (local_j / nb) * nb * (npcol - 1);
  }

  int translate_g2l_row(const int& global_i) const {
    const int local_offset_block = global_i / mb;
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

  void set_local(int local_i, int local_j, double value) {
    array[get_array_index(local_i, local_j)] = value;
  }

  void set_global(int global_i, int global_j, double value) {
    if ((is_gindex_myrow(global_i)) && (is_gindex_mycol(global_j)))
      set_local(translate_g2l_row(global_i), translate_g2l_col(global_j), value);
  }

  void print() const {
    /* each proc prints it's local_array out, in order */
    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == myrank) {
	printf("Rank = %d  myrow=%d mycol=%d\n", myrank, myrow, mycol);
	printf("Local Matrix:\n");
	for (int local_i=0; local_i<m_local; ++local_i) {
	  for (int local_j=0; local_j<n_local; ++local_j) {
	    printf("%e ",array[get_array_index(local_i, local_j)]);
	  }
	  printf("\n");
	}
	printf("\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  int m_global, n_global;
  double* array;
  int mb, nb;
  int m_local, n_local;
  // variables of class Grid
  int myrank, nprocs;
  int myrow, mycol;
  int nprow, npcol;
  int lld;

  const grid_base& g;

private:
  int info;
};


template<>
inline int distributed_matrix<rokko::matrix_row_major>::get_lld() const {
  return m_local;
}

template<>
inline int distributed_matrix<rokko::matrix_col_major>::get_lld() const {
  return n_local;
}

template<>
inline int distributed_matrix<rokko::matrix_row_major>::get_array_index(int local_i, int local_j) const {
  return local_i * lld + local_j;
}


template<>
inline int distributed_matrix<rokko::matrix_col_major>::get_array_index(int local_i, int local_j) const {
  return  local_i + local_j * lld;
}

template<typename MATRIX_MAJOR>
void print_matrix(const rokko::distributed_matrix<MATRIX_MAJOR>& mat) {
  // each proc prints it's local_array out, in order
  for (int proc=0; proc<mat.nprocs; ++proc) {
    if (proc == mat.myrank) {
      printf("Rank = %d  myrow=%d mycol=%d\n", mat.myrank, mat.myrow, mat.mycol);
      printf("Local Matrix:\n");
      for (int local_i=0; local_i<mat.m_local; ++local_i) {
	for (int local_j=0; local_j<mat.n_local; ++local_j) {
	  printf("%e ", mat.array[mat.get_array_index(local_i, local_j)]);
	}
	printf("\n");
      }
      printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

} // namespace rokko

#endif // ROKKO_DISTRIBUTED_H
