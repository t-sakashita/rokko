#ifndef ROKKO_SCALAPACK_DISTRIBUTED_MATRIX_HPP
#define ROKKO_SCALAPACK_DISTRIBUTED_MATRIX_HPP


#include <cstdlib>
//#include <mpi.h>
#include <rokko/grid.hpp>

#include <rokko/distributed_matrix.hpp>

namespace rokko {

template<>
class distributed_matrix<rokko::scalapack>
{
public:
  distributed_matrix(int m_global, int n_global, const grid<rokko::scalapack>& g)
    : m_global(m_global), n_global(n_global), g(g), myrank(g.myrank), nprocs(g.nprocs), myrow(g.myrow), mycol(g.mycol), nprow(g.nprow), npcol(g.npcol), ictxt(g.ictxt)
  {
    // ローカル行列の形状を指定
    mb = m_global / g.nprow;
    if (mb == 0) mb = 1;
    //mb = 10;
    nb = n_global / g.npcol;
    if (nb == 0) nb = 1;
    //nb = 10;
    // mbとnbを最小値にそろえる．（注意：pdsyevではmb=nbでなければならない．）
    mb = min(mb, nb);
    nb = mb;

    const int ZERO=0, ONE=1;
    m_local = numroc_( m_global, mb, myrow, ZERO, nprow );
    n_local = numroc_( n_global, nb, mycol, ZERO, npcol );

    lld = m_local;

    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == g.myrank) {
	std::cout << "proc=" << proc << std::endl;
	std::cout << "  mb=" << mb << "  nb=" << nb << std::endl;
	std::cout << "  mA=" << m_local << "  nprow=" << g.nprow << std::endl;
	std::cout << "  nA=" << n_local << "  npcol=" << g.npcol << std::endl;
	std::cout << " m_local=" << m_local << " n_local=" << n_local << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    array = new double[m_local * n_local];
    if (array == NULL) {
      cerr << "failed to allocate array." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 3);
    }
  }

  distributed_matrix(int m_global, int n_global, const grid<rokko::scalapack>& g, bool mb)
    : m_global(m_global), n_global(n_global), g(g), myrank(g.myrank), nprocs(g.nprocs), myrow(g.myrow), mycol(g.mycol), nprow(g.nprow), npcol(g.npcol), ictxt(g.ictxt)
  {
    // ローカル行列の形状を指定
    /*
    mb = m_global / g.nprow;
    if (mb == 0) mb = 1;
    //mb = 10;
    nb = n_global / g.npcol;
    if (nb == 0) nb = 1;
    //nb = 10;
    // mbとnbを最小値にそろえる．（注意：pdsyevではmb=nbでなければならない．）
    mb = min(mb, nb);
    */
    mb = 1;
    nb = mb;

    const int ZERO=0, ONE=1;
    m_local = numroc_( m_global, mb, myrow, ZERO, nprow );
    n_local = numroc_( n_global, nb, mycol, ZERO, npcol );

    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == g.myrank) {
	std::cout << "proc=" << proc << std::endl;
	std::cout << "  mb=" << mb << "  nb=" << nb << std::endl;
	std::cout << "  mA=" << m_local << "  nprow=" << g.nprow << std::endl;
	std::cout << "  nA=" << n_local << "  npcol=" << g.npcol << std::endl;
	std::cout << " m_local=" << m_local << " n_local=" << n_local << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    array = new double[m_local * n_local];
    if (array == NULL) {
      cerr << "failed to allocate array." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 3);
    }
  }

  ~distributed_matrix()
  {
    //std::cout << "Destructor ~Distributed_Matrix()" << std::endl;
    delete[] array;
    array = NULL;
  }

  double* get_array()
  {
    return array;
  }

  int translate_l2g_row(const int& local_i) const
  {
    return (g.myrow * mb) + local_i + (local_i / mb) * mb * (g.nprow - 1);
  }

  int translate_l2g_col(const int& local_j) const
  {
    return (g.mycol * nb) + local_j + (local_j / nb) * nb * (g.npcol - 1);
  }

  int translate_g2l_row(const int& global_i) const
  {
    const int local_offset_block = global_i / mb;
    return (local_offset_block - g.myrow) / g.nprow * mb + global_i % mb;
  }

  int translate_g2l_col(const int& global_j) const
  {
    const int local_offset_block = global_j / nb;
    return (local_offset_block - g.mycol) / g.npcol * nb + global_j % nb;
  }

  bool is_gindex_myrow(const int& global_i) const
  {
    int local_offset_block = global_i / mb;
    return (local_offset_block % g.nprow) == g.myrow;
  }

  bool is_gindex_mycol(const int& global_j) const
  {
    int local_offset_block = global_j / nb;
    return (local_offset_block % g.npcol) == g.mycol;
  }

  void set_local(int local_i, int local_j, double value)
  {
    array[local_j * m_local + local_i] = value;
  }

  void set_global(int global_i, int global_j, double value)
  {
    if ((is_gindex_myrow(global_i)) && (is_gindex_mycol(global_j)))
      set_local(translate_g2l_row(global_i), translate_g2l_col(global_j), value);
  }

  void calculate_grid(int proc_rank, int& proc_row, int& proc_col) const
  {
    proc_row = proc_rank / npcol;
    proc_col = proc_rank % npcol;
  }

  int calculate_grid_row(int proc_rank) const
  {
    return proc_rank / npcol;
  }

  int calculate_grid_col(int proc_rank) const
  {
    return proc_rank % npcol;
  }

  void print() const
  {
    /* each proc prints it's local_array out, in order */
    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == myrank) {
	printf("Rank = %d  myrow=%d mycol=%d\n", myrank, myrow, mycol);
	printf("Local Matrix:\n");
	for (int ii=0; ii<m_local; ++ii) {
	  for (int jj=0; jj<n_local; ++jj) {
	    printf("%3.2f ",array[jj * m_local + ii]);
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
  int ictxt, nprow, npcol;
  const grid<rokko::scalapack>& g;
  int lld;

private:
  int info;
};

template<>
void print_matrix(const rokko::distributed_matrix<rokko::scalapack>& mat)
{
  // each proc prints it's local_array out, in order
  for (int proc=0; proc<mat.nprocs; ++proc) {
    if (proc == mat.myrank) {
      printf("Rank = %d  myrow=%d mycol=%d\n", mat.myrank, mat.myrow, mat.mycol);
      printf("Local Matrix:\n");
      for (int ii=0; ii<mat.m_local; ++ii) {
	for (int jj=0; jj<mat.n_local; ++jj) {
	  printf("%3.2f ", mat.array[jj * mat.m_local + ii]);
	}
	printf("\n");
      }
      printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

} // namespace rokko

#endif // ROKKO_SCALAPACK_DISTRIBUTED_MATRIX_HPP

