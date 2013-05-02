#ifndef ROKKO_EIGEN_S_DISTRIBUTED_MATRIX_HPP
#define ROKKO_EIGEN_S_DISTRIBUTED_MATRIX_HPP

#include <mpi.h>
#include <cstdlib>

#include <rokko/scalapack/blacs.hpp>
#include <rokko/eigen_s/eigen_s.hpp>
#include <rokko/eigen_s/grid.hpp>
#include "../distributed_matrix.hpp"

namespace rokko {

template<>
class distributed_matrix<rokko::eigen_s>
{
public:
  // torus-wrap distribution used in EigenK and Elemental
  distributed_matrix(int m_global, int n_global, const grid<rokko::eigen_s>& g_in)
    : m_global(m_global), n_global(n_global), g(g_in), myrank(g_in.myrank), nprocs(g_in.nprocs), myrow(g_in.myrow), mycol(g_in.mycol), nprow(g_in.nprow), npcol(g_in.npcol)
  {
    int n = m_global;
    int nx = ((n-1)/nprow+1);
    int i1 = 6, i2 = 16*4, i3 = 16*4*2, nm;
    CSTAB_get_optdim(nx, i1, i2, i3, nm);  // return an optimized (possiblly) leading dimension of local block-cyclic matrix to nm.
    //int para_int = 0;   eigen_free_wrapper(para_int);

    int NB  = 64+32;
    int nmz = ((n-1)/nprow+1);
    nmz = ((nmz-1)/NB+1)*NB+1;
    int nmw = ((n-1)/npcol+1);
    nmw = ((nmw-1)/NB+1)*NB+1;
    cout << "nm=" << nm << endl;
    cout << "nmz=" << nmz << endl;
    cout << "nmw=" << nmw << endl;
    int larray = std::max(nmz,nm) * nmw;

    // calculate sizes of my proc's local part of distributed matrix
    mb = 1;
    nb = 1;
    const int ZERO=0, ONE=1;
    lld = nm;
    m_local = numroc_( m_global, mb, myrow, ZERO, nprow );  //(m_global + nprow - myrow) / nprow;
    n_local = numroc_( n_global, nb, mycol, ZERO, npcol );  //(n_global + npcol - mycol) / npcol;
    //m_local = nm;
    //n_local = (larray + (nm-1)) / nm;

    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == g.myrank) {
	cout << "proc=" << proc << endl;
	//cout << "  mb=" << mb << "  nb=" << nb << endl;
	//cout << "  mA=" << m_local << "  nprow=" << g.nprow << endl;
	//cout << "  nA=" << n_local << "  npcol=" << g.npcol << endl;
	cout << "nm(lld)=" << lld << endl;
	cout << " m_local=" << m_local << " n_local=" << n_local << endl;
     }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    array = new double[larray];  //m_local * n_local];
    if (array == NULL) {
      cerr << "failed to allocate array." << endl;
      MPI_Abort(MPI_COMM_WORLD, 3);
    }
    for (int ii=0; ii<larray; ++ii)
      array[ii] = -3;
  }

  ~distributed_matrix()
  {
    //cout << "Destructor ~Distributed_Matrix_EigenK()" << endl;
    delete[] array;
    array = NULL;
  }

  double* get_array()
  {
    return array;
  }

  int translate_l2g_row(const int& local_i) const
  {
    return local_i * g.nprow + g.myrow;
    //return (g.myrow * mb) + local_i + (local_i / mb) * mb * (g.nprow - 1);
  }

  int translate_l2g_col(const int& local_j) const
  {
    return local_j * g.npcol + g.mycol;
    //return (g.mycol * nb) + local_j + (local_j / nb) * nb * (g.npcol - 1);
  }

  int translate_g2l_row(const int& global_i) const
  {
    return (global_i - g.myrow) / g.nprow;
    //const int local_offset_block = global_i / mb;
    //return (local_offset_block - g.myrow) / g.nprow * mb + global_i % mb;
  }

  int translate_g2l_col(const int& global_j) const
  {
    return (global_j - g.mycol) / g.npcol;
    //const int local_offset_block = global_j / nb;
    //return (local_offset_block - g.mycol) / g.npcol * nb + global_j % nb;
  }

  bool is_gindex_myrow(const int& global_i) const
  {
    return (global_i % g.nprow) == g.myrow;
    //int local_offset_block = global_i / mb;
    //return (local_offset_block % g.nprow) == g.myrow;
  }

  bool is_gindex_mycol(const int& global_j) const
  {
    return (global_j % g.npcol) == g.mycol;
    //int local_offset_block = global_j / nb;
    //return (local_offset_block % g.npcol) == g.mycol;
  }

  void set_local(int local_i, int local_j, double value)
  {
    //array[local_j * lld + local_i] = value;
    array[local_j + local_i * lld] = value;
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
	    //printf("%3.2f ", array[jj * lld + ii]);
            //printf("%3.2f ", array[jj + ii * lld]);
            printf("%e ", array[jj + ii * lld]);
	  }
	  printf("\n");
	}
	printf("\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
}

  int m_global, n_global;
  int lld;
  double* array;
  int mb, nb;
  int m_local, n_local;
  // variables of class Grid
  int myrank, nprocs;
  int myrow, mycol;
  int nprow, npcol;
  const grid<rokko::eigen_s>& g;

private:
  int info;
};


template<>
void print_matrix(const rokko::distributed_matrix<rokko::eigen_s>& mat)
{
  // each proc prints it's local_array out, in order
  for (int proc=0; proc<mat.nprocs; ++proc) {
    if (proc == mat.myrank) {
      printf("Rank = %d  myrow=%d mycol=%d\n", mat.myrank, mat.myrow, mat.mycol);
      printf("Local Matrix:\n");
      for (int ii=0; ii<mat.m_local; ++ii) {
	for (int jj=0; jj<mat.n_local; ++jj) {
	  //printf("%3.2f ", mat.array[jj * mat.lld + ii]);
          printf("%3.2f ", mat.array[jj + ii * mat.lld]);
	}
	printf("\n");
      }
      printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

} // namespace rokko

#endif // ROKKO_EIGEN_S_DISTRIBUTED_MATRIX_HPP


