#ifndef ROKKO_DISTRIBUTED_H
#define ROKKO_DISTRIBUTED_H


#include <cstdlib>
//#include <mpi.h>
#include <rokko/grid.hpp>

#include "elemental.hpp"

namespace rokko {

class distributed_matrix
{
public:
  // torus-wrap distribution used in EigenK and Elemental
  distributed_matrix(int m_global, int n_global, const grid& g)
    : m_global(m_global), n_global(n_global), g(g), myrank(g.myrank), nprocs(g.nprocs), myrow(g.myrow), mycol(g.mycol), nprow(g.nprow), npcol(g.npcol), mat(m_global, n_global, *(g.elem_grid))
  {
    //m_local = LocalLength(m_global, myrow, nprow);
    //n_local = LocalLength(n_global, mycol, npcol);
    array = mat.Buffer();
    m_local = mat.LocalHeight();
    n_local = mat.LocalWidth();
    lld = mat.LDim();
    mb = 1;
    nb = 1;
  }

  ~distributed_matrix()
  {
    //mat.~DistMatrix();
    //cout << "Destructor ~DistMatrix" << endl;
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
    //array[local_j + local_i * lld] = value;
    mat.SetLocal(local_i, local_j, value);
  }

  void set_global(int global_i, int global_j, double value)
  {
    //if ((is_gindex_myrow(global_i)) && (is_gindex_mycol(global_j)))
    //  set_local(translate_g2l_row(global_i), translate_g2l_col(global_j), value);
    mat.Set(global_i, global_j, value); /// Todo::fix
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

  void print()
  {
    /* each proc prints it's local_array out, in order */
    for (int proc=0; proc<nprocs; ++proc) {
      if (proc == myrank) {
	printf("Rank = %d  myrow=%d mycol=%d\n", myrank, myrow, mycol);
	printf("Local Matrix:\n");
	//elem::Matrix<double> Y = mat.LocalMatrix();
	mat.Matrix().Print("");
	/*
	for (int ii=0; ii<m_local; ++ii) {
	  for (int jj=0; jj<n_local; ++jj) {
	    //printf("%3.2f ", array[jj * lld + ii]);
            //printf("%3.2f ", array[jj + ii * lld]);
            printf("%e ", array[jj * lld + ii]);
	  }
	  printf("\n");
	}
	*/
	printf("\n");
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  elem::DistMatrix<double> mat;
  int m_global, n_global;
  int lld;
  double* array;
  int mb, nb;
  int m_local, n_local;
  // variables of class Grid
  int myrank, nprocs;
  int myrow, mycol;
  int ictxt, nprow, npcol;
  const grid& g;

private:
  int info;
};


void print_matrix(const rokko::distributed_matrix& mat)
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

#endif // ROKKO_DISTRIBUTED_H

