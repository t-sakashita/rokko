#ifndef ROKKO_SCALAPACK_DIAGONALIZE_H
#define ROKKO_SCALAPACK_DIAGONALIZE_H

#include <mpi.h>
#include <rokko/scalapack/scalapack.hpp>

namespace rokko {

template<typename T>
int diagonalize(rokko::distributed_matrix<T>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<T>& eigvecs)
{
}

template<typename T>
int diagonalize(rokko::distributed_matrix<T>& mat, double* eigvals, rokko::distributed_matrix<T>& eigvecs)
{
}

template<>
int diagonalize(rokko::distributed_matrix<rokko::scalapack>& mat, double* eigvals, rokko::distributed_matrix<rokko::scalapack>& eigvecs)
{
  int dim = mat.m_global;
  //cout << "pdsyev_dim=" << dim << endl;

  int ictxt = mat.g.ictxt;
  const int ZERO=0, ONE=1;
  int desc[9];
  int info;

  int lld = mat.m_local;
  cout << "lld=" << lld << endl;
  if (lld == 0) lld = 1;
  descinit_(desc, mat.m_global, mat.n_global, mat.mb, mat.nb, ZERO, ZERO, ictxt, lld, info);
  if (info) {
    cerr << "error " << info << " at descinit function of descA " << "mA=" << mat.m_local << "  nA=" << mat.n_local << "  lld=" << lld << "." << endl;
    MPI_Abort(MPI_COMM_WORLD, 89);
  }

  for (int proc=0; proc<mat.g.nprocs; ++proc) {
    if (proc == mat.g.myrank) {
      cout << "pdsyev:proc=" << proc << " m_global=" << mat.m_global << "  n_global=" << mat.n_global << "  mb=" << mat.mb << "  nb=" << mat.nb << " lld=" << lld << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  double* work = new double[1];
  long lwork = -1;

  // work配列のサイズの問い合わせ
  pdsyev_( "V",  "U",  dim,  mat.array, ONE,  ONE,  desc, &eigvals[0], eigvecs.array, ONE, ONE,
 	   desc, work, lwork, info );

  /*
  for (int i=0; i<mat.mb*mat.nb; ++i)
    cout << (mat.array)[i] << " ";
  cout << endl;
  */

  lwork = work[0];
  delete[] work;
  work = new double [lwork];
  if (work == NULL) {
    cerr << "failed to allocate work. info=" << info << endl;
    return info;
  }

  // 固有値分解
  pdsyev_( "V",  "U",  dim,  mat.array,  ONE,  ONE,  desc, eigvals, eigvecs.array, ONE, ONE,
  	   desc, work, lwork, info );
  /*
  pdsyev_( "V",  "U",  dim,  mat.array,  ONE,  ONE,  desc, eigvals.data(), eigvecs.array, ONE, ONE,
	   desc, work, lwork, info );
  */

  if (info) {
    cerr << "error at pdsyev function. info=" << info  << endl;
    exit(1);
  }

  delete[] work;
  return info;
}


template<>
int diagonalize(rokko::distributed_matrix<rokko::scalapack>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<rokko::scalapack>& eigvecs)
{
  return diagonalize(mat, &eigvals[0], eigvecs);
}


} // namespace rokko

#endif // ROKKO_SCALAPACK_DIAGONALIZE_H
