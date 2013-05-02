#ifndef ROKKO_SCALAPACK_DIAGONALIZE_H
#define ROKKO_SCALAPACK_DIAGONALIZE_H

#include <mpi.h>
#include <rokko/scalapack/scalapack.hpp>

namespace rokko {

/*
template<typename T, typename GRID_MAJOR>
int diagonalize(rokko::distributed_matrix<T, GRID_MAJOR>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<T, GRID_MAJOR>& eigvecs)
//int diagonalize(rokko::distributed_matrix& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix& eigvecs)
{
  cout << "pp100" << endl;
}
*/
 /*
template<typename GRID_MAJOR>
int diagonalize(rokko::distributed_matrix<rokko::scalapack, GRID_MAJOR>& mat, double* eigvals, rokko::distributed_matrix<rokko:scalapack, GRID_MAJOR>& eigvecs)
{
  cout << "pp2"<< endl;
}

template<typename GRID_MAJOR = rokko::R>
int diagonalize(rokko::distributed_matrix<GRID_MAJOR>& mat, Eigen::VectorXd& eigvals)
{
  cout << "pp3"<< endl;
}
 */
/*
template<typename GRID_MAJOR = rokko::R>
int diagonalize(rokko::distributed_matrix<rokko:scalapack, GRID_MAJOR>& mat, double* eigvals)
{
  cout << "pp4"<< endl;
}
*/


 //template<typename GRID_MAJOR = rokko::R>
int diagonalize(rokko::distributed_matrix<rokko::scalapack>& mat, double* eigvals)
{
  cout << "pp4"<< endl;

  int dim = mat.m_global;
  //cout << "pdsyev_dim=" << dim << endl;

  int ictxt = mat.ictxt;
  //int ictxt = grid<rokko::scalapack, rokko::grid_major_type<rokko::scalapack> >(mat.g).ictxt;

  const int ZERO=0, ONE=1;
  int desc[9];
  int info;

  int lld = mat.m_local;
  double * double_null_ptr = NULL;

  cout << "lld=" << lld << endl;
  if (lld == 0) lld = 1;
  descinit_(desc, mat.m_global, mat.n_global, mat.mb, mat.nb, ZERO, ZERO, ictxt, lld, info);
  if (info) {
    cerr << "error " << info << " at descinit function of descA " << "mA=" << mat.m_local << "  nA=" << mat.n_local << "  lld=" << lld << "." << endl;
    MPI_Abort(MPI_COMM_WORLD, 89);
  }

  for (int proc=0; proc<mat.nprocs; ++proc) {
    if (proc == mat.myrank) {
      cout << "pdsyev:proc=" << proc << " m_global=" << mat.m_global << "  n_global=" << mat.n_global << "  mb=" << mat.mb << "  nb=" << mat.nb << " lld=" << lld << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  double* work = new double[1];
  long lwork = -1;

  // work配列のサイズの問い合わせ
  pdsyev_( "N",  "U",  dim,  mat.array, ONE,  ONE,  desc, eigvals, double_null_ptr, ONE, ONE,
 	   desc, work, lwork, info );

  lwork = work[0];
  delete[] work;
  work = new double [lwork];
  if (work == NULL) {
    cerr << "failed to allocate work. info=" << info << endl;
    return info;
  }

  // 固有値分解
  pdsyev_( "N",  "U",  dim,  mat.array,  ONE,  ONE,  desc, eigvals, double_null_ptr, ONE, ONE,
  	   desc, work, lwork, info );

  if (info) {
    cerr << "error at pdsyev function. info=" << info  << endl;
    exit(1);
  }

  delete[] work;
  return info;
}

//template<typename GRID_MAJOR = rokko::R>
int diagonalize(rokko::distributed_matrix<rokko::scalapack>& mat, double* eigvals, rokko::distributed_matrix<rokko::scalapack>& eigvecs)
{
  cout << "pp6"<< endl;

  int dim = mat.m_global;
  //cout << "pdsyev_dim=" << dim << endl;

  int ictxt = mat.ictxt;
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

  for (int proc=0; proc<mat.nprocs; ++proc) {
    if (proc == mat.myrank) {
      cout << "pdsyev:proc=" << proc << " m_global=" << mat.m_global << "  n_global=" << mat.n_global << "  mb=" << mat.mb << "  nb=" << mat.nb << " lld=" << lld << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  double* work = new double[1];
  long lwork = -1;

  // work配列のサイズの問い合わせ
  pdsyev_( "V",  "U",  dim,  mat.array, ONE,  ONE,  desc, eigvals, eigvecs.array, ONE, ONE,
 	   desc, work, lwork, info );

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

  if (info) {
    cerr << "error at pdsyev function. info=" << info  << endl;
    exit(1);
  }

  delete[] work;
  return info;
}

//template<class GRID_MAJOR = rokko::R>
int diagonalize(rokko::distributed_matrix<rokko::scalapack>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<rokko::scalapack>& eigvecs)
{
  cout << "pp7"<< endl;

  return diagonalize(mat, &eigvals[0], eigvecs);
}

/*
template<typename GRID_MAJOR>  // = rokko::R>
int diagonalize(rokko::distributed_matrix<rokko::scalapack, GRID_MAJOR>& mat, Eigen::VectorXd& eigvals, rokko::distributed_matrix<rokko::scalapack, GRID_MAJOR>& eigvecs)
{
  cout << "pp7"<< endl;

  return diagonalize(mat, &eigvals[0], eigvecs);
}
*/

 //template<typename GRID_MAJOR = rokko::R>
int diagonalize(rokko::distributed_matrix<rokko::scalapack>& mat, Eigen::VectorXd& eigvals)
{
  cout << "pp8"<< endl;

  return diagonalize(mat, &eigvals[0]);
}


} // namespace rokko

#endif // ROKKO_SCALAPACK_DIAGONALIZE_H
