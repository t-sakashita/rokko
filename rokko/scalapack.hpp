#ifndef ROKKO_SCALAPACK_H
#define ROKKO_SCALAPACK_H

#include <mpi.h>

extern "C" {
  /* BLACS */
  void blacs_pinfo_( int& mypnum, int& nprocs);
  void blacs_get_( const int& context, const int& request, int& value );
  void blacs_gridinfo_(const int& ictxt, int& nprow, int& npcol, int& myrow, int& mycol);
  void blacs_gridinit_(const int& ictxt, char* order, int& nprow, int& npcol );
  void blacs_gridexit_(int* ictxt);
  void blacs_exit_(const int& conti);
  void blacs_barrier_(const int& ictxt, const char* score);
  /* ScaLAPACK */
  void sl_init_(int, int, int);
  void descinit_(int* desc, const int& m, const int& n, const int& mb, const int& nb,
                 const int& irsrc, const int& icsrc, const int& ixtxt, const int& lld, int& info);
  void pdelset_(double* A, const int& ia, const int& ja, const int* descA, const double& alpha);
  void pdelget_(char* scope, char* top, double& alpha, const double* A, const int& ia, const int& ja, const int* descA);
  int numroc_(const int& n, const int& nb, const int& iproc, const int& isrcproc, const int& nprocs);
  void pdsyev_(char* jobz, char* uplo, const int& n,
               double* a, const int& ia, const int& ja, const int* descA,
               double* w, double* z, const int& iz, const int& jz, const int* descZ,
               double* work, const int& lwork, int& info);

  void pdsyevr_(const char* jobz, const char* uplo, const int& n,
                double* A, const int& ia, const int& ja, const int* descA,
                const double& vl, const double& vu, const int& il, const int& iu,
                const int& m, const int& nz, double* w,
                double* Z, const int& iz, const int& jz, const int* descZ,
                double* work, const int& lwork, int* iwork, const int& liwork, int& info);
  void pdsyevd_(const char* jobz, const char* uplo, const int& n,
                double* A, const int& ia, const int& ja, const int* descA,
                const double& vl, const double& vu, const int& il, const int& iu,
                const int& m, const int& nz, double* w,
                double* Z, const int& iz, const int& jz, const int* descZ,
                double* work, const int& lwork, int* iwork, const int& liwork, int& info);

  void pdlaprnt_(const int& m, const int& n, const double* A, const int& ia, const int& ja, const int* descA, const int& irprnt, const int& icprnt, char* cmatnm, const int& nout, double* work);
}

namespace rokko {
namespace scalapack {

template <class MATRIX, class VECTOR>
int diagonalize(MATRIX& mat, VECTOR& eigvals, MATRIX& eigvecs)
{
  int dim = mat.m_global;
  cout << "pdsyev_dim=" << dim << endl;

  int ictxt = mat.g.ictxt;
  const int ZERO=0, ONE=1;
  int desc[9];
  int info;

  int lld = mat.m_local;
  cout << "lld=" << lld << endl;
  if (lld == 0) lld = 1;
  descinit_(desc, mat.m_global, mat.n_global, mat.mb, mat.nb, ZERO, ZERO, ictxt, lld, info);

  if (info) {
    cerr << "error " << info << endl;
    //<< " at descinit function of descA " << "mA=" << m_local << "  nA=" << n_local << "  lld=" << lld << "." << endl;
    exit(1);
  }

  double* work = new double[1];
  long lwork = -1;

  // work配列のサイズの問い合わせ
  pdsyev_( "V",  "U",  dim,  mat.array, ONE,  ONE,  desc, &eigvals[0], eigvecs.array, ONE, ONE,
 	   desc, work, lwork, info );

  for (int i=0; i<mat.mb*mat.nb; ++i)
    cout << (mat.array)[i] << " ";
  cout << endl;

  lwork = work[0];
  delete[] work;
  work = new double [lwork];
  if (work == NULL) {
    cerr << "failed to allocate work. info=" << info << endl;
    return info;
  }
  //info = 0;

  // 固有値分解
  pdsyev_( "V",  "U",  dim,  mat.array,  ONE,  ONE,  desc, &eigvals[0], eigvecs.array, ONE, ONE,
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

} // namespace scalapack
} // namespace rokko

#endif // ROKKO_SCALAPACK_H
