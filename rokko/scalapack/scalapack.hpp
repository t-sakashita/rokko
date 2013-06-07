#ifndef ROKKO_SCALAPACK_HPP
#define ROKKO_SCALAPACK_HPP

extern "C" {
  /* ScaLAPACK */
  void sl_init_(int, int, int);
  void descinit_(int* desc, const int& m, const int& n, const int& mb, const int& nb,
                 const int& irsrc, const int& icsrc, const int& ixtxt, const int& lld, int& info);
  void pdelset_(double* A, const int& ia, const int& ja, const int* descA, const double& alpha);
  void pdelget_(char* scope, char* top, double& alpha, const double* A, const int& ia, const int& ja, const int* descA);
  int numroc_(const int& n, const int& nb, const int& iproc, const int& isrcproc, const int& nprocs);
  void pdsyev_(const char* jobz, const char* uplo, const int& n,
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
                //                const double& vl, const double& vu, const int& il, const int& iu,
                double* w,
                double* Z, const int& iz, const int& jz, const int* descZ,
                double* work, const int& lwork, int* iwork, const int& liwork, int& info);

  void pdsyevx_(const char* jobz, const char* range, const char* uplo, const int& n,
                double* A, const int& iA, const int& jA, const int* descA,
                const double& vl, const double& vu, const int& il, const int& iu,
                const double& abstol, const int& m, const int& nZ, double* w, const double& orfac,
                double* Z, const int& iZ, const int& jZ, const int* descZ,
                double* work, const int& lwork, int* iwork, const int& liwork,
                int* ifail, int* iclustr, double* gap, int& info);

  void pdlaprnt_(const int& m, const int& n, const double* A, const int& ia, const int& ja, const int* descA, const int& irprnt, const int& icprnt, char* cmatnm, const int& nout, double* work);
  double pdlamch_(const int&, const char&);

}


#endif // ROKKO_SCALAPACK_HPP

