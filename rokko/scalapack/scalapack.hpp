#ifndef ROKKO_SCALAPACK_HPP
#define ROKKO_SCALAPACK_HPP

extern "C" {

void pdelset_(double* A, const int* ia, const int* ja, const int* descA,
              const double* alpha);

inline void ROKKO_pdelset(double* A, int ia, int ja, const int* descA, double alpha) {
  pdelset_(A, &ia, &ja, descA, &alpha);
}

void pdelget_(const char* scope, const char* top, double* alpha,
              const double* A, const int* ia, const int* ja, const int* descA);

inline void ROKKO_pdelget(char scope, char top, double* alpha,
                          const double* A, int ia, int ja, const int* descA) {
  pdelget_(&scope, &top, alpha, A, &ia, &ja, descA);
}

void pdsyev_(const char* jobz, const char* uplo, const int* n,
             double* A, const int* ia, const int* ja, const int* descA,
             double* w, double* Z, const int* iz, const int* jz, const int* descZ,
             double* work, const int* lwork, int* info);

inline void ROKKO_pdsyev(char jobz, char uplo, int n,
                         double* A, int ia, int ja, const int* descA,
                         double* w, double* Z, int iz, int jz, const int* descZ,
                         double* work, int lwork, int* info) {
  pdsyev_(&jobz, &uplo, &n, A, &ia, &ja, descA, w, Z, &iz, &jz, descZ,
          work, &lwork, info);
}

void pdsyevr_(const char* jobz, const char* uplo, const int* n,
              double* A, const int* ia, const int* ja, const int* descA,
              const double* vl, const double* vu, const int* il, const int* iu,
              int* m, int* nz, double* w,
              double* Z, const int* iz, const int* jz, const int* descZ,
              double* work, const int* lwork, int* iwork, const int* liwork, int* info);

inline void ROKKO_pdsyevr(char jobz, char uplo, int n,
                          double* A, int ia, int ja, const int* descA,
                          double vl, double vu, int il, int iu,
                          int* m, int* nz, double* w,
                          double* Z, int iz, int jz, const int* descZ,
                          double* work, int lwork, int* iwork, int liwork, int* info) {
  pdsyevr_(&jobz, &uplo, &n, A, &ia, &ja, descA, &vl, &vu, &il, &iu, m, nz, w,
           Z, &iz, &jz, descZ, work, &lwork, iwork, &liwork, info);
}
  
void pdsyevd_(const char* jobz, const char* uplo, const int* n,
              double* A, const int* ia, const int* ja, const int* descA,
              double* w,
              double* Z, const int* iz, const int* jz, const int* descZ,
              double* work, const int* lwork, int* iwork, const int* liwork, int* info);

inline void ROKKO_pdsyevd(char jobz, char uplo, int n,
                          double* A, int ia, int ja, const int* descA,
                          double* w, double* Z, int iz, int jz, const int* descZ,
                          double* work, int lwork, int* iwork, int liwork, int* info) {
  pdsyevd_(&jobz, &uplo, &n, A, &ia, &ja, descA, w, Z, &iz, &jz, descZ,
           work, &lwork, iwork, &liwork, info);
}

void pdsyevx_(const char* jobz, const char* range, const char* uplo, const int* n,
              double* A, const int* iA, const int* jA, const int* descA,
              const double* vl, const double* vu, const int* il, const int* iu,
              const double* abstol, int* m, int* nZ, double* w, const double* orfac,
              double* Z, const int* iZ, const int* jZ, const int* descZ,
              double* work, const int* lwork, int* iwork, const int* liwork,
              int* ifail, int* iclustr, double* gap, int* info);

inline void ROKKO_pdsyevx(char jobz, char range, char uplo, int n,
                          double* A, int iA, int jA, const int* descA,
                          double vl, double vu, int il, int iu,
                          double abstol, int* m, int* nZ, double* w, double orfac,
                          double* Z, int iZ, int jZ, const int* descZ,
                          double* work, int lwork, int* iwork, int liwork,
                          int* ifail, int* iclustr, double* gap, int* info) {
  pdsyevx_(&jobz, &range, &uplo, &n, A, &iA, &jA, descA, &vl, &vu, &il, &iu,
           &abstol, m, nZ, w, &orfac, Z, &iZ, &jZ, descZ, work, &lwork, iwork, &liwork,
           ifail, iclustr, gap, info);
}

void pdlaprnt_(const int* m, const int* n, const double* A, const int* ia, const int* ja, const int* descA, const int* irprnt, const int* icprnt, const char* cmatnm, const int* nout, double* work);

inline void ROKKO_pdlaprnt(int m, int n, const double* A, int ia, int ja, const int* descA, int irprnt, int icprnt, const char* cmatnm, int nout, double* work) {
  pdlaprnt_(&m, &n, A, &ia, &ja, descA, &irprnt, &icprnt, cmatnm, &nout, work);
}

double pdlamch_(const int* icnt, const char* cmch);

inline double ROKKO_pdlamch(int icnt, char cmch) {
  return pdlamch_(&icnt, &cmch);
}

}

#endif // ROKKO_SCALAPACK_HPP
