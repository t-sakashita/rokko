#ifndef ROKKO_GRID_H
#define ROKKO_GRID_H

#include <cmath>
#include <mpi.h>
#include <boost/noncopyable.hpp>

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

  template<typename T>
  class grid : private boost::noncopyable
{
public:
  grid(MPI_Comm& comm)
  {
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    const int ZERO=0, ONE=1;
    long MINUS_ONE = -1;
    blacs_pinfo_( myrank, nprocs );
    blacs_get_( MINUS_ONE, ZERO, ictxt );

    npcol = int(sqrt(nprocs+0.5));
    while (1) {
      if ( npcol == 1 ) break;
      if ( (nprocs % npcol) == 0 ) break;
      npcol = npcol - 1;
    }
    nprow = nprocs / npcol;
    blacs_gridinit_( ictxt, "R", nprow, npcol ); // ColがMPI_Comm_createと互換
    //blacs_gridinit_( ictxt, "Row", nprow, npcol );
    blacs_gridinfo_( ictxt, nprow, npcol, myrow, mycol );

    if (myrank == 0) {
      cout << "gridinfo nprow=" << nprow << "  npcol=" << npcol << "  ictxt=" << ictxt << endl;
    }
  }

  ~grid()
  {
    //blacs_gridexit_(&ictxt);
  }

  int myrank, nprocs;
  int myrow, mycol;
  int nprow, npcol;
  int ictxt;
private:

  int info;

};

} // namespace rokko

#endif // ROKKO_GRID_H
