#ifndef ROKKO_EIGEN_S_GRID_H
#define ROKKO_EIGEN_S_GRID_H

#include <cmath>
#include <mpi.h>
#include <boost/noncopyable.hpp>

#include <rokko/eigen_s/eigen_s.hpp>
#include <rokko/scalapack/blacs.hpp>

#include <rokko/grid.hpp>

namespace rokko {

template<>
class grid<rokko::eigen_s> : private boost::noncopyable
{
public:
  grid(MPI_Comm& comm)
  {
    int size_of_col_local, size_of_row_local;
    int ndims = 2;
    eigen_init_wrapper(ndims, size_of_col_local, size_of_row_local);

    /*
    int nprow = cycl2d_.size_of_row;
    int npcol = cycl2d_.size_of_col;
    int myrow = cycl2d_.my_row;
    int mycol = cycl2d_.my_col;
    //cout << "NPROW=" << cycl2d_.size_of_row << "  NPCOL=" << cycl2d_.size_of_col  << endl;
    */

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

#endif // ROKKO_EIGEN_S_GRID_H
