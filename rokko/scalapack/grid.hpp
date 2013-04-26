#ifndef ROKKO_BLACS_GRID_H
#define ROKKO_BLACS_GRID_H

#include <cmath>
#include <mpi.h>
#include <boost/noncopyable.hpp>

#include <rokko/scalapack/blacs.hpp>

#include <rokko/grid.hpp>

namespace rokko {

template<>
class grid<rokko::scalapack> : private boost::noncopyable
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

#endif // ROKKO_BLACS_GRID_H
