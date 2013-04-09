#include <mpi.h>

namespace rokko {

class grid
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
    
    nprow = int(sqrt(nprocs+0.5));
    while (1) {
      if ( nprow == 1 ) break;
      if ( (nprocs % nprow) == 0 ) break;
      nprow = nprow - 1;
    }                                                                                                                                                                           
    npcol = nprocs / nprow;
    blacs_gridinit_( ictxt, "C", nprow, npcol ); // ColがMPI_Comm_createと互換
    //blacs_gridinit_( ictxt, "Row", nprow, npcol );
    blacs_gridinfo_( ictxt, nprow, npcol, myrow, mycol );
    
    if (myrank == 0) {
      cout << "gridinfo nprow=" << nprow << "  npcol=" << npcol << "  ictxt=" << ictxt << endl;
    }
  }

  /*
  ~grid()
  {
    //blacs_gridexit_(&ictxt);
  }
  */

  int myrank, nprocs;
  int myrow, mycol;
  int nprow, npcol;
  int ictxt;
private:

  int info;

};

} // namespace rokko
