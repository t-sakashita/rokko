#ifndef ROKKO_GRID_ELEMENTAL_H
#define ROKKO_GRID_ELEMENTAL_H

#include <cmath>
//#include <mpi.h>
#include <new>

#include "elemental.hpp"

#include <boost/noncopyable.hpp>


#include "grid.hpp"

typedef struct
{
} elemental;

namespace rokko {

  template<>
  class grid<elemental> : private boost::noncopyable
{
public:
  grid(MPI_Comm& comm) //: elem_grid(comm)
  {
    MPI_Comm_size(comm, &nprocs);
    //nprocs = MPI::COMM_WORLD.Get_rank();
    MPI_Comm_rank(comm, &myrank);

    npcol = int(sqrt(nprocs + 0.5));
    while (1) {
      if ( npcol == 1 ) break;
      if ( (nprocs % npcol) == 0 ) break;
      npcol = npcol - 1;
    }
    nprow = nprocs / npcol;

    elem_grid = new elem::Grid(comm, nprow, npcol);
    myrow = elem_grid->Row();
    mycol = elem_grid->Col();
  }

  ~grid()
  {
  }

  int myrank, nprocs;
  int myrow, mycol;
  int nprow, npcol;
  elem::Grid* elem_grid;
private:

  int info;

};

} // namespace rokko

#endif // ROKKO_GRID_ELEMENTAL_H
