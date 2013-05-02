#ifndef ROKKO_GRID_ELEMENTAL_H
#define ROKKO_GRID_ELEMENTAL_H

#include <cmath>
#include <mpi.h>

#include <rokko/elemental/elemental.hpp>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "../grid.hpp"


namespace rokko {

namespace elemental {

  rokko::Initialize

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

    elem_grid.reset(new elem::Grid(comm, nprow, npcol));
    myrow = elem_grid->Row();
    mycol = elem_grid->Col();
  }

  boost::shared_ptr<const elem::Grid> get_elem_grid() const
  {
    return elem_grid;
  }

  ~grid()
  {
  }

  int myrank, nprocs;
  int myrow, mycol;
  int nprow, npcol;
private:
  boost::shared_ptr<elem::Grid> elem_grid;

  int info;

};

} // namespace rokko

#endif // ROKKO_GRID_ELEMENTAL_H
