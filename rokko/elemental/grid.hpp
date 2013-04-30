#ifndef ROKKO_GRID_ELEMENTAL_H
#define ROKKO_GRID_ELEMENTAL_H

#include <cmath>
#include <mpi.h>

#include <rokko/elemental/elemental.hpp>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "../grid.hpp"


namespace rokko {

  //namespace elemental {

template<>
class grid<rokko::elemental> : private boost::noncopyable
{
public:
  /*
  grid(MPI_Comm& comm, int nprow, int npcol) : comm_(comm), nprow(nprow), npcol(npcol)
  {
    //MPI_Comm_size(comm_, &nprocs);
    //MPI_Comm_rank(comm_, &myrank);

    //npcol = int(sqrt(nprocs + 0.5));
    //while (1) {
    //  if ( npcol == 1 ) break;
    //  if ( (nprocs % npcol) == 0 ) break;
    //  npcol = npcol - 1;
    //}
    //nprow = nprocs / npcol;


    elem_grid.reset(new elem::Grid(comm_, nprow, npcol));
    myrow = elem_grid->Row();
    mycol = elem_grid->Col();
  }
  */

  grid(MPI_Comm& comm) : comm_(comm)
  {
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myrank);

    npcol = int(sqrt(nprocs + 0.5));
    while (1) {
      if ( npcol == 1 ) break;
      if ( (nprocs % npcol) == 0 ) break;
      npcol = npcol - 1;
    }
    nprow = nprocs / npcol;

    //grid(comm, nprow, npcol);

    elem_grid.reset(new elem::Grid(comm, nprow, npcol));
    if(!elem_grid) {
      std::cerr << "error: create grid" << std::endl;
      MPI_Finalize();
      exit(1);
    }

    myrow = elem_grid->Row();
    mycol = elem_grid->Col();

  }

  /*
  grid()
  {  MPI_Comm comm2 = MPI_COMM_WORLD;

    grid(comm2); //MPI_COMM_WORLD);
  }
  */

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
  MPI_Comm comm_;
  boost::shared_ptr<elem::Grid> elem_grid;

  int info;
};

} // namespace rokko

#endif // ROKKO_GRID_ELEMENTAL_H
