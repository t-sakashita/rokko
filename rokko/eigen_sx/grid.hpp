#ifndef ROKKO_EIGEN_SX_GRID_H
#define ROKKO_EIGEN_SX_GRID_H

#include <cmath>
#include <mpi.h>
#include <boost/noncopyable.hpp>

#include <rokko/eigen_sx/eigen_sx.hpp>

#include <rokko/grid.hpp>

namespace rokko {

template<>
class grid<rokko::eigen_sx> : private boost::noncopyable
{
public:
  grid(MPI_Comm& comm)
  {
    int size_of_col_local, size_of_row_local;
    int ndims = 2;
    eigen_init_wrapper(ndims, size_of_col_local, size_of_row_local);

    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    nprow = cycl2d_.size_of_row;
    npcol = cycl2d_.size_of_col;
    myrow = cycl2d_.my_row - 1;
    mycol = cycl2d_.my_col - 1;

    if (myrank == 0) {
      cout << "gridinfo nprow=" << nprow << "  npcol=" << npcol << endl;
    }
  }

  ~grid()
  {
  }

  int myrank, nprocs;
  int myrow, mycol;
  int nprow, npcol;

private:
  int info;

};

} // namespace rokko

#endif // ROKKO_EIGEN_SX_GRID_H
