#ifndef ROKKO_GRID_H
#define ROKKO_GRID_H

#include <cmath>
#include <boost/noncopyable.hpp>

namespace rokko {

template<typename T>
class grid : private boost::noncopyable
{
public:
  grid(MPI_Comm& comm)
  {
  }

  ~grid()
  {
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
