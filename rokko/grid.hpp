#ifndef ROKKO_GRID_H
#define ROKKO_GRID_H

#include <mpi.h>
#include <cmath>
#include <boost/noncopyable.hpp>
#include <boost/type_traits/is_same.hpp>

namespace rokko {

struct grid_row_major {};

struct grid_col_major {};

class grid_base {
public:
  //grid_base(MPI_Comm& comm_in)  {} //: comm(comm_in) {}
  virtual int calculate_grid_row(int proc_rank) const = 0;
  virtual int calculate_grid_col(int proc_rank) const = 0;
  virtual ~grid_base() {}
  virtual bool is_row_major() const = 0;
  virtual bool is_col_major() const = 0;
  virtual MPI_Comm get_comm() const = 0;

  MPI_Comm comm;
};

template<typename GRID_MAJOR = rokko::grid_row_major>
class grid : public grid_base, private boost::noncopyable {
public:
  grid(MPI_Comm comm_in = MPI_COMM_WORLD) { //: comm(comm_in) {
    comm = comm_in;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &myrank);

    npcol = int(std::sqrt(nprocs + 0.5));
    while (1) {
      if ( npcol == 1 ) break;
      if ( (nprocs % npcol) == 0 ) break;
      npcol = npcol - 1;
    }
    nprow = nprocs / npcol;
    myrow = calculate_grid_row(myrank);
    mycol = calculate_grid_col(myrank);
  }

  ~grid()  {}

  int calculate_grid_row(int proc_rank) const;
  int calculate_grid_col(int proc_rank) const;

  int get_nprocs() const { return nprocs; }
  int get_nprow() const { return nprow; }
  int get_npcol() const { return npcol; }

  int get_myrank() const { return myrank; }
  int get_myrow() const { return myrow; }
  int get_mycol() const { return mycol; }

  bool is_row_major() const {
    return boost::is_same<GRID_MAJOR, grid_row_major>::value;
  }
  bool is_col_major() const {
    return boost::is_same<GRID_MAJOR, grid_col_major>::value;
  }

  MPI_Comm get_comm() const { return comm; }

  //MPI_Comm comm;

  int myrank, nprocs;
  int myrow, mycol;
  int nprow, npcol;
};

template<>
inline int grid<rokko::grid_row_major>::calculate_grid_row(int proc_rank) const {
  return proc_rank / npcol;
}

template<>
inline int grid<rokko::grid_row_major>::calculate_grid_col(int proc_rank) const {
  return proc_rank % npcol;
}

template<>
inline int grid<rokko::grid_col_major>::calculate_grid_row(int proc_rank) const {
  return proc_rank % nprow;
}

template<>
inline int grid<rokko::grid_col_major>::calculate_grid_col(int proc_rank) const {
  return proc_rank / nprow;
}

} // namespace rokko

#endif // ROKKO_GRID_H
