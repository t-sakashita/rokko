#ifndef ROKKO_BLACS_GRID_H
#define ROKKO_BLACS_GRID_H

#include <cmath>
#include <mpi.h>
#include <boost/noncopyable.hpp>

#include <rokko/scalapack/blacs.hpp>
#include <rokko/scalapack/scalapack.hpp>

//#include <rokko/grid.hpp>

namespace rokko {

template<typename T>
struct grid_row_major
{
};


template<typename T>
struct grid_col_major
{
};


template<>
struct grid_row_major<rokko::scalapack>
{
public:
  grid_row_major<rokko::scalapack>() : major_char('R') {}
  char major_char;
};

template<>
struct grid_col_major<rokko::scalapack>
{
public:
  grid_col_major<rokko::scalapack>() : major_char('C') {}
  char major_char;
};

template<typename T, typename GRID_MAJOR = rokko::grid_row_major<T> >
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
private:

  int info;

};

template<typename GRID_MAJOR> // = rokko::R>
class grid<rokko::scalapack, GRID_MAJOR> : private boost::noncopyable
{
public:
  //template<typename GRID_MAJOR>
  //grid<rokko::scalapack, GRID_MAJOR>(MPI_Comm& comm)
  //{
  //}

  //template<typename MAJOR>
  //struct char_major
  //{
  //  char major;
  //}

  /*
  template<typename MAJOR>
  struct return_char_major<rokko::R>
  {
  public:
    return_char_major(void) : char_major('R') {}
    char char_major;
  };
  */


  //template<>
  //grid<rokko::scalapack, rokko::R>(MPI_Comm& comm)
  grid<rokko::scalapack, GRID_MAJOR>(MPI_Comm& comm)
  {
    GRID_MAJOR my_grid_major;
    cout << "GRID_MAJOR_char=" << my_grid_major.major_char << endl;

    //return_char_major<GRID_MAJOR> my_return_char_major;
    //char major2 = my_return_char_major.char_major;

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
    blacs_gridinit_( ictxt, &my_grid_major.major_char, nprow, npcol ); // ColがMPI_Comm_createと互換
    //blacs_gridinit_( ictxt, "Row", nprow, npcol );
    blacs_gridinfo_( ictxt, nprow, npcol, myrow, mycol );

    if (myrank == 0) {
      cout << "gridinfo nprow=" << nprow << "  npcol=" << npcol << "  ictxt=" << ictxt << endl;
    }
  }

  /*
  //template<>
  grid<rokko::scalapack, rokko::C>(MPI_Comm& comm)
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

    blacs_gridinit_( ictxt, "C", nprow, npcol );
    // ColがMPI_Comm_createと互換
    //blacs_gridinit_( ictxt, "Row", nprow, npcol );
    blacs_gridinfo_( ictxt, nprow, npcol, myrow, mycol );
    if (myrank == 0) {
      cout << "gridinfo nprow=" << nprow << "  npcol=" << npcol << "  ictxt=" << ictxt << endl;
    }
  }
  */

  ~grid()
  {
    blacs_gridexit_(&ictxt);
  }

  int myrank, nprocs;
  int myrow, mycol;
  int nprow, npcol;
  int ictxt;
  typedef GRID_MAJOR grid_major_type;

private:

  int info;

};

} // namespace rokko

#endif // ROKKO_BLACS_GRID_H
