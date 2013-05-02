#ifndef ROKKO_GRID_H
#define ROKKO_GRID_H

#include <cmath>
//#include <mpi.h>
//##include "mpi.h"

#include "elemental.hpp"
//#include "elemental-lite.hpp"

#include <boost/noncopyable.hpp>


namespace rokko {

class grid : private boost::noncopyable
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

    /*
    elem::Grid grid(elem::mpi::COMM_WORLD);
    elem::DistMatrix<double> H(10, 10, grid);
    elem::DistMatrix<double,elem::VR,elem::STAR> w(grid);
    elem::DistMatrix<double> X(grid);
    //elem::DistMatrix<double,elem::VR,elem::STAR> w(10, 1, grid);
    //elem::DistMatrix<double> X(10, 10, grid);

    const int colShift = H.ColShift(); // first row we own
    const int rowShift = H.RowShift(); // first col we own
    const int colStride = H.ColStride();
    const int rowStride = H.RowStride();
    const int localHeight = H.LocalHeight();
    const int localWidth = H.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
      {
	for( int iLocal=0; iLocal<localHeight; ++iLocal )
	  {
	    // Our process owns the rows colShift:colStride:n,
	    //           and the columns rowShift:rowStride:n
	    const int i = colShift + iLocal*colStride;
	    const int j = rowShift + jLocal*rowStride;
	    H.SetLocal( iLocal, jLocal, i+j );
	  }
      }

    elem::HermitianEig(elem::LOWER, H, w, X); // only access lower half of H
    */

    //new( &elem_grid ) elem::Grid(comm, nprow, npcol);
    elem_grid = new elem::Grid(comm, nprow, npcol);
    myrow = elem_grid->Row();
    mycol = elem_grid->Col();

    /*
    if (myrank == 0) {
      cout << "gridinfo nprow=" << nprow << "  npcol=" << npcol << endl;
    }

    elem::DistMatrix<double> H(10, 10, *elem_grid);
    elem::DistMatrix<double,elem::VR,elem::STAR> w(*elem_grid);
    elem::DistMatrix<double> X(*elem_grid);
    const int colShift = H.ColShift(); // first row we own
    const int rowShift = H.RowShift(); // first col we own
    const int colStride = H.ColStride();
    const int rowStride = H.RowStride();
    const int localHeight = H.LocalHeight();
    const int localWidth = H.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
      {
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
          {
            // Our process owns the rows colShift:colStride:n,
            //           and the columns rowShift:rowStride:n
            const int i = colShift + iLocal*colStride;
            const int j = rowShift + jLocal*rowStride;
            H.SetLocal( iLocal, jLocal, i+j );
          }
      }

    elem::HermitianEig(elem::LOWER, H, w, X); // only access lower half of H
    */

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

#endif // ROKKO_GRID_H
