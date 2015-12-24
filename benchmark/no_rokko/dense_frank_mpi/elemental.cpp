/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include <rokko/utility/machine_info.hpp>

//using namespace std;
//using namespace El;

int main( int argc, char* argv[] ) {
  double init_tick, initend_tick, gen_tick, diag_tick, end_tick;
  
  // This detects whether or not you have already initialized MPI and 
  // does so if necessary. The full routine is El::Initialize.
  init_tick = MPI_Wtime();
  El::Initialize( argc, argv );
  initend_tick = MPI_Wtime();
  
  // Surround the Elemental calls with try/catch statements in order to 
  // safely handle any exceptions that were thrown during execution.
  try {
    const El::Int dim = 100; //30000;
    El::ProcessInput();
    El::PrintInputReport();
    
    // Create a 2d process grid from a communicator. In our case, it is
    // MPI_COMM_WORLD. There is another constructor that allows you to 
    // specify the grid dimensions, Grid g( comm, r ), which creates a
    // grid of height r.
    El::Grid g( El::mpi::COMM_WORLD );
    
    // Create an n x n complex distributed matrix, 
    // We distribute the matrix using grid 'g'.
    //
    // There are quite a few available constructors, including ones that 
    // allow you to pass in your own local buffer and to specify the 
    // distribution alignments (i.e., which process row and column owns the
    // top-left element)
    El::DistMatrix<double> H( dim, dim, g );
    
    // Fill the matrix since we did not pass in a buffer. 
    //
    // We will fill entry (i,j) with the complex value (i+j,i-j) so that 
    // the global matrix is Hermitian. However, only one triangle of the 
    // matrix actually needs to be filled, the symmetry can be implicit.
    //
    const El::Int localHeight = H.LocalHeight();
    const El::Int localWidth = H.LocalWidth();
    for( El::Int jLoc=0; jLoc<localWidth; ++jLoc ) {
      // Our process owns the rows colShift:colStride:n,
      //           and the columns rowShift:rowStride:n
      const El::Int j = H.GlobalCol(jLoc);
      for( El::Int iLoc=0; iLoc<localHeight; ++iLoc ) {
	const El::Int i = H.GlobalRow(iLoc);
	H.SetLocal( iLoc, jLoc, dim - std::max(i, j) );
      }
    }
    gen_tick = MPI_Wtime();
    
    // Make a backup of H before we overwrite it within the eigensolver
    auto HCopy( H );
    
    // Call the eigensolver. We first create an empty complex eigenvector 
    // matrix, X[MC,MR], and an eigenvalue column vector, w[VR,* ]
    //
    // Optional: set blocksizes and algorithmic choices here. See the 
    //           'Tuning' section of the README for details.
    diag_tick = MPI_Wtime();
    El::DistMatrix<double,El::VR,El::STAR> w( g );
    El::DistMatrix<double> X( g );
    El::HermitianEig( El::LOWER, H, w, X, El::ASCENDING ); 
    end_tick = MPI_Wtime();        
    
    // Check the residual, || H X - Omega X ||_F
    const double frobH = El::HermitianFrobeniusNorm( El::LOWER, HCopy );
    auto E( X );
    El::DiagonalScale( El::RIGHT, El::NORMAL, w, E );
    El::Hemm( El::LEFT, El::LOWER, double(-1), HCopy, X, double(1), E );
    const double frobResid = El::FrobeniusNorm( E );
    
    // Check the orthogonality of X
    El::Identity( E, dim, dim );
    El::Herk( El::LOWER, El::NORMAL, double(-1), X, double(1), E );
    const double frobOrthog = El::HermitianFrobeniusNorm( El::LOWER, E );
    
    double* eigvals;
    eigvals = new double[dim];
    for (int i = 0; i < dim; ++i) {
      eigvals[i] = w.Get(i,0);
    }
    
    if( El::mpi::WorldRank() == 0 ) {
      std::cout << "|| H ||_F = " << frobH << "\n"
		<< "|| H X - X Omega ||_F / || A ||_F = " 
		<< frobResid / frobH << "\n"
		<< "|| X X^H - I ||_F = " << frobOrthog / frobH
		<< "\n" << std::endl;
      std::cout << "init_time = " << initend_tick - init_tick << std::endl
		<< "gen_time = " << diag_tick - gen_tick << std::endl
		<< "diag_time = " << end_tick - diag_tick << std::endl;
      rokko::machine_info();
      bool sorted = true;
      for (unsigned int i = 1; i < dim; ++i) sorted &= (eigvals[i-1] <= eigvals[i]);
      if (!sorted) std::cout << "Warning: eigenvalues are not sorted in ascending order!\n";
      std::cout << "largest eigenvalues:";
      for (int i = 0; i < std::min(dim, 100); ++i)
	std::cout << ' ' << eigvals[i];
      std::cout << std::endl;
    }
    
  }
  catch( std::exception& e ) { El::ReportException(e); }
  
  El::Finalize();
  return 0;
}

