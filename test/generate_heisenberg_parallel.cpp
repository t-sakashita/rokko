/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <iostream>
#include <vector>
#include <rokko/utility/heisenberg_hamiltonian_parallel.hpp>
#include <rokko/localized_matrix.hpp>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  int L = 5;
  
  std::vector<std::pair<int, int> > lattice;
  for (int i=0; i<L-1; ++i) {
    lattice.push_back(std::make_pair(i, i+1));
  }

  int myrank, nproc;
  MPI_Status status;
  int ierr;

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int n = nproc;
  int p = -1;
  do {
    n /= 2;
    ++p;
  } while (n > 0);

  //std::cerr << "nproc=" << nproc << " p=" << p << std::endl;
  if (nproc != (1 << p)) {    
    if ( myrank == 0 ) {
      std::cout << "This program can be run only for powers of 2" << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  int N = 1 << (L-p);

  std::vector<double> buffer(N);
  for (int proc = 0; proc<nproc; ++proc) {
    for (int i=0; i<N; ++i) {
      std::vector<double> v, w;
      v.assign(N, 0);
      if (myrank == proc) v[i] = 1;
      w.assign(N, 0);
      rokko::heisenberg_hamiltonian::multiply(MPI_COMM_WORLD, L, lattice, v, w, buffer);
      //std::cout << "w=";
      for (int proc = 0; proc<nproc; ++proc) {
	if (proc == 0) {
	  if (myrank == proc) {
	    for (int j=0; j<N; ++j) {
	      std::cout << w[j] << " ";
	    }
	  }
	} else {
	  if (myrank == proc) {
	    MPI_Send(&w[0], N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	  }
	  if (myrank == 0) {
            MPI_Recv(&buffer[0], N, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &status);
	    for (int j=0; j<N; ++j) {
	      std::cout << buffer[j] << " ";
	    }
	  }
	}
      }
      if (myrank == 0) std::cout << std::endl;
    }
  }

  MPI_Finalize();
  return 0;
}


