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
#include <fstream>
#include <vector>
#include <boost/tuple/tuple.hpp>

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_matrix.hpp>
#include <rokko/localized_vector.hpp>
#include <rokko/utility/xyz_hamiltonian.hpp>
#include <rokko/utility/xyz_hamiltonian_mpi.hpp>

#include <rokko/collective.hpp>

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  std::string solver_name("scalapack");
  rokko::solver solver(solver_name);
  solver.initialize(argc, argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid g(comm, rokko::grid_col_major);
  const int root = 0;

  if (argc <= 1) {
    std::cerr << "error: " << argv[0] << " path_to_xyz.ip" << std::endl;
    exit(1);
  }

  std::ifstream ifs(argv[1]);
  if (!ifs) {
    std::cout << "can't open file" << std::endl;
    exit(2);
  }

  int L, num_bonds;
  std::vector<std::pair<int, int> > lattice;
  std::vector<boost::tuple<double, double, double> > coupling;
  ifs >> L >> num_bonds;
  for (int i=0; i<num_bonds; ++i) {
    int j, k;
    ifs >> j >> k;
    lattice.push_back(std::make_pair(j, k));
  }
  
  for (int i=0; i<num_bonds; ++i) {
    double jx, jy, jz;
    ifs >> jx >> jy >> jz;
    coupling.push_back(boost::make_tuple(jx, jy, jz));
  }
  
  std::cout << "L=" << L << " num_bonds=" << num_bonds << std::endl;
  for (int i=0; i<num_bonds; ++i) {
    std::cout << lattice[i].first << " " << lattice[i].second << " " << coupling[i].get<0>() << " " << coupling[i].get<1>() << " " << coupling[i].get<2>() << std::endl;
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

  // creating column vectors which forms a heisenberg hamiltonian.
  int N_global = 1 << L;
  rokko::localized_matrix<rokko::matrix_col_major> lmat0(N_global, N_global);
  std::vector<double> buffer(N);
  for (int i=0; i<N; ++i) {
    std::vector<double> v, w;
    v.assign(N, 0);
    v[i] = 1;
    w.assign(N, 0);
    rokko::xyz_hamiltonian::multiply(MPI_COMM_WORLD, L, lattice, coupling, v, w, buffer);
    for (int j=0; j<N; ++j) {
      //lmat0(myrank * N + j, myrank * N + i) = w[j];
      std::cout << w[j] << " ";
    }
    std::cout << std::endl;
  }

  // create a hamiltonian as a distributed_matrix and gather to compare a matrix genereated by serial version
  rokko::distributed_matrix<rokko::matrix_col_major> mat(N_global, N_global, g, solver);
  rokko::xyz_hamiltonian::generate(L, lattice, coupling, mat);

  rokko::localized_matrix<rokko::matrix_col_major> lmat1(N_global, N_global);
  rokko::gather(mat, lmat1, root);

  // create a matrix genereated by serail version
  rokko::localized_matrix<rokko::matrix_col_major> lmat2(N_global, N_global);
  rokko::xyz_hamiltonian::generate(L, lattice, coupling, lmat2);

  //mat.print();
  if (myrank == root) {
    if (lmat1 == lmat2) {
      std::cout << "OK: matrix by MPI version 'generate' equals to a matrix by serial version 'genertate'." << std::endl;
    } else {
      std::cout << "ERROR: matrix by MPI version 'generate' is differnet from a matrix by serial version 'genertate'." << std::endl;
      std::cout << "lmat1=" << std::endl << lmat1 << std::endl;
      std::cout << "lmat2=" << std::endl << lmat2 << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  MPI_Finalize();
  return 0;
}


