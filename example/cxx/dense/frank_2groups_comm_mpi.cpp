/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2016 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
#include <rokko/mapping_bc.hpp>
#include <rokko/collective.hpp>
#include <rokko/utility/solver_name.hpp>
#include <rokko/utility/frank_matrix.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>


typedef rokko::matrix_col_major matrix_major;

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  
  MPI_Group group_world, even_group, odd_group;
  MPI_Comm comm, even_comm, odd_comm;
  int rank, p, Neven, Nodd, ierr;
  int *even_members, *odd_members;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_group(MPI_COMM_WORLD, &group_world);

  Neven = (p + 1) / 2;
  Nodd = p - Neven;
  even_members = new int[Neven];
  odd_members = new int[Nodd];
  for (int i=0; i<Neven; ++i)
    even_members[i] = 2 * i;
  for (int i=0; i<Nodd; ++i)
    odd_members[i] = 2 * i + 1;
  
  ierr = MPI_Group_incl(group_world, Neven, even_members, &even_group);
  ierr = MPI_Group_incl(group_world, Nodd, odd_members, &odd_group);
  ierr = MPI_Comm_create(MPI_COMM_WORLD, even_group, &even_comm);
  ierr = MPI_Comm_create(MPI_COMM_WORLD, odd_group, &odd_comm);
  comm = ((rank % 2) == 0) ? even_comm : odd_comm;
  MPI_Group_free(&group_world);
  MPI_Group_free(&even_group);
  MPI_Group_free(&odd_group);
  
  if (comm == MPI_COMM_NULL) {
    std::cout << "orig_rank=" << rank << " is COMM_NULL" << std::endl;
  } else {
    std::cout << "orig_rank=" << rank << " is NOT COMM_NULL" << std::endl;
  }
  rokko::grid g(comm);
  int myrank = g.get_myrank();

  if (comm != MPI_COMM_NULL) {
    std::string library_routine(rokko::parallel_dense_ev::default_solver());
    std::string library, routine;
    int dim = 10;
    if (argc >= 2) library_routine = argv[1];
    if (argc >= 3) dim = boost::lexical_cast<int>(argv[2]);
    rokko::split_solver_name(library_routine, library, routine);
    
    std::cout.precision(5);
    
    rokko::parallel_dense_ev solver(library);
    solver.initialize(argc, argv);
    if (myrank == 0)
      std::cout << "Eigenvalue decomposition of Frank matrix" << std::endl
		<< "library:routine = " << library_routine << std::endl
		<< "num_procs = " << g.get_nprocs() << std::endl
                #ifdef _OPENMP
		<< "num_threads per process = " << omp_get_max_threads() << std::endl
                #endif
		<< "routine = " << routine << std::endl
		<< "dimension = " << dim << std::endl;
    
    rokko::mapping_bc<matrix_major> map = solver.default_mapping(dim, g);
    rokko::distributed_matrix<double, matrix_major> mat(map);
    rokko::frank_matrix::generate(mat);
    rokko::localized_matrix<double, matrix_major> mat_loc(dim, dim);
    rokko::gather(mat, mat_loc, 0);
    
    rokko::localized_vector<double> eigval(dim);
    rokko::distributed_matrix<double, matrix_major> eigvec(map);
    rokko::parameters params;
    params.set("routine", routine);
    try {
      solver.diagonalize(mat, eigval, eigvec, params);
    } catch (const char *e) {
      if (myrank == 0) std::cout << "Exception : " << e << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 22);
    }
    
    rokko::localized_matrix<double, matrix_major> eigvec_loc(dim, dim);
    rokko::gather(eigvec, eigvec_loc, 0);
    if (myrank == 0) {
      bool sorted = true;
      for (unsigned int i = 1; i < dim; ++i) sorted &= (eigval(i-1) <= eigval(i));
      if (!sorted) std::cout << "Warning: eigenvalues are not sorted in ascending order!\n";
      
      std::cout << "largest eigenvalues:";
      for (int i = 0; i < std::min(dim, 10); ++i) std::cout << ' ' << eigval(dim - 1 - i);
      std::cout << std::endl;
      std::cout << "residual of the largest eigenvalue/vector: |x A x - lambda| = "
		<< std::abs(eigvec_loc.col(dim - 1).transpose() * mat_loc * eigvec_loc.col(dim - 1)
			    - eigval(dim - 1))
		<< std::endl;
    }    
    solver.finalize();
    MPI_Comm_free(&comm);
  }
  MPI_Finalize();
}

