#include <mpi.h>
#include <iostream>

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/localized_vector.hpp>

#include <rokko/utility/frank_matrix.hpp>
#include <boost/lexical_cast.hpp>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  typedef rokko::grid_row_major grid_major;
  //typedef rokko::grid_col_major grid_major;
  //typedef rokko::matrix_col_major matrix_major;
  typedef rokko::matrix_row_major matrix_major;

  if (argc <= 2) {
    std::cerr << "error: " << argv[0] << " solver_name matrix_size" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 34);
  }

  std::string solver_name(argv[1]);
  unsigned int dim = boost::lexical_cast<unsigned int>(argv[2]);

  rokko::solver solver(solver_name);
  solver.initialize(argc, argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid<grid_major> g(comm);
  int myrank = g.myrank;

  const int root = 0;

  rokko::distributed_matrix<matrix_major> mat(dim, dim, g, solver);
  rokko::generate_frank_matrix(mat);

  rokko::localized_vector w(dim);
  rokko::distributed_matrix<matrix_major> Z(dim, dim, g, solver);

  double start, end;
  try {
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    solver.diagonalize(mat, w, Z);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();
  }
  catch (const char *e) {
    std::cout << "Exception : " << e << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 22);
  }

  std::cout.precision(20);

  if (myrank == root) {
    std::cout << "Computed Eigenvalues= " << w.transpose() << std::endl;
  }

  if (myrank == 0) {
    double time = end - start;
    std::cout << "time = " << time << std::endl;
  }

  solver.finalize();
  MPI_Finalize();
  return 0;
}
