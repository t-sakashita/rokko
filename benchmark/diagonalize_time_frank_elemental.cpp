#include <mpi.h>
#include <iostream>
#include <fstream>

#include <rokko/solver.hpp>
#include <rokko/grid.hpp>
#include <rokko/distributed_matrix.hpp>
#include <rokko/collective.hpp>

#include <rokko/utility/frank_matrix.hpp>
#include <rokko/utility/sort_eigenpairs.hpp>

int main (int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  rokko::solver solver("elemental");
  solver.initialize(argc, argv);

  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid<> g(comm);
  int myrank = g.myrank;

  std::ofstream ofs;
  if (myrank == 0) {
   ofs.open("elemental_time.txt");
   if (!ofs) {
     MPI_Abort(MPI_COMM_WORLD, 22) ;
   }
  }

  const int root = 0;
  const int dim = 100;

  rokko::distributed_matrix<rokko::matrix_col_major> mat(dim, dim, g);
  rokko::generate_frank_matrix(mat);

  rokko::distributed_matrix<rokko::matrix_col_major> Z(dim, dim, g);
  Eigen::VectorXd w(dim);

  // Solve the problem
  double start, end;

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
  solver.diagonalize(mat, w, Z);
  MPI_Barrier(MPI_COMM_WORLD);
  end = MPI_Wtime();

  double time;
  if (myrank == 0) {
    time = end - start;
    std::cout << "time=" << time << std::endl;
    ofs << "time=" << time << std::endl;
  }

  solver.finalize();
  MPI_Finalize();
  return 0;
}
