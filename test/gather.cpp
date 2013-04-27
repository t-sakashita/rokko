#include <mpi.h>

#include <Eigen/Dense>
#include <iostream>
using namespace std;


#include <rokko/scalapack/grid.hpp>
#include <rokko/scalapack/distributed_matrix.hpp>

#include <rokko/collective.hpp>
#include <rokko/utility/frank_matrix.hpp>


int main(int argc, char* argv[])
{
  typedef rokko::scalapack solver;

  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  rokko::grid<solver> g(comm);

  int myrank, nprocs;
  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm, &nprocs);

  int dim = 10;
  int root = 0;

  rokko::distributed_matrix<solver> mat(dim, dim, g);
  rokko::generate_frank_matrix(mat);
  //rokko::generate_frank_matrix_global(mat);
  Eigen::MatrixXd mat_global(dim, dim);
  for(int i=0; i<dim; ++i) {
    for(int j=0; j<dim; ++j) {
      mat_global(i, j) = dim * j + i;
    }
  }
  rokko::scatter(mat, mat_global, root);

  mat.print();
  if (myrank == root)
    cout << "global_mat:" << endl << mat_global << endl;

  MPI_Finalize();
}
