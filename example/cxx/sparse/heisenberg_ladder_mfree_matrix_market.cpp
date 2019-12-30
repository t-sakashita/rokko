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
#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>
#include <rokko/utility/lattice.hpp>
#include <rokko/utility/math.hpp>
#include <rokko/distributed_mfree_to_crs.hpp>

class heisenberg_op : public rokko::distributed_mfree {
public:
  heisenberg_op(int L, const std::vector<std::pair<int, int>>& lattice)
    : L_(L), lattice_(lattice) {
    comm_ = MPI_COMM_WORLD;
    int size, rank;
    MPI_Comm_size(comm_, &size);
    MPI_Comm_rank(comm_, &rank);
    const int p = rokko::find_power_of_two(size);
    dim_ = 1 << L;
    num_local_rows_ = 1 << (L-p);
    local_offset_ = num_local_rows_ * rank;
    buffer_.assign(num_local_rows_, 0);
  }
  ~heisenberg_op() {}

  void multiply(const double* x, double* y) const {
    rokko::heisenberg_hamiltonian::multiply(comm_, L_, lattice_, x, y, buffer_.data());
  }
  int get_dim() const { return dim_; }
  int get_local_offset() const { return local_offset_; }
  int get_num_local_rows() const { return num_local_rows_; }

private:
  MPI_Comm comm_;
  int L_;
  std::vector<std::pair<int, int>> lattice_;
  int dim_, local_offset_, num_local_rows_;
  mutable std::vector<double> buffer_;
};

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  std::string library(rokko::parallel_sparse_ev::default_solver());
  if (argc >= 2) library = argv[1];
  
  int len_ladder = (argc >= 3) ? boost::lexical_cast<int>(argv[2]) : 5;
  int L = 2 * len_ladder;
  int dim = 1 << L;
  std::vector<std::pair<int, int>> lattice;
  rokko::create_ladder_lattice_1dim(len_ladder, lattice);

  rokko::parameters params;
  params.set("Block Size", 5);
  params.set("Maximum Iterations", 500);
  params.set("Convergence Tolerance", 1.0e-8);
  params.set("num_eigenvalues", 10);
  rokko::parallel_sparse_ev solver(library);
  heisenberg_op op(L, lattice);
  rokko::distributed_crs_matrix mat(dim, dim, solver);
  rokko::distributed_mfree_to_crs(op, mat);
  mat.output_matrix_market();
  //mat.print();
  
  solver.finalize();
  MPI_Finalize();
}
