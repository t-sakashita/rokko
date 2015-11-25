/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2015 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#include <rokko/rokko.hpp>
//#include <rokko/utility/heisenberg_hamiltonian.hpp>
#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>

void distributed_mfree_to_crs(rokko::distributed_mfree const& op, rokko::distributed_crs_matrix& mat) {
  std::vector<double> x(op.get_num_local_rows()), y(op.get_num_local_rows());
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  for (int global_row=0; global_row<op.get_dim(); ++global_row) {
    x.assign(op.get_num_local_rows(), 0.);
    y.assign(op.get_num_local_rows(), 0.);
    if ((global_row >= op.get_local_offset()) && (global_row < (op.get_local_offset() + op.get_num_local_rows()))) {
      x[global_row - op.get_local_offset()] = 1.;
      std::cout << "myrank=" << myrank << " has row" << global_row << "idx=" << global_row - op.get_local_offset() << std::endl; 
    }
    MPI_Barrier(MPI_COMM_WORLD);
    op.multiply(&x[0], &y[0]);
    MPI_Barrier(MPI_COMM_WORLD);
    for (int local_col=0; local_col<op.get_num_local_rows(); ++local_col) {
      if (y[local_col] != 0) {
	int global_col = local_col + op.get_local_offset();
	std::cout << "myrank=" << myrank << " global_row=" << global_row << " global_col=" << global_col << " value=" << y[local_col] << std::endl;
	mat.insert(global_row, 1, &global_col, &y[local_col]);
      }
    }
  }
  mat.complete();
}


class heisenberg_op : public rokko::distributed_mfree {
public:
  heisenberg_op(int L, const std::vector<std::pair<int, int> >& lattice)
    : L_(L), lattice_(lattice) {
    comm_ = MPI_COMM_WORLD;
    int size, rank;
    MPI_Comm_size(comm_, &size);
    MPI_Comm_rank(comm_, &rank);
    int n = size;
    int p = -1;
    do {
      n /= 2;
      ++p;
    } while (n > 0);
    dim_ = 1 << L;
    num_local_rows_ = 1 << (L-p);
    local_offset_ = num_local_rows_ * rank;
    buffer_.assign(num_local_rows_, 0);
  }
  ~heisenberg_op() {}

  void multiply(const double* x, double* y) const {
    rokko::heisenberg_hamiltonian::multiply(comm_, L_, lattice_, x, y, &(buffer_[0]));
  }
  int get_dim() const { return dim_; }
  int get_local_offset() const { return local_offset_; }
  int get_num_local_rows() const { return num_local_rows_; }

private:
  MPI_Comm comm_;
  int L_;
  std::vector<std::pair<int, int> > lattice_;
  int dim_, local_offset_, num_local_rows_;
  mutable std::vector<double> buffer_;
};

int main(int argc, char *argv[]) {
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::vector<std::string> solvers;
  if (argc >= 2) {
    solvers.push_back(argv[1]);
  } else {
    solvers = rokko::parallel_sparse_ev::solvers();
  }

  int L = (argc >= 3) ? boost::lexical_cast<int>(argv[2]) : 10;
  int dim = 1 << L;
  std::vector<std::pair<int, int> > lattice;
  for (int i = 0; i < L; ++i) lattice.push_back(std::make_pair(i, (i+1) % L));

  rokko::parameters params;
  params.set("Block Size", 5);
  params.set("Maximum Iterations", 500);
  params.set("Convergence Tolerance", 1.0e-8);
  params.set("num_eigenvalues", 10);
  BOOST_FOREACH(std::string const& name, solvers) {
    rokko::parallel_sparse_ev solver(name);
    heisenberg_op op(L, lattice);
    rokko::distributed_crs_matrix mat(dim, dim, solver);
    distributed_mfree_to_crs(op, mat);
    //mat2.print();
    solver.finalize();
  }
  MPI_Finalize();
}
