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
#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>

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
    solvers = rokko::parallel_sparse_solver::solvers();
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
    rokko::parallel_sparse_solver solver(name);
    heisenberg_op mat(L, lattice);
    if (rank == 0)
      std::cout << "Eigenvalue decomposition of antiferromagnetic Heisenberg chain" << std::endl
                << "solver = " << name << std::endl
                << "L = " << L << std::endl
                << "dimension = " << dim << std::endl;

    rokko::parameters info = solver.diagonalize(mat, params);

    int num_conv = info.get<int>("num_conv");
    if (num_conv == 0) MPI_Abort(MPI_COMM_WORLD, -1);
    std::vector<double> eigvec;
    solver.eigenvector(0, eigvec);
    if (rank == 0) {
      std::cout << "number of converged eigenpairs = " << num_conv << std::endl;
      std::cout << "smallest eigenvalues: ";
      for (int i = 0; i < num_conv; ++i) std::cout << ' ' << solver.eigenvalue(i);
      std::cout << std::endl;
      std::cout << "smallest eigenvector: ";
      for (int j = 0; j < eigvec.size(); ++j) std::cout << eigvec[j] << ' ';
      std::cout << std::endl;
    }
  }

  MPI_Finalize();
}
