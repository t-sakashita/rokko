/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_UTILITY_MPI_VECTOR_HPP
#define ROKKO_UTILITY_MPI_VECTOR_HPP

#include <vector>
#include <rokko/eigen3.hpp>

namespace rokko {

class mpi_vector {
public:
  mpi_vector(int dim, MPI_Comm comm_in = MPI_COMM_WORLD, int root_in = 0) : comm(comm_in), root(root_in) {
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    divisible = (dim % nprocs) == 0;
    set_num_local_rows(dim);
    set_counts_displacements(dim);
  }

  int get_num_local_rows(int dim, int rank) const {
    return (dim + nprocs - rank - 1) / nprocs;
  }

  int get_num_local_rows(int dim) const {
    return (dim + nprocs - myrank - 1) / nprocs;
  }

  void set_num_local_rows(int dim) {
    num_local_rows = get_num_local_rows(dim);
  }

  void set_counts_displacements(int dim) {
    counts.resize(nprocs);
    displacements.resize(nprocs);
    int offset = 0;
    for (int proc=0; proc<nprocs; ++proc) {
      displacements[proc] = offset;
      counts[proc] = get_num_local_rows(dim, proc);
      offset += counts[proc];
    }
  }

  int get_count() const {
    return counts[myrank];
  }

  int get_displacement() const {
    return displacements[myrank];
  }

  template<int SIZE>
  void scatter(const Eigen::Vector<double, SIZE>& source,
               Eigen::Vector<double, SIZE>& target) const {
    if (divisible)
      MPI_Scatter(source.data(), num_local_rows, MPI_DOUBLE, target.data(), num_local_rows, MPI_DOUBLE,
                  root, comm);
    else
      MPI_Scatterv(source.data(), counts.data(), displacements.data(), MPI_DOUBLE,
                   target.data(), num_local_rows, MPI_DOUBLE,
                   root, comm);
  }

  template<int SIZE>
  void gather(const Eigen::Vector<double, SIZE>& source,
              Eigen::Vector<double, SIZE>& target) const {
    if (divisible)
      MPI_Gather(source.data(), num_local_rows, MPI_DOUBLE, target.data(), num_local_rows, MPI_DOUBLE,
                 root, comm);
    else
      MPI_Gatherv(source.data(), num_local_rows, MPI_DOUBLE,
                  target.data(), counts.data(), displacements.data(), MPI_DOUBLE,
                  root, comm);
  }

  // for distributed_vector
  void gather(const rokko::distributed_vector<double>& source, double* target,
              int root) const {
    MPI_Gatherv(source.get_storage(), source.size_local(), MPI_DOUBLE,
                target, counts.data(), displacements.data(), MPI_DOUBLE,
                root, comm);
  }

  template<int SIZE>
  void gather(const rokko::distributed_vector<double>& source, Eigen::Vector<double, SIZE>& target, int root) const {
    gather(source, target.data(), root);
  }

  void scatter(const double* source, rokko::distributed_vector<double>& target,
               int root) const {
    MPI_Scatterv(source, counts.data(), displacements.data(), MPI_DOUBLE,
                 target.get_storage(), target.size_local(), MPI_DOUBLE,
                 root, comm);
  }

  template<int SIZE>
  void scatter(const Eigen::Vector<double, SIZE>& source, rokko::distributed_vector<double>& target, int root) const {
    scatter(source.data(), target, root);
  }

private:
  MPI_Comm comm;
  int myrank, nprocs;
  int num_local_rows;
  std::vector<int> counts, displacements;
  bool divisible;
  int root;
};

} // namespace rokko

#endif // ROKKO_UTILITY_MPI_VECTOR_HPP
