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
#include <rokko/mpi/mpi_type_traits.hpp>

namespace rokko {

class mpi_vector {
public:
  mpi_vector(int dim, MPI_Comm comm_in = MPI_COMM_WORLD) : comm(comm_in) {
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

  template<typename T, int SIZE>
  void scatter(const Eigen::Vector<T, SIZE>& source,
               Eigen::Vector<T, SIZE>& target, int root) const {
    if (divisible)
      MPI_Scatter(source.data(), num_local_rows, rokko::mpi_type<T>,
                  target.data(), num_local_rows, rokko::mpi_type<T>,
                  root, comm);
    else
      MPI_Scatterv(source.data(), counts.data(), displacements.data(), rokko::mpi_type<T>,
                   target.data(), num_local_rows,rokko::mpi_type<T>,
                   root, comm);
  }

  template<typename T, int SIZE>
  void gather(const Eigen::Vector<T, SIZE>& source,
              Eigen::Vector<T, SIZE>& target, int root) const {
    if (divisible)
      MPI_Gather(source.data(), num_local_rows, rokko::mpi_type<T>,
                 target.data(), num_local_rows, rokko::mpi_type<T>,
                 root, comm);
    else
      MPI_Gatherv(source.data(), num_local_rows, rokko::mpi_type<T>,
                  target.data(), counts.data(), displacements.data(), rokko::mpi_type<T>,
                  root, comm);
  }

  // for distributed_vector
  template<typename T>
  void gather(const rokko::distributed_vector<T>& source, T* target,
              int root) const {
    MPI_Gatherv(source.get_storage(), source.size_local(), rokko::mpi_type<T>,
                target, counts.data(), displacements.data(), rokko::mpi_type<T>,
                root, comm);
  }

  template<typename T, int SIZE>
  void gather(const rokko::distributed_vector<T>& source, Eigen::Vector<T, SIZE>& target, int root) const {
    gather(source, target.data(), root);
  }

  template<typename T>
  void scatter(const T* source, rokko::distributed_vector<T>& target,
               int root) const {
    MPI_Scatterv(source, counts.data(), displacements.data(), rokko::mpi_type<T>,
                 target.get_storage(), target.size_local(), rokko::mpi_type<T>,
                 root, comm);
  }

  template<typename T, int SIZE>
  void scatter(const Eigen::Vector<T, SIZE>& source, rokko::distributed_vector<T>& target, int root) const {
    scatter(source.data(), target, root);
  }

private:
  MPI_Comm comm;
  int myrank, nprocs;
  int num_local_rows;
  std::vector<int> counts, displacements;
  bool divisible;
};

} // namespace rokko

#endif // ROKKO_UTILITY_MPI_VECTOR_HPP
