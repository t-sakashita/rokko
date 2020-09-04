/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MPI_VECTOR_REDUCTION_HPP
#define ROKKO_MPI_VECTOR_REDUCTION_HPP

#include <vector>
#include <rokko/eigen3.hpp>
#include <rokko/mpi/mpi_type_traits.hpp>

namespace rokko {

// for a scalar variable
template<typename T>
T all_reduce(T local_sum, MPI_Comm comm) {

  T sum;
  MPI_Allreduce(&local_sum, &sum, 1, rokko::mpi_type<T>,
                MPI_SUM, comm);

  return sum;
}

} // namespace rokko

#endif // ROKKO_MPI_VECTOR_REDUCTION_HPP
