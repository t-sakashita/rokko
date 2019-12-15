/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2013-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_MPI_TYPE_TRAITS_HPP
#define ROKKO_MPI_TYPE_TRAITS_HPP

#include <mpi.h>

namespace rokko {

namespace detail {

template <typename T>
struct mpi_type_traits {
  static inline MPI_Datatype get_type();
};


#define DEF_MPI_TRAITS(type, mpi_type) \
  template<> \
  inline MPI_Datatype mpi_type_traits<type>::get_type() { \
    return mpi_type; \
  }

DEF_MPI_TRAITS(int,  MPI_INT);
DEF_MPI_TRAITS(char,  MPI_CHAR);
DEF_MPI_TRAITS(float,  MPI_FLOAT);
DEF_MPI_TRAITS(double,  MPI_DOUBLE);

#undef DEF_MPI_TRAITS

} // end namespace detail

template <typename T>
static const auto mpi_type = detail::mpi_type_traits<T>::get_type();

} // namespace rokko

#endif // ROKKO_MPI_TYPE_TRAITS_HPP
