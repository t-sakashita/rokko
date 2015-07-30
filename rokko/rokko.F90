!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

#include <rokko/config.h>

module rokko
  use iso_c_binding
  use solver_name_utility
  use parameters
  use rokko_dense
#ifdef ROKKO_HAVE_PARALLEL_SPARSE_SOLVER
  use rokko_sparse
#endif ROKKO_HAVE_PARALLEL_SPARSE_SOLVER
  implicit none
end module rokko
