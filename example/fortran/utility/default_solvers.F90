!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2020 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/license_1_0.txt)
!
!*****************************************************************************

#include <rokko/config.h>

program default_solver
  use rokko
  implicit none

  print *, rokko_serial_dense_ev_default_solver()

#ifdef ROKKO_HAVE_PARALLEL_DENSE_SOLVER
  print *, rokko_parallel_dense_ev_default_solver()
#endif

#ifdef ROKKO_HAVE_PARALLEL_SPARSE_SOLVER
  print *, rokko_parallel_sparse_ev_default_solver()
#endif

end program default_solver
