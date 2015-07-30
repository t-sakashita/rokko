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

module rokko_dense
  use iso_c_binding
  use rokko_serial_dense
#ifdef ROKKO_HAVE_PARALLEL_DENSE_SOLVER
  use rokko_parallel_dense
#endif ROKKO_HAVE_PARALLEL_DENSE_SOLVER
  implicit none
end module rokko_dense

