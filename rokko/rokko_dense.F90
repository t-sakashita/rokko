!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

#include <rokko/config.h>

module rokko_dense
  use iso_c_binding
  use rokko_serial_dense_ev_mod
#ifdef ROKKO_HAVE_PARALLEL_DENSE_SOLVER
  use rokko_parallel_dense_mod
  use collective
#endif
  implicit none
end module rokko_dense

