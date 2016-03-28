!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/license_1_0.txt)
!
!*****************************************************************************

program solvers
  use rokko
  implicit none

  integer :: i
  type(string), allocatable :: names(:)

  call rokko_serial_dense_ev_solvers(names)
  do i = 1, size(names)
     print *, names(i)%str
  enddo

end program solvers

