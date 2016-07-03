!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2014 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/license_1_0.txt)
!
!*****************************************************************************

program frank_matrix
  use solver_name_utility
  implicit none

  character(len=:), allocatable :: solver_name
  character(len=:), allocatable :: library, routine

  solver_name = "lapack:dsyev"
!  solver_name = "lapack"
  call rokko_split_solver_name(solver_name, library, routine)

  print*, "library=", library
  print*, "routine=", routine

end program frank_matrix
