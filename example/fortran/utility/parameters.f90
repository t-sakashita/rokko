!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/license_1_0.txt)
!
!*****************************************************************************

program frank_matrix
  use parameters
  implicit none

  type(rokko_parameters) :: params
  integer :: i
  double precision :: d
  character :: c
  character(255) :: val_fixed
  character(len=:), allocatable :: val
  type(string), allocatable :: keys(:)

  call rokko_parameters_construct(params)

  call rokko_parameters_set_int(params, "INTE", 3)
  call rokko_parameters_get_int(params, "INTE", i)
  print*, "T=", i

  call rokko_parameters_set_double(params, "DABURU", 1.2D0)
  call rokko_parameters_get_double(params, "DABURU", d)
  print*, "D=", d

  call rokko_parameters_set_char(params, "cyara", 'A')
  call rokko_parameters_get_char(params, "cyara", c)
  print*, "cyara=", c

  call rokko_parameters_set_string(params, "solver", "ansazi")
  call rokko_parameters_get_string(params, "solver", val)
  print*, "solver=", val
  val_fixed = rokko_parameters_get_string_fixed(params, "solver")
  print*, "solver=", trim(val_fixed)

  call rokko_parameters_keys(params, keys)
  do i = 1, size(keys)
     print*, "i=", i, "str=", keys(i)%str
  enddo
  call rokko_parameters_destruct(params)

end program frank_matrix
