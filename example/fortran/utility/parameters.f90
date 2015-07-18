!/*****************************************************************************
!*
!* Rokko: Integrated Interface for libraries of eigenvalue decomposition
!*
!* Copyright (C) 2012-2014 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
!*                            Synge Todo <wistaria@comp-phys.org>,
!*                            Tsuyoshi Okubo <t-okubo@issp.u-tokyo.ac.jp>
!*
!* Distributed under the Boost Software License, Version 1.0. (See accompanying
!* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!*
!*****************************************************************************/

program frank_matrix
  use iso_c_binding
  use parameters
  implicit none

  type(rokko_parameters) :: params
  integer :: i
  double precision :: d
  character :: c
  character(255) :: val_fixed
  character(len=:), allocatable :: val
  character(len=:), allocatable, dimension(:) :: strs

  character(10) :: key

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

  call rokko_parameters_keys(params)
!  allocate(character :: strs(3))
  
!  strs(3) = "auto"
!  strs(2) = "auto"
!  print*, "strs(3) = ", trim(strs(3))
!  print*, "end"

!  strs=(/'1. This is the first element','uoiuoi'/)
  !allocate(strs, source=["1. This is the first element", &
  !     "2. This is the second element."])

!  strs = (/ character(len=len(strs)):: 'jonas', 'sssssssssssssss','q' /)
  
  call rokko_parameters_destruct(params)

end program frank_matrix
