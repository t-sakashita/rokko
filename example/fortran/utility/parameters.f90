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
  logical :: is
  character :: c
  character(len=255) :: str_fixed
  character(len=255), allocatable :: str_alloc_fixed
  character(len=:), pointer :: str_ptr
  character(len=:), allocatable :: str
  type(string), allocatable :: keys(:)

  call rokko_parameters_construct(params)

  call rokko_parameters_set_int(params, "INTE", 3)
  call rokko_parameters_get_int(params, "INTE", i)
  print*, "T=", i

  call rokko_parameters_set_double(params, "DABURU", 1.2D0)
  call rokko_parameters_get_double(params, "DABURU", d)
  print*, "D=", d

  call rokko_parameters_set_logical(params, "is", .true.)
  call rokko_parameters_get_logical(params, "is", is)
  print*, "is=", is

  call rokko_parameters_set_char(params, "cyara", 'A')
  call rokko_parameters_get_char(params, "cyara", c)
  print*, "cyara=", c

  call rokko_parameters_set_string(params, "solver", "ansazi")
  call rokko_get(params, "solver", str)
  print*, "solver=", str

  print*, "By function:"
  str = rokko_parameters_get_string(params, "solver")
  print*, "  deferred length & allocatable=", str
  str_fixed = rokko_parameters_get_string(params, "solver")
  print*, "  fixed length=", trim(str_fixed)
  str_alloc_fixed = rokko_parameters_get_string(params, "solver")
  print*, "  fixed length & allocatable=", trim(str_alloc_fixed)
  str_ptr => rokko_parameters_get_string(params, "solver")
  print*, "  deferred length & pointer=", str_ptr

  call rokko_parameters_keys(params, keys)
  do i = 1, size(keys)
     print*, "i=", i, "str=", keys(i)%str
  enddo
  call rokko_parameters_destruct(params)

end program frank_matrix
