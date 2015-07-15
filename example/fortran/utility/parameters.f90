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
!  integer(c_int) :: i
  integer :: i
  double precision :: d
  !real(c_double) :: d
  character :: c
!  character(:) :: str
!  character(44) :: str(:)
  character(len=:), allocatable :: val

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
  
  call rokko_parameters_destruct(params)

end program frank_matrix
