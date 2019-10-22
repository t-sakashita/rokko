!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/license_1_0.txt)
!
!*****************************************************************************

program generate_eigen_vector
  use rokko
  implicit none
  integer :: dim
  type(rokko_eigen_vector) :: vec

  dim = 5
  print *,"dimension = ", dim

  call rokko_construct(vec, dim)

  print *, "vec (0-started index):"
  call rokko_generate(vec, func)
  call rokko_print(vec)
  print *, "vec (1-started index):"
  call rokko_generate_f(vec, func_f)
  call rokko_print(vec)

  call rokko_destruct(vec)

contains

  function func(i) bind(c)
    use iso_c_binding
    implicit none
    real(c_double) :: func
    integer(c_int), value, intent(in) :: i
    func = dble(2*(i+1))
  end function func

  function func_f(i) bind(c)
    use iso_c_binding
    implicit none
    real(c_double) :: func_f
    integer(c_int), value, intent(in) :: i
    func_f = dble(2*i)
  end function func_f

end program generate_eigen_vector
