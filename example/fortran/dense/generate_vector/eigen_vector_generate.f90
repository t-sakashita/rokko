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
  call rokko_generate0(vec, func0)
  call rokko_print(vec)
  print *, "vec (1-started index):"
  call rokko_generate(vec, func)
  call rokko_print(vec)

  call rokko_destruct(vec)

contains

  function func0(i)
    implicit none
    double precision :: func0
    integer, intent(in) :: i
    func0 = dble(2*(i+1))
  end function func0

  function func(i)
    implicit none
    double precision :: func
    integer, intent(in) :: i
    func = dble(2*i)
  end function func

end program generate_eigen_vector
