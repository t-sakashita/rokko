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

  integer :: i, j

  dim = 10
  print *,"dimension = ", dim

  call rokko_construct(vec, dim)

  do i = 1, dim
     do j = 1, dim
        call rokko_set_elem_f(vec, i, dble(i))
     enddo
  enddo
  call rokko_print(vec)

  call rokko_destruct(vec)

end program generate_eigen_vector
