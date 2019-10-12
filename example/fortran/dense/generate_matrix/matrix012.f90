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

program matrix012
  use rokko
  implicit none
  integer :: dim
  type(rokko_eigen_matrix) :: mat

  dim = 10
  print *,"dimension = ", dim

  call rokko_construct(mat, dim, dim, rokko_matrix_col_major)

  call rokko_matrix012_generate(mat)
  call rokko_print(mat)

  call rokko_destruct(mat)

end program matrix012
