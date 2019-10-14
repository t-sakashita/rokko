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

program frank_matrix
  use rokko
  implicit none
  integer :: dim
  type(rokko_eigen_matrix) :: mat
  double precision, allocatable :: array(:,:)
  double precision :: val
  integer :: i, j

  dim = 10
  print *,"dimension = ", dim
  allocate( array(dim, dim) )

  call rokko_construct(mat, array, rokko_matrix_col_major)

  ! generate frank matrix
  do i = 1, dim
     do j = 1, dim
        val = dble(dim+1 - max(i, j))
        array(i, j) = val
     enddo
  enddo
  call rokko_print(mat)

  call rokko_destruct(mat)

end program frank_matrix
