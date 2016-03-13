!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

subroutine rokko_all_gather(matrix, array)
  use iso_c_binding
  use rokko_parallel_dense_classes
!  use rokko_distributed_matrix_mod
  use rokko_parallel_dense, only: rokko_distributed_matrix, rokko_distributed_matrix_get_nprocs, rokko_gather
  implicit none
  type(rokko_distributed_matrix), value, intent(in) :: matrix
  double precision, intent(in), target :: array(:,:)
  integer(c_int) :: root, nprocs, ierr
  double precision, pointer :: parray
  nprocs = rokko_distributed_matrix_get_nprocs(matrix)
  parray => array(1, 1)
  do root = 0, nprocs - 1
    ierr = rokko_gather(matrix, c_loc(parray), root)
  end do
end subroutine rokko_all_gather
