!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2016 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

module collective
  use iso_c_binding
  use rokko_distributed_matrix_mod
  implicit none

  interface
     function rokko_gather(matrix, array, root) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_gather
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       type(c_ptr), value, intent(in) :: array
       integer(c_int), value ::root
     end function rokko_gather

     function rokko_scatter(array, matrix, root) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_scatter
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       type(c_ptr), value, intent(in) :: array
       integer(c_int), value :: root
     end function rokko_scatter
  end interface

contains

  subroutine rokko_all_gather(matrix, array)
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

end module collective
