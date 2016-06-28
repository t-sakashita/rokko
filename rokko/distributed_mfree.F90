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

module rokko_distributed_mfree_mod
  use iso_c_binding
  implicit none

  type, bind(c) :: rokko_distributed_mfree
     type(c_ptr) :: ptr
  end type rokko_distributed_mfree

  ! generic names
  interface rokko_construct
     procedure rokko_distributed_mfree_construct
  end interface rokko_construct
  
  interface rokko_destruct
     procedure rokko_distributed_mfree_destruct
  end interface rokko_destruct
  
  interface rokko_get_num_local_rows
     procedure rokko_distributed_mfree_num_local_rows
  end interface rokko_get_num_local_rows

  interface
     subroutine rokko_distributed_mfree_f_construct(matrix, func, dim, num_local_rows) bind(c)
       use iso_c_binding
       import rokko_distributed_mfree
       implicit none
       type(rokko_distributed_mfree), intent(out) :: matrix
       type(c_funptr), intent(in), value :: func
       integer(c_int), value, intent(in) :: dim, num_local_rows
     end subroutine rokko_distributed_mfree_f_construct
    
     integer(c_int) function rokko_distributed_mfree_num_local_rows(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_mfree
       implicit none
       type(rokko_distributed_mfree), value, intent(in) :: matrix
     end function rokko_distributed_mfree_num_local_rows
     
     integer(c_int) function rokko_distributed_mfree_dim(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_mfree
       implicit none
       type(rokko_distributed_mfree), value, intent(in) :: matrix
     end function rokko_distributed_mfree_dim
     
     subroutine rokko_distributed_mfree_destruct(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_mfree
       implicit none
       type(rokko_distributed_mfree), intent(inout) :: matrix
     end subroutine rokko_distributed_mfree_destruct     
  end interface

contains

  subroutine rokko_distributed_mfree_construct(mat, multiply_in, dim, num_local_rows)
    use, intrinsic :: iso_c_binding
    type(rokko_distributed_mfree), intent(inout) :: mat
    integer(c_int), intent(in) :: dim, num_local_rows
    type(c_funptr) :: cproc
    interface
       subroutine multiply_in (n, x, y) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int), intent(in), value :: n
         real(c_double), intent(in) :: x(n)
         real(c_double), intent(out) :: y(n)
       end subroutine multiply_in
    end interface
    ! get c procedure pointer.
    cproc = c_funloc(multiply_in)
    ! call wrapper written in c.
    call rokko_distributed_mfree_f_construct(mat, cproc, dim, num_local_rows)
  end subroutine rokko_distributed_mfree_construct

end module rokko_distributed_mfree_mod

