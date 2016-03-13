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

module rokko_mapping_bc_mod
  use rokko_parallel_dense_classes
  use iso_c_binding
  implicit none

  interface

     subroutine rokko_mapping_bc_construct(map, global_dim, grid, solver) &
          bind(c)
       use iso_c_binding
       import rokko_grid, rokko_mapping_bc, rokko_parallel_dense_ev
       implicit none
       type(rokko_mapping_bc), intent(inout) :: map
       integer(c_int), value, intent(in) :: global_dim
       type(rokko_grid), value, intent(in) :: grid
       type(rokko_parallel_dense_ev), value, intent(in) :: solver
     end subroutine rokko_mapping_bc_construct
     
     subroutine rokko_mapping_bc_destruct(map) bind(c)
       use iso_c_binding
       import rokko_mapping_bc
       implicit none
       type(rokko_mapping_bc), intent(inout) :: map
     end subroutine rokko_mapping_bc_destruct

  end interface
  
end module rokko_mapping_bc_mod

     
module rokko_distributed_matrix_mod
  use rokko_parallel_dense_classes
  use iso_c_binding
  implicit none

  !
  ! rokko_distributed_matrix
  !

  interface
     subroutine rokko_distributed_matrix_construct(matrix, map) &
          bind(c)
       use iso_c_binding
       import rokko_mapping_bc, rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       type(rokko_mapping_bc), value, intent(in) :: map
     end subroutine rokko_distributed_matrix_construct

     subroutine rokko_distributed_matrix_construct_solver(matrix, dim1, dim2, grid, solver) &
          bind(c)
       use iso_c_binding
       import rokko_grid, rokko_parallel_dense_ev, rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(out) :: matrix
       integer(c_int), value, intent(in) :: dim1, dim2
       type(rokko_grid), value, intent(in) :: grid
       type(rokko_parallel_dense_ev), value, intent(in) :: solver
     end subroutine rokko_distributed_matrix_construct_solver
     
     subroutine rokko_distributed_matrix_destruct(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), intent(inout) :: matrix
     end subroutine rokko_distributed_matrix_destruct
     
     subroutine rokko_distributed_matrix_generate_function_c(matrix, cproc) &
          bind(c,name='rokko_distributed_matrix_generate_function')
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       type(c_funptr), value, intent(in) :: cproc
     end subroutine rokko_distributed_matrix_generate_function_c
     
     subroutine rokko_distributed_matrix_print(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end subroutine rokko_distributed_matrix_print
     
     subroutine rokko_distributed_matrix_set_local(matrix, local_i, local_j, value) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix       
       integer(c_int), value, intent(in) :: local_i, local_j
       real(c_double), value, intent(in) :: value
     end subroutine rokko_distributed_matrix_set_local
     
     function rokko_distributed_matrix_get_local(matrix, local_i,local_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       real(c_double) :: rokko_distributed_matrix_get_local
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: local_i,local_j
     end function rokko_distributed_matrix_get_local
     
     subroutine rokko_distributed_matrix_set_global(matrix, global_i, global_j, value) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: global_i, global_j
       real(c_double), value, intent(in) :: value
     end subroutine rokko_distributed_matrix_set_global
     
     function rokko_distributed_matrix_get_global(matrix, global_i, global_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       real(c_double) :: rokko_distributed_matrix_get_global
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in):: global_i, global_j
     end function rokko_distributed_matrix_get_global
     
     function rokko_distributed_matrix_get_m_local(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_m_local
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_m_local

     function rokko_distributed_matrix_get_n_local(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_n_local
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_n_local

     function rokko_distributed_matrix_get_m_global(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_m_global
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_m_global

     function rokko_distributed_matrix_get_n_global(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_n_global
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_n_global

     function rokko_distributed_matrix_get_nprocs(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_nprocs
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_nprocs

     function rokko_distributed_matrix_get_myrank(matrix) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_get_myrank
       type(rokko_distributed_matrix), value, intent(in) :: matrix
     end function rokko_distributed_matrix_get_myrank
 
     function rokko_distributed_matrix_translate_l2g_row(matrix, local_i) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_l2g_row
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int),value,intent(in) :: local_i
     end function rokko_distributed_matrix_translate_l2g_row

     function rokko_distributed_matrix_translate_l2g_col(matrix, local_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_l2g_col
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in)::local_j
     end function rokko_distributed_matrix_translate_l2g_col

     function rokko_distributed_matrix_translate_g2l_row(matrix, global_i) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_g2l_row
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in)::global_i
     end function rokko_distributed_matrix_translate_g2l_row

     function rokko_distributed_matrix_translate_g2l_col(matrix, global_j) bind(c)
       use iso_c_binding
       import rokko_distributed_matrix
       implicit none
       integer(c_int) :: rokko_distributed_matrix_translate_g2l_col
       type(rokko_distributed_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: global_j
     end function rokko_distributed_matrix_translate_g2l_col
  end interface

contains
  
  subroutine rokko_distributed_matrix_generate_array(matrix, array)
    type(rokko_distributed_matrix), intent(out) :: matrix
    double precision, intent(in) :: array(:,:)
    integer :: m_local, n_local, local_i, local_j, global_i, global_j
    m_local = rokko_distributed_matrix_get_m_local(matrix)
    n_local = rokko_distributed_matrix_get_n_local(matrix)
    do local_i = 0, m_local - 1 
       do local_j = 0, n_local - 1
          global_i = rokko_distributed_matrix_translate_l2g_row(matrix, local_i)
          global_j = rokko_distributed_matrix_translate_l2g_col(matrix, local_j)
          call rokko_distributed_matrix_set_local(matrix, local_i, local_j, &
               array(global_i + 1, global_j + 1))
       enddo
    enddo
  end subroutine rokko_distributed_matrix_generate_array

  subroutine rokko_distributed_matrix_generate_function(matrix, func_in)
    use iso_c_binding
    type(rokko_distributed_matrix), value, intent(in) :: matrix
    type(c_funptr) :: cproc
    interface
       double precision function func_in (i, j) bind(c)
         use, intrinsic :: iso_c_binding
         integer(c_int), value, intent(in) :: i, j
       end function func_in
    end interface
    ! get c procedure pointer.
    cproc = c_funloc(func_in)
    ! call wrapper written in c.
    call rokko_distributed_matrix_generate_function_c(matrix, cproc)
  end subroutine rokko_distributed_matrix_generate_function

end module rokko_distributed_matrix_mod



    
