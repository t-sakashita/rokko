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

module rokko_eigen_matrix_mod
  use iso_c_binding
  implicit none

  enum, bind(c)
     enumerator :: rokko_matrix_col_major = 3, rokko_matrix_row_major = 4
  end enum

  type, bind(c) :: rokko_eigen_matrix
     type(c_ptr) :: ptr
     integer(c_int) :: major
  end type rokko_eigen_matrix

  ! generic names
  interface rokko_construct
     procedure rokko_eigen_matrix_construct
     procedure rokko_eigen_matrix_construct_array_sizes
     module procedure rokko_eigen_matrix_construct_array
  end interface rokko_construct

  interface rokko_destruct
     procedure rokko_eigen_matrix_destruct
  end interface rokko_destruct

  interface rokko_get_elem0
     procedure rokko_eigen_matrix_get0
  end interface rokko_get_elem0

  interface rokko_get_elem
     procedure rokko_eigen_matrix_get
  end interface rokko_get_elem

  interface rokko_set_elem0
     procedure rokko_eigen_matrix_set0
  end interface rokko_set_elem0

  interface rokko_set_elem
     procedure rokko_eigen_matrix_set
  end interface rokko_set_elem

  interface rokko_get_m
     procedure rokko_eigen_matrix_get_m
  end interface rokko_get_m

  interface rokko_get_n
     procedure rokko_eigen_matrix_get_n
  end interface rokko_get_n

  interface rokko_is_row_major
     procedure rokko_eigen_matrix_is_row_major
  end interface rokko_is_row_major

  interface rokko_is_col_major
     procedure rokko_eigen_matrix_is_col_major
  end interface rokko_is_col_major

  interface rokko_generate0
     module procedure rokko_eigen_matrix_generate_function0
  end interface rokko_generate0

  interface rokko_generate
     module procedure rokko_eigen_matrix_generate_function1
     module procedure rokko_eigen_matrix_generate_from_array
  end interface rokko_generate

  interface rokko_get_array_pointer
     module procedure rokko_eigen_matrix_get_array_pointer
  end interface rokko_get_array_pointer

  interface rokko_print
     procedure rokko_eigen_matrix_print
  end interface rokko_print

  interface
     subroutine rokko_eigen_matrix_construct(matrix, dim1, dim2, matrix_major) bind(c)
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       type(rokko_eigen_matrix), intent(out) :: matrix
       integer(c_int), value, intent(in) :: dim1, dim2
       integer(c_int), value, intent(in) :: matrix_major
     end subroutine rokko_eigen_matrix_construct

     subroutine rokko_eigen_matrix_construct_array_sizes(matrix, dim1, dim2, array, matrix_major) bind(c)
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       type(rokko_eigen_matrix), intent(out) :: matrix
       integer(c_int), value, intent(in) :: dim1, dim2
       double precision, intent(in) :: array(dim1, dim2)
       integer(c_int), value, intent(in) :: matrix_major
     end subroutine rokko_eigen_matrix_construct_array_sizes

     subroutine rokko_eigen_matrix_destruct(matrix) bind(c)
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       type(rokko_eigen_matrix), intent(inout) :: matrix
     end subroutine rokko_eigen_matrix_destruct

     subroutine rokko_eigen_matrix_print(matrix) bind(c)
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       type(rokko_eigen_matrix), value, intent(in) :: matrix
     end subroutine rokko_eigen_matrix_print

     subroutine rokko_eigen_matrix_set0(matrix, i, j, value) &
          & bind(c,name="rokko_eigen_matrix_set")
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       type(rokko_eigen_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: i, j
       real(c_double), value, intent(in) :: value
     end subroutine rokko_eigen_matrix_set0

     subroutine rokko_eigen_matrix_set(matrix, i, j, value) &
          & bind(c,name="rokko_eigen_matrix_set1")
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       type(rokko_eigen_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: i, j
       real(c_double), value, intent(in) :: value
     end subroutine rokko_eigen_matrix_set

     function rokko_eigen_matrix_get0(matrix, i, j) &
          & bind(c,name="rokko_eigen_matrix_get")
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       real(c_double) :: rokko_eigen_matrix_get0
       type(rokko_eigen_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: i, j
     end function rokko_eigen_matrix_get0

     function rokko_eigen_matrix_get(matrix, i, j) &
          & bind(c,name="rokko_eigen_matrix_get1")
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       real(c_double) :: rokko_eigen_matrix_get
       type(rokko_eigen_matrix), value, intent(in) :: matrix
       integer(c_int), value, intent(in) :: i, j
     end function rokko_eigen_matrix_get

     function rokko_eigen_matrix_get_m(matrix) bind(c)
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       integer(c_int) :: rokko_eigen_matrix_get_m
       type(rokko_eigen_matrix), value, intent(in) :: matrix
     end function rokko_eigen_matrix_get_m

     function rokko_eigen_matrix_get_n(matrix) bind(c)
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       integer(c_int) :: rokko_eigen_matrix_get_n
       type(rokko_eigen_matrix), value, intent(in) :: matrix
     end function rokko_eigen_matrix_get_n

     function rokko_eigen_matrix_is_row_major(matrix) bind(c)
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       logical(c_bool) :: rokko_eigen_matrix_is_row_major
       type(rokko_eigen_matrix), value, intent(in) :: matrix
     end function rokko_eigen_matrix_is_row_major

     function rokko_eigen_matrix_is_col_major(matrix) bind(c)
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       logical(c_bool) :: rokko_eigen_matrix_is_col_major
       type(rokko_eigen_matrix), value, intent(in) :: matrix
     end function rokko_eigen_matrix_is_col_major

     subroutine rokko_eigen_matrix_generate_function0_p(matrix, cproc) &
       & bind(c,name="rokko_eigen_matrix_generate_function_p")
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       type(rokko_eigen_matrix), value, intent(in) :: matrix
       type(c_funptr), value, intent(in) :: cproc
     end subroutine rokko_eigen_matrix_generate_function0_p

     subroutine rokko_eigen_matrix_generate_function1_p(matrix, cproc) bind(c)
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       type(rokko_eigen_matrix), value, intent(in) :: matrix
       type(c_funptr), value, intent(in) :: cproc
     end subroutine rokko_eigen_matrix_generate_function1_p

     type(c_ptr) function rokko_eigen_matrix_get_array_pointer_c(matrix) &
          & bind(c,name='rokko_eigen_matrix_get_array_pointer')
       use iso_c_binding
       import rokko_eigen_matrix
       implicit none
       type(rokko_eigen_matrix), value, intent(in) :: matrix
       type(c_ptr) :: c_array_ptr
     end function rokko_eigen_matrix_get_array_pointer_c

  end interface

contains

  subroutine rokko_eigen_matrix_generate_function0(matrix, func_in)
    type(rokko_eigen_matrix), value, intent(in) :: matrix
    type(c_funptr) :: cproc
    interface
       function func_in (i, j)
         double precision :: func_in
         integer, intent(in) :: i, j
       end function func_in
    end interface
    ! get c procedure pointer.
    cproc = c_funloc(func_in)
    ! call wrapper written in c.
    call rokko_eigen_matrix_generate_function0_p(matrix, cproc)
  end subroutine rokko_eigen_matrix_generate_function0

  subroutine rokko_eigen_matrix_generate_function1(matrix, func_in)
    type(rokko_eigen_matrix), value, intent(in) :: matrix
    type(c_funptr) :: cproc
    interface
       function func_in (i, j)
         double precision :: func_in
         integer, intent(in) :: i, j
       end function func_in
    end interface
    ! get c procedure pointer.
    cproc = c_funloc(func_in)
    ! call wrapper written in c.
    call rokko_eigen_matrix_generate_function1_p(matrix, cproc)
  end subroutine rokko_eigen_matrix_generate_function1

  subroutine rokko_eigen_matrix_generate_from_array(matrix, array)
    type(rokko_eigen_matrix), value, intent(in) :: matrix
    double precision, intent(in) :: array(:,:)
    integer :: m, n, i, j
    m = rokko_eigen_matrix_get_m(matrix)
    n = rokko_eigen_matrix_get_n(matrix)
    do i = 0, m-1
       do j = 0, n-1
          call rokko_eigen_matrix_set(matrix, i, j, array(i+1, j+1))
       enddo
    enddo
  end subroutine rokko_eigen_matrix_generate_from_array

  subroutine rokko_eigen_matrix_construct_array(matrix, array, matrix_major) bind(c)
    use iso_c_binding
    implicit none
    type(rokko_eigen_matrix), intent(out) :: matrix
    double precision, intent(in) :: array(:, :)
    integer(c_int), value, intent(in) :: matrix_major
    integer :: sizes(2)

    sizes = shape(array)
    call rokko_eigen_matrix_construct_array_sizes(matrix, sizes(1), sizes(2), array, matrix_major)
  end subroutine rokko_eigen_matrix_construct_array

  subroutine rokko_eigen_matrix_get_array_pointer(matrix, f_array_ptr)
    type(rokko_eigen_matrix), value, intent(in) :: matrix
    double precision, pointer, dimension(:,:), intent(out) :: f_array_ptr
    type(c_ptr) :: c_array_ptr
    integer(c_int) :: m, n
    c_array_ptr = rokko_eigen_matrix_get_array_pointer_c(matrix)
    m = rokko_eigen_matrix_get_m(matrix)
    n = rokko_eigen_matrix_get_n(matrix)
    call c_f_pointer(c_array_ptr, f_array_ptr, (/m,n/) )
  end subroutine rokko_eigen_matrix_get_array_pointer

end module rokko_eigen_matrix_mod
