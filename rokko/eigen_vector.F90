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

module rokko_eigen_vector_mod
  use iso_c_binding
  implicit none

  type, bind(c) :: rokko_eigen_vector
     type(c_ptr) :: ptr
  end type rokko_eigen_vector

  ! generic names
  interface rokko_construct
     procedure rokko_eigen_vector_construct
     procedure rokko_eigen_vector_construct_array_size
     module procedure rokko_eigen_vector_construct_array
  end interface rokko_construct

  interface rokko_destruct
     procedure rokko_eigen_vector_destruct
  end interface rokko_destruct

  interface rokko_get_elem
     procedure rokko_eigen_vector_get
  end interface rokko_get_elem

  interface rokko_get_elem_f
     procedure rokko_eigen_vector_get_f
  end interface rokko_get_elem_f

  interface rokko_set_elem
     procedure rokko_eigen_vector_set
  end interface rokko_set_elem

  interface rokko_set_elem_f
     procedure rokko_eigen_vector_set_f
  end interface rokko_set_elem_f

  interface rokko_print
     procedure rokko_eigen_vector_print
  end interface rokko_print

  interface rokko_get_array_pointer
     module procedure rokko_eigen_vector_get_array_pointer
  end interface rokko_get_array_pointer

  interface rokko_generate
     module procedure rokko_eigen_vector_generate_function
  end interface rokko_generate

  interface rokko_generate_f
     module procedure rokko_eigen_vector_generate_function_f
  end interface rokko_generate_f

  interface
     subroutine rokko_eigen_vector_construct(vec, dim) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       type(rokko_eigen_vector), intent(out) :: vec
       integer(c_int), value, intent(in) :: dim
     end subroutine rokko_eigen_vector_construct

     subroutine rokko_eigen_vector_construct_array_size(vector, dim, array) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       type(rokko_eigen_vector), intent(out) :: vector
       integer(c_int), value, intent(in) :: dim
       double precision, intent(in) :: array(dim)
     end subroutine rokko_eigen_vector_construct_array_size

     subroutine rokko_eigen_vector_destruct(vec) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       type(rokko_eigen_vector), intent(inout) :: vec
     end subroutine rokko_eigen_vector_destruct

     function rokko_eigen_vector_get_dim(vector) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       integer(c_int) :: rokko_eigen_vector_get_dim
       type(rokko_eigen_vector), value, intent(in) :: vector
     end function rokko_eigen_vector_get_dim

     subroutine rokko_eigen_vector_print(vec) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       type(rokko_eigen_vector), value, intent(in) :: vec
     end subroutine rokko_eigen_vector_print

     type(c_ptr) function rokko_eigen_vector_get_array_pointer_c(vector) &
          & bind(c,name='rokko_eigen_vector_get_array_pointer')
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       type(rokko_eigen_vector), value, intent(in) :: vector
       type(c_ptr) :: c_array_ptr
     end function rokko_eigen_vector_get_array_pointer_c

     subroutine rokko_eigen_vector_set(vec, i, value) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       type(rokko_eigen_vector), value, intent(in) :: vec
       integer(c_int), value, intent(in) :: i
       real(c_double), value, intent(in) :: value
     end subroutine rokko_eigen_vector_set

     subroutine rokko_eigen_vector_set_f(vec, i, value) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       type(rokko_eigen_vector), value, intent(in) :: vec
       integer(c_int), value, intent(in) :: i
       real(c_double), value, intent(in) :: value
     end subroutine rokko_eigen_vector_set_f

     function rokko_eigen_vector_get(vec, i) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       real(c_double) :: rokko_eigen_vector_get
       type(rokko_eigen_vector), value, intent(in) :: vec
       integer(c_int), value, intent(in) :: i
     end function rokko_eigen_vector_get

     function rokko_eigen_vector_get_f(vec, i) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       real(c_double) :: rokko_eigen_vector_get_f
       type(rokko_eigen_vector), value, intent(in) :: vec
       integer(c_int), value, intent(in) :: i
     end function rokko_eigen_vector_get_f

     subroutine rokko_eigen_vector_generate_function_p(vec, cproc) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       type(rokko_eigen_vector), value, intent(in) :: vec
       type(c_funptr), value, intent(in) :: cproc
     end subroutine rokko_eigen_vector_generate_function_p

     subroutine rokko_eigen_vector_generate_function_f_p(vec, cproc) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       type(rokko_eigen_vector), value, intent(in) :: vec
       type(c_funptr), value, intent(in) :: cproc
     end subroutine rokko_eigen_vector_generate_function_f_p

  end interface

contains

  subroutine rokko_eigen_vector_generate_function(vector, func_in)
    type(rokko_eigen_vector), value, intent(in) :: vector
    type(c_funptr) :: cproc
    interface
       function func_in(i)
         double precision :: func_in
         integer, intent(in) :: i
       end function func_in
    end interface
    ! get c procedure pointer.
    cproc = c_funloc(func_in)
    ! call wrapper written in c.
    call rokko_eigen_vector_generate_function_p(vector, cproc)
  end subroutine rokko_eigen_vector_generate_function

  subroutine rokko_eigen_vector_generate_function_f(vector, func_in)
    type(rokko_eigen_vector), value, intent(in) :: vector
    type(c_funptr) :: cproc
    interface
       function func_in(i)
         double precision :: func_in
         integer, intent(in) :: i
       end function func_in
    end interface
    ! get c procedure pointer.
    cproc = c_funloc(func_in)
    ! call wrapper written in c.
    call rokko_eigen_vector_generate_function_f_p(vector, cproc)
  end subroutine rokko_eigen_vector_generate_function_f

  subroutine rokko_eigen_vector_construct_array(vector, array) bind(c)
    use iso_c_binding
    implicit none
    type(rokko_eigen_vector), intent(out) :: vector
    double precision, intent(in) :: array(:)

    call rokko_eigen_vector_construct_array_size(vector, size(array), array)
  end subroutine rokko_eigen_vector_construct_array

  subroutine rokko_eigen_vector_get_array_pointer(vector, f_array_ptr)
    type(rokko_eigen_vector), value, intent(in) :: vector
    double precision, pointer, dimension(:), intent(out) :: f_array_ptr
    type(c_ptr) :: c_array_ptr
    integer(c_int) :: dim
    c_array_ptr = rokko_eigen_vector_get_array_pointer_c(vector)
    dim = rokko_eigen_vector_get_dim(vector)
    call c_f_pointer(c_array_ptr, f_array_ptr, (/dim/) )
  end subroutine rokko_eigen_vector_get_array_pointer

end module rokko_eigen_vector_mod
