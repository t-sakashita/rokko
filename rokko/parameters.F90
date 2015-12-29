!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2015 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/license_1_0.txt)
!
!*****************************************************************************

#include <rokko/config.h>

module parameters
  use iso_c_binding
  use rokko_string
  implicit none

  !
  ! classes
  !

  type, bind(c) :: rokko_parameters
     type(c_ptr) ptr
  end type rokko_parameters

  ! generic_interface
  interface rokko_parameters_get
     module procedure rokko_parameters_get_int
     module procedure rokko_parameters_get_double
     !module procedure rokko_parameters_get_char
     module procedure rokko_parameters_get_logical
     module procedure rokko_parameters_get_string
  end interface rokko_parameters_get

  interface rokko_parameters_set
     module procedure rokko_parameters_set_int
     module procedure rokko_parameters_set_double
     !module procedure rokko_parameters_set_char
     module procedure rokko_parameters_set_logical
     module procedure rokko_parameters_get_string
  end interface rokko_parameters_set

  interface
  
     subroutine rokko_parameters_construct(params) bind(c)
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), intent(in) :: params
     end subroutine rokko_parameters_construct

     subroutine rokko_parameters_destruct(params) bind(c)
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), intent(in) :: params
     end subroutine rokko_parameters_destruct

     subroutine rokko_parameters_set_int_c(params, key, val) &
          bind(c,name='rokko_parameters_set_int')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
       character(c_char) :: key(*)
       integer(c_int), value, intent(in) :: val
     end subroutine rokko_parameters_set_int_c

     subroutine rokko_parameters_set_true_c(params, key) &
          bind(c,name='rokko_parameters_set_true')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
       character(c_char) :: key(*)
     end subroutine rokko_parameters_set_true_c

     subroutine rokko_parameters_set_false_c(params, key) &
          bind(c,name='rokko_parameters_set_false')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
       character(c_char) :: key(*)
     end subroutine rokko_parameters_set_false_c

     integer(c_int) function rokko_parameters_get_key_size_c(params, key) result(n) &
          bind(c,name='rokko_parameters_get_key_size')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
       character(c_char) :: key(*)
     end function rokko_parameters_get_key_size_c
    
     subroutine rokko_parameters_set_double_c(params, key, val) &
          bind(c,name='rokko_parameters_set_double')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
       character(c_char) :: key(*)
       real(c_double), value, intent(in) :: val
     end subroutine rokko_parameters_set_double_c

     
     subroutine rokko_parameters_set_char_c(params, key, val) &
          bind(c,name='rokko_parameters_set_char')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
       character(c_char) :: key(*)
       character(c_char), value, intent(in) :: val   
     end subroutine rokko_parameters_set_char_c

     subroutine rokko_parameters_set_string_c(params, key, val) &
          bind(c,name='rokko_parameters_set_string')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
       character(c_char), intent(in) :: key(*)
       character(c_char), intent(in) :: val(*)
     end subroutine rokko_parameters_set_string_c
     
     integer(c_int) function rokko_parameters_get_int_c (params, key) &
          bind(c,name='rokko_parameters_get_int')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
       character(c_char) :: key(*)
     end function rokko_parameters_get_int_c

     integer(c_int) function rokko_parameters_get_logical_c (params, key) &
          bind(c,name='rokko_parameters_get_logicalint')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
       character(c_char) :: key(*)
     end function rokko_parameters_get_logical_c

     real(c_double) function rokko_parameters_get_double_c (params, key) &
          bind(c,name='rokko_parameters_get_double')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
       character(c_char) :: key(*)
     end function rokko_parameters_get_double_c

     character(c_char) function rokko_parameters_get_char_c (params, key) &
          bind(c,name='rokko_parameters_get_char')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
       character(c_char) :: key(*)
     end function rokko_parameters_get_char_c

     type(c_ptr) function rokko_parameters_get_string_c (params, key) &
          bind(c,name='rokko_parameters_get_string')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
       character(c_char) :: key(*)
     end function rokko_parameters_get_string_c
     
     integer(c_int) function rokko_parameters_defined_c(params, key) &
          bind(c,name='rokko_parameters_defined')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
       character(c_char) :: key(*)
     end function rokko_parameters_defined_c

     type(c_ptr) function rokko_parameters_keys_c (params) &
          bind(c,name='rokko_parameters_keys')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
     end function rokko_parameters_keys_c

     integer(c_int) function rokko_parameters_size_c (params) &
          bind(c,name='rokko_parameters_size')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(rokko_parameters), value, intent(in) :: params
     end function rokko_parameters_size_c

     type(c_ptr) function rokko_string_i_c (ptr, i) &
          bind(c,name='rokko_string_i')
       use iso_c_binding
       import rokko_parameters
       implicit none
       type(c_ptr), value, intent(in) :: ptr
       integer(c_int), value, intent(in) :: i
     end function rokko_string_i_c
     
  end interface

contains
   
  subroutine rokko_parameters_get_int (params, key, val)
    use iso_c_binding
    implicit none
    type(rokko_parameters), value, intent(in) :: params
    character(*), intent(in) :: key
    integer, intent(out) :: val
    val =  rokko_parameters_get_int_c (params, trim(key)//c_null_char)
  end subroutine rokko_parameters_get_int
  
  subroutine rokko_parameters_get_double (params, key, val)
    use iso_c_binding
    implicit none
    type(rokko_parameters), value, intent(in) :: params
    character(*), intent(in) :: key
    double precision, intent(out) :: val
    val = rokko_parameters_get_double_c (params, trim(key)//c_null_char)
  end subroutine rokko_parameters_get_double

  subroutine rokko_parameters_get_logical (params, key, val)
    use iso_c_binding
    implicit none
    type(rokko_parameters), value, intent(in) :: params
    character(*), intent(in) :: key
    logical, intent(out) :: val
    integer(c_int) :: tmp
    tmp = rokko_parameters_get_logical_c (params, trim(key)//c_null_char)
    if (tmp == 0) then
       val = .false.
    else
       val = .true.
    endif
  end subroutine rokko_parameters_get_logical
  
  subroutine rokko_parameters_get_char (params, key, val)
    use iso_c_binding
    implicit none
    type(rokko_parameters), value, intent(in) :: params
    character(*), intent(in) :: key
    character, intent(out) :: val
    val =  rokko_parameters_get_char_c (params, trim(key)//c_null_char)
  end subroutine rokko_parameters_get_char

  subroutine rokko_parameters_get_string (params, key, val)
    use iso_c_binding
    implicit none
    type(rokko_parameters), value, intent(in) :: params
    character(*), intent(in) :: key
    character(len=:), allocatable, intent(out) :: val
    type(c_ptr) :: ptr
    character, pointer, dimension(:) :: tmp_array
    character*255 :: tmp
    integer :: i
    integer(c_int) :: n
    n = rokko_parameters_get_key_size_c (params, trim(key)//c_null_char)
    ptr = rokko_parameters_get_string_c (params, trim(key)//c_null_char)
    call c_f_pointer(ptr, tmp_array, (/n/) )
    do i=1, n
       tmp(i:i) = tmp_array(i)
    enddo
    call free_c(ptr)
    val = trim(tmp(1:n))  ! automatically allocating suitable size
  end subroutine rokko_parameters_get_string

  subroutine rokko_parameters_keys (params, keys)
    use iso_c_binding
    implicit none
    type(rokko_parameters), value, intent(in) :: params
    type(string), allocatable, intent(out) :: keys(:)
    type(c_ptr) :: ptr, ptr_i
    integer :: i, size
    character(len=:), allocatable :: str
    ptr = rokko_parameters_keys_c (params)
    size = rokko_parameters_size_c (params)
    allocate(keys(size))
    do i = 1, size
       ptr_i = rokko_string_i_c (ptr, i-1)
       call rokko_get_string(ptr_i, str)
       keys(i)%str = str
    enddo
  end subroutine rokko_parameters_keys
  
  function rokko_parameters_get_string_fixed (params, key) result(val)
    use iso_c_binding
    implicit none
    type(rokko_parameters), value, intent(in) :: params
    character(*), intent(in) :: key
    character*255 :: val
    type(c_ptr) :: ptr!(*)
    character, pointer, dimension(:) :: tmp_array
    character*255 :: tmp
    integer :: i
    integer(c_int) :: n
    n = rokko_parameters_get_key_size_c (params, trim(key)//c_null_char)
    ptr = rokko_parameters_get_string_c (params, trim(key)//c_null_char)
    call c_f_pointer(ptr, tmp_array, (/n/) )
    do i=1, n
       tmp(i:i) = tmp_array(i)
    enddo
    call free_c(ptr)
    val = trim(tmp(1:n))  ! the rest of letters of val is not changed.
    !print*, "val=", val
  end function rokko_parameters_get_string_fixed
  
  subroutine rokko_parameters_set_int (params, key, val)
    use iso_c_binding
    implicit none
    type(rokko_parameters), value, intent(in) :: params
    character(*), intent(in) :: key
    integer, value, intent(in) :: val
    call rokko_parameters_set_int_c (params, trim(key)//c_null_char, val)
  end subroutine rokko_parameters_set_int
  
  subroutine rokko_parameters_set_double (params, key, val)
    use iso_c_binding
    implicit none
    type(rokko_parameters), value, intent(in) :: params
    character(*), intent(in) :: key
    double precision, value, intent(in) :: val
    call rokko_parameters_set_double_c (params, trim(key)//c_null_char, val)
  end subroutine rokko_parameters_set_double

  subroutine rokko_parameters_set_logical (params, key, val)
    use iso_c_binding
    implicit none
    type(rokko_parameters), value, intent(in) :: params
    character(*), intent(in) :: key
    logical, value, intent(in) :: val
    if (val .eqv. .true.) then
       call rokko_parameters_set_true_c (params, trim(key)//c_null_char)
    else
       call rokko_parameters_set_false_c (params, trim(key)//c_null_char)
    endif
  end subroutine rokko_parameters_set_logical
  
  subroutine rokko_parameters_set_char (params, key, val)
    use iso_c_binding
    implicit none
    type(rokko_parameters), value, intent(in) :: params
    character(*), intent(in) :: key
    character, value, intent(in) :: val
    call rokko_parameters_set_char_c (params, trim(key)//c_null_char, val)
  end subroutine rokko_parameters_set_char

  subroutine rokko_parameters_set_string (params, key, val)
    use iso_c_binding
    implicit none
    type(rokko_parameters), value, intent(in) :: params
    character(*), intent(in) :: key
    character(*), intent(in) :: val
    call rokko_parameters_set_string_c (params, trim(key)//c_null_char, val//c_null_char)
  end subroutine rokko_parameters_set_string

end module parameters

