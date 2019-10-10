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
  end interface rokko_construct

  interface rokko_destruct
     procedure rokko_eigen_vector_destruct
  end interface rokko_destruct

  interface rokko_print
     procedure rokko_eigen_vector_print
  end interface rokko_print
  
  interface
     subroutine rokko_eigen_vector_construct(vec, dim1) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       type(rokko_eigen_vector), intent(out) :: vec
       integer(c_int), value, intent(in) :: dim1
     end subroutine rokko_eigen_vector_construct

     subroutine rokko_eigen_vector_destruct(vec) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       type(rokko_eigen_vector), intent(inout) :: vec
     end subroutine rokko_eigen_vector_destruct

     subroutine rokko_eigen_vector_print(vec) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       type(rokko_eigen_vector), value, intent(in) :: vec
     end subroutine rokko_eigen_vector_print
  end interface

  interface rokko_eigen_vector_get
     function rokko_eigen_vector_get_f(vec, i) bind(c)
       use iso_c_binding
       import rokko_eigen_vector
       implicit none
       real(c_double) :: rokko_eigen_vector_get_f
       type(rokko_eigen_vector), value, intent(in) :: vec
       integer(c_int), value, intent(in) :: i
     end function rokko_eigen_vector_get_f
  end interface rokko_eigen_vector_get

end module rokko_eigen_vector_mod
