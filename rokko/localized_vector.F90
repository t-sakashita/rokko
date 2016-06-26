module localized_vector_mod
  use iso_c_binding
  implicit none

  type, bind(c) :: rokko_localized_vector
     type(c_ptr) ptr
  end type rokko_localized_vector

  ! generic names
  interface construct
     procedure rokko_localized_vector_construct
  end interface construct

  interface destruct
     procedure rokko_localized_vector_destruct
  end interface destruct

!  interface print
!     procedure rokko_localized_vector_print
!  end interface print
  
  interface
     subroutine rokko_localized_vector_construct(vec, dim1) bind(c)
       use iso_c_binding
       import rokko_localized_vector
       implicit none
       type(rokko_localized_vector), intent(out) :: vec
       integer(c_int), value, intent(in) :: dim1
     end subroutine rokko_localized_vector_construct

     subroutine rokko_localized_vector_destruct(vec) bind(c)
       use iso_c_binding
       import rokko_localized_vector
       implicit none
       type(rokko_localized_vector), intent(inout) :: vec
     end subroutine rokko_localized_vector_destruct
  end interface

  interface rokko_localized_vector_get
     function rokko_localized_vector_get_f(vec, i) bind(c)
       use iso_c_binding
       import rokko_localized_vector
       implicit none
       real(c_double) :: rokko_localized_vector_get_f
       type(rokko_localized_vector), value, intent(in) :: vec
       integer(c_int), value, intent(in) :: i
     end function rokko_localized_vector_get_f
  end interface rokko_localized_vector_get

end module localized_vector_mod
