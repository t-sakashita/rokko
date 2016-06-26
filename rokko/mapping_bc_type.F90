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

module rokko_mapping_bc_type
  use iso_c_binding
  implicit none

  type, bind(c) :: rokko_mapping_bc
     type(c_ptr) ptr
     integer(c_int) major
  end type rokko_mapping_bc
end module rokko_mapping_bc_type

!module rokko_parallel_dense_classes
!  use iso_c_binding

  !
  ! classes (circular dependent types)
  !

!end module rokko_parallel_dense_classes
