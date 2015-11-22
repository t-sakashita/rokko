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

module rokko_parallel_dense_classes
  use iso_c_binding

  !
  ! classes
  !

  type, bind(c) :: rokko_grid
     type(c_ptr) ptr
     integer(c_int) major
  end type rokko_grid

  type, bind(c) :: rokko_parallel_dense_ev
     type(c_ptr) ptr
  end type rokko_parallel_dense_ev

  type, bind(c) :: rokko_distributed_matrix
     type(c_ptr) ptr
     integer(c_int) major
  end type rokko_distributed_matrix
  
end module rokko_parallel_dense_classes
