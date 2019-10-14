!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/license_1_0.txt)
!
!*****************************************************************************

program eigen_vector_construct_by_array
  use rokko
  implicit none
  integer :: dim
  double precision, allocatable :: array(:)
  type(rokko_eigen_vector) :: vector
  integer :: i

  dim = 5
  print *,"dimension = ", dim

  allocate( array(dim) )
  call rokko_construct(vector, array)
  do i=1, dim
     array(i) = i
  enddo

  call rokko_print(vector)

  call rokko_destruct(vector)

end program eigen_vector_construct_by_array
