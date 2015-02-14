!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2013 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
!                            Synge Todo <wistaria@comp-phys.org>
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

subroutine eigen_init_wrap(comm, order)
  use eigen_libs, only : eigen_init
  implicit none
  integer, intent(in), optional :: comm
  character*(*), intent(in), optional :: order
  call eigen_init( comm, order )
  return
end subroutine eigen_init_wrap

subroutine eigen_free_wrap(flag)
  use eigen_libs, only : eigen_free
  implicit none
  integer, intent(in), optional ::  flag
  call eigen_free(flag)
  return
end subroutine eigen_free_wrap

subroutine eigen_get_matdims_wrap(NPROW, NPCOL, n, nx, ny)
  use eigen_libs, only : eigen_NB
  use eigen_devel, only : x_nnod
  implicit none
  integer, intent(in) :: NPROW, NPCOL, n
  integer, intent(out) :: nx, ny

  integer :: NB
  integer :: n1, nm, nmz, nmw, nn, larray

  n1 = ((n-1)/NPROW+1)
  call CSTAB_get_optdim( n1, 6, 16*4, 16*4*2, nm )

  NB  = eigen_NB

  nmz = ((n-1)/NPROW+1)
  nmz = ((nmz-1)/NB+1)*NB+1
  nn  = nmz
  nmz = (n-1)/NB+1
  nmz = ((nmz-1)/NPROW+1)*NB
  ! Fix on version 2.2b
  ! to avoid unexpected SIGSEGV,
  ! use the maximum of nn and nmz.
  nmz = MAX(nn, nmz)
  
  nmw = ((n-1)/NPCOL+1)
  nmw = ((nmw-1)/NB+1)*NB+1
  nn  = nmw
  nmw = (n-1)/NB+1
  nmw = ((nmw-1)/NPCOL+1)*NB
  ! Fix on version 2.2b
  ! to avoid unexpected SIGSEGV,
  ! use the maximum of nn and nmz.
  nmw = MAX(nn, nmw)
  
  larray = MAX(nmz, nm)*nmw
  
  nx = nm
  ny = (larray-1)/nm+1

  return
end subroutine eigen_get_matdims_wrap
