!*****************************************************************************
!
! Rokko: Integrated Interface for libraries of eigenvalue decomposition
!
! Copyright (C) 2012-2019 by Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>,
!                            Synge Todo <wistaria@comp-phys.org>
!
! Distributed under the Boost Software License, Version 1.0. (See accompanying
! file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
!
!*****************************************************************************

subroutine eigen_init_wrap0()
  use eigen_libs_mod, only : eigen_init
  implicit none
  call eigen_init()
  return
end subroutine eigen_init_wrap0

subroutine eigen_init_wrap1(comm)
  use eigen_libs_mod, only : eigen_init
  implicit none
  integer, intent(in) :: comm
  call eigen_init(comm)
  return
end subroutine eigen_init_wrap1

subroutine eigen_init_wrap2(comm, order)
  use eigen_libs_mod, only : eigen_init
  implicit none
  integer, intent(in) :: comm
  character*(*), intent(in) :: order
  call eigen_init(comm, order)
  return
end subroutine eigen_init_wrap2

subroutine eigen_free_wrap0()
  use eigen_libs_mod, only : eigen_free
  implicit none
  call eigen_free()
  return
end subroutine eigen_free_wrap0

subroutine eigen_free_wrap1(flag)
  use eigen_libs0_mod, only : eigen_free0
  implicit none
  integer, intent(in) :: flag
  call eigen_free0(flag)
  return
end subroutine eigen_free_wrap1

subroutine eigen_get_matdims_wrap(NPROW, NPCOL, n, nx, ny)
  use eigen_libs_mod, only : eigen_NB
  use CSTAB_mod, only : CSTAB_get_optdim
  implicit none
  integer, intent(in) :: NPROW, NPCOL, n
  integer, intent(out) :: nx, ny

  integer :: NB
  integer :: n1, nm, nmz, nmw, nn, larray

  if (n <= 0) then
     nx = -1; ny = -1
     return
  end if

  n1 = ((n-1)/NPROW+1)
  call CSTAB_get_optdim(n1, 6, 16*4, 16*4*2, nm)

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

subroutine eigen_get_procs_wrap(procs, x_procs, y_procs)
  use eigen_libs_mod, only : eigen_get_procs
  implicit none
  integer, intent(out) :: procs, x_procs, y_procs
  call eigen_get_procs(procs, x_procs, y_procs)
  return
end subroutine eigen_get_procs_wrap

subroutine eigen_get_id_wrap(id, x_id, y_id)
  use eigen_libs_mod, only : eigen_get_id
  implicit none
  integer, intent(out) :: id, x_id, y_id
  call eigen_get_id(id, x_id, y_id)
  return
end subroutine eigen_get_id_wrap

function eigen_loop_start_wrap(id, x_id, y_id)
  use eigen_libs_mod, only : eigen_loop_start
  implicit none
  integer :: eigen_loop_start_wrap
  integer, intent(in) :: id, x_id, y_id
  eigen_loop_start_wrap = eigen_loop_start(id, x_id, y_id)
end function eigen_loop_start_wrap

function eigen_loop_end_wrap(id, x_id, y_id)
  use eigen_libs_mod, only : eigen_loop_end
  implicit none
  integer :: eigen_loop_end_wrap
  integer, intent(in) :: id, x_id, y_id
  eigen_loop_end_wrap = eigen_loop_end(id, x_id, y_id)
end function eigen_loop_end_wrap

function eigen_translate_l2g_wrap(ictr, nnod, inod)
  use eigen_libs_mod, only : eigen_translate_l2g
  implicit none
  integer :: eigen_translate_l2g_wrap
  integer, intent(in) :: ictr, nnod, inod
  eigen_translate_l2g_wrap = eigen_translate_l2g(ictr, nnod, inod)
end function eigen_translate_l2g_wrap

function eigen_translate_g2l_wrap(ictr, nnod, inod)
  use eigen_libs_mod, only : eigen_translate_g2l
  implicit none
  integer :: eigen_translate_g2l_wrap
  integer, intent(in) :: ictr, nnod, inod
  eigen_translate_g2l_wrap = eigen_translate_g2l(ictr, nnod, inod)
end function eigen_translate_g2l_wrap
