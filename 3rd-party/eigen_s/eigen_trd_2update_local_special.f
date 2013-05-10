      subroutine eigen_trd_2update_local_special(
     &           w, nm,
     &           ux, uy, vx, vy, nv, uxx, anorm2,
     &           i_base, i)
!----
      use communication_s, only : get_loop_start
     &                          , get_loop_end
!----
      implicit none
!----
      integer, intent(in)    ::  nm, nv
      real(8), intent(inout) ::  w(1:nm)
      real(8), intent(in)    ::  ux(1:nv),uy(1:nv)
      real(8), intent(in)    ::  vx(1:nv),vy(1:nv)
      real(8), intent(out)   ::  uxx(1:nv),anorm2
      integer, intent(in)    ::  i_base, i
!----
      real(8), parameter     ::  zero = 0d0, one = 1d0
!----
      integer                ::  k_1
      integer                ::  l
!----
      integer                ::  i_1, i_2, i_3
      integer                ::  j_1, j_2, j_3
!----
      include 'mpif.h'
      include 'trd.h'
      include 'CSTAB.h'
!----
      real(8)                :: u_x, v_x
      real(8)                :: uy0, vy0
      real(8)                :: w0

      k_1 = i - i_base
      if ( k_1 <= 1 ) return
!----
      l = i - 1
      i_2 = get_loop_start(l, size_of_row,my_row)
      i_3 = get_loop_end  (l, size_of_row,my_row)
      if ( i_2 > i_3 ) return

      j_2 = get_loop_start(1, size_of_col,my_col)
      j_3 = get_loop_end  (l-1, size_of_col,my_col)

      i_1=i_2

      uy0 = uy(i_1+0)
      vy0 = vy(i_1+0)

      anorm2=zero
!dir$ ivdep
!dir$ vector always
      do j_1=j_2,j_3
         u_x = ux(j_1)
         v_x = vx(j_1)
         w0 = w(j_1)
     &        -(u_x*vy0)
     &        -(v_x*uy0)
         anorm2 = anorm2 + w0**2
         w(j_1) = w0
         uxx(j_1) = w0
      end do! j_1

      j_2 = get_loop_start(l, size_of_col,my_col)
      j_3 = get_loop_end  (l, size_of_col,my_col)
!dir$ ivdep
!dir$ vector always
      do j_1=j_2,j_3
         u_x = ux(j_1)
         v_x = vx(j_1)
         w0 = w(j_1)
     &        -(u_x*vy0)
     &        -(v_x*uy0)
         w(j_1) = w0
      end do! j_1

      return
      end subroutine ! eigen_trd_2update_local_special

