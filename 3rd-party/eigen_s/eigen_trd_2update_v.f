      subroutine eigen_trd_2update_v(
     &           u_x, v_x,
     &           ux, vx, nv,
     &           u_t, v_t,
     &           i, i_base, m0)
!----
      use communication_s, only : reduce_dbl
     &                          , allgather_dbl
     &                          , get_loop_end
     &                          , translate_g2l
!----
      implicit none
!----
      integer, intent(in)    ::  nv, i, i_base, m0
      real(8), intent(in)    ::  u_x(*)
      real(8), intent(out)   ::  v_x(*)
      real(8), intent(in)    ::  ux(1:nv,*)
      real(8), intent(in)    ::  vx(1:nv,*)
      real(8), intent(inout) ::  u_t(*)
      real(8), intent(inout) ::  v_t(*)
!----
      real(8), parameter     ::  zero = 0d0, one = 1d0
!----
      integer                ::  j, k, l, n, ll
      integer                ::  k_1, k_2, k_3
!----
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  l_1, l_2, l_3, l_4
      integer                ::  jj_1, jj_2, jj_3, jj_4
      integer                ::  lx
!----
      integer, parameter     ::  vtol = 2048
!----
      include 'trd.h'
      include 'CSTAB.h'
!----
      real(8)                ::  w0
      real(8)                ::  u0, v0
      real(8)                ::  u1, v1
      real(8)                ::  u2, v2
      real(8)                ::  ux0, vx0
      real(8)                ::  ux1, vx1
      real(8)                ::  ux2, vx2


      k_1 = i - i_base
      k_2 = m0

      if ( k_2 <= k_1 ) return

      l = i-1
      n  = translate_g2l(l, size_of_col,my_col)

!---- for attention to unexpected overflow or nan
      j_3 = get_loop_end(l, size_of_col,my_col)
      if ( j_3 < n ) then
         v_x(j_3+1:n) = zero ! in case
      end if
!----
!---- v=v-(uv+vu)u
!----
      l_2 = k_2-k_1
      do j=1,l_2*2
         u_t(j) = zero
      end do

      l_4 = mod(k_2-k_1, 3)+k_1+1
      lx  = l1_lsize*l1_way/16

      ll = (n-1)/size_of_row+1
      ll = ((ll-1)/2+1)*2

      if ( n > vtol ) then
         jj_2 = 1+ll*(my_row-1)
         jj_3 = min(n, ll*my_row)
      else
         jj_2 = 1
         jj_3 = n
      endif

      k_3 = get_loop_end(l, size_of_col,my_col)

      do jj_1=jj_2,jj_3,lx
         j_2 = jj_1; j_3 = min(jj_1+lx-1, jj_3)

         if(l_4-1==k_1+1)then
            l_1 = k_1+1                           ! 0

            j = l_1-k_1
            u0 = u_t(2*(j+0)-1)
            v0 = u_t(2*(j+0)-0)
!dir$ ivdep
!dir$ vector always
            do j_1=j_2,j_3
               w0 = u_x(j_1)
               u0 = u0+vx(j_1,l_1+0)*w0
               v0 = v0+ux(j_1,l_1+0)*w0
            end do! j_1
            u_t(2*(j+0)-1) = u0
            u_t(2*(j+0)-0) = v0
         end if
         if(l_4-2==k_1+1)then
            l_1 = k_1+1                           ! 1

            j = l_1-k_1
            u0 = u_t(2*(j+0)-1)
            v0 = u_t(2*(j+0)-0)
            u1 = u_t(2*(j+1)-1)
            v1 = u_t(2*(j+1)-0)
!dir$ ivdep
!dir$ vector always
            do j_1=j_2,j_3
               w0 = u_x(j_1)
               u0 = u0+vx(j_1,l_1+0)*w0
               v0 = v0+ux(j_1,l_1+0)*w0
               u1 = u1+vx(j_1,l_1+1)*w0
               v1 = v1+ux(j_1,l_1+1)*w0
            end do! j_1
            u_t(2*(j+0)-1) = u0
            u_t(2*(j+0)-0) = v0
            u_t(2*(j+1)-1) = u1
            u_t(2*(j+1)-0) = v1
         end if
         do l_1=l_4,k_2,3                  ! 2

            j = l_1-k_1
            u0 = u_t(2*(j+0)-1)
            v0 = u_t(2*(j+0)-0)
            u1 = u_t(2*(j+1)-1)
            v1 = u_t(2*(j+1)-0)
            u2 = u_t(2*(j+2)-1)
            v2 = u_t(2*(j+2)-0)
!dir$ ivdep
!dir$ vector always
            do j_1=j_2,j_3
            w0 = u_x(j_1)
               u0 = u0+vx(j_1,l_1+0)*w0
               v0 = v0+ux(j_1,l_1+0)*w0
               u1 = u1+vx(j_1,l_1+1)*w0
               v1 = v1+ux(j_1,l_1+1)*w0
               u2 = u2+vx(j_1,l_1+2)*w0
               v2 = v2+ux(j_1,l_1+2)*w0
            end do! j_1
            u_t(2*(j+0)-1) = u0
            u_t(2*(j+0)-0) = v0
            u_t(2*(j+1)-1) = u1
            u_t(2*(j+1)-0) = v1
            u_t(2*(j+2)-1) = u2
            u_t(2*(j+2)-0) = v2
         end do! l_1

      end do! jj_1

      if ( n > vtol ) then
         call reduce_dbl(u_t, v_t, (k_2-k_1)*2, 1, mpi_comm_row)
      end if
      call reduce_dbl(u_t, v_t, (k_2-k_1)*2, 1, mpi_comm_col)

      if ( n > vtol ) then
         jj_2 = 1+ll*(my_row-1)
         jj_3 = min(n, ll*my_row)
      else
         jj_2 = 1
         jj_3 = n
      endif
      do jj_1=jj_2,jj_3,lx
         j_2 = jj_1; j_3 = min(jj_1+lx-1, jj_3)

         if(l_4-1==k_1+1)then
            l_1 = k_1+1                           ! 0

            j = l_1-k_1

            u0 = u_t(2*(j+0)-1)
            v0 = u_t(2*(j+0)-0)
!dir$ ivdep
!dir$ vector always
            do j_1=j_2,j_3
               w0 = v_x(j_1)
               ux0 = ux(j_1,l_1+0)
               vx0 = vx(j_1,l_1+0)
               w0 = w0
     &              -ux0*u0
     &              -vx0*v0
               v_x(j_1) = w0
            end do! j_1
         end if
         if(l_4-2==k_1+1)then
            l_1 = k_1+1                           ! 1

            j = l_1-k_1

            u0 = u_t(2*(j+0)-1)
            v0 = u_t(2*(j+0)-0)
            u1 = u_t(2*(j+1)-1)
            v1 = u_t(2*(j+1)-0)
!dir$ ivdep
!dir$ vector always
            do j_1=j_2,j_3
               w0 = v_x(j_1)
               ux0 = ux(j_1,l_1+0)
               vx0 = vx(j_1,l_1+0)
               w0 = w0
     &              -ux0*u0
     &              -vx0*v0
               ux1 = ux(j_1,l_1+1)
               vx1 = vx(j_1,l_1+1)
               w0 = w0
     &              -ux1*u1
     &              -vx1*v1
               v_x(j_1) = w0
            end do! j_1
         end if
         do l_1=l_4,k_2,3                  ! 2

            j = l_1-k_1

            u0 = u_t(2*(j+0)-1)
            v0 = u_t(2*(j+0)-0)
            u1 = u_t(2*(j+1)-1)
            v1 = u_t(2*(j+1)-0)
            u2 = u_t(2*(j+2)-1)
            v2 = u_t(2*(j+2)-0)
!dir$ ivdep
!dir$ vector always
            do j_1=j_2,j_3
               w0 = v_x(j_1)
               ux0 = ux(j_1,l_1+0)
               vx0 = vx(j_1,l_1+0)
               w0 = w0
     &              -ux0*u0
     &              -vx0*v0
               ux1 = ux(j_1,l_1+1)
               vx1 = vx(j_1,l_1+1)
               w0 = w0
     &              -ux1*u1
     &              -vx1*v1
               ux2 = ux(j_1,l_1+2)
               vx2 = vx(j_1,l_1+2)
               w0 = w0
     &              -ux2*u2
     &              -vx2*v2
               v_x(j_1) = w0
            end do! j_1
         end do! l_1

      end do! jj_1

      if ( n > vtol ) then
!----    do j_1=1,jj_2-1
!----       vx(j_1,k_1) = zero
!----    enddo
!----    do j_1=jj_3+1,n
!----       vx(j_1,k_1) = zero
!----    enddo
!----    call reduce_dbl(vx(1,k_1), v_t, n, 1, mpi_comm_row)
         call allgather_dbl(v_x(jj_2), v_t, ll, mpi_comm_row)
         j_3 = get_loop_end(l, size_of_col,my_col)
         v_x(1:j_3)=v_t(1:j_3)
      endif

!---- for attention to unexpected overflow or nan
      j_3 = get_loop_end(l, size_of_col,my_col)
      if ( j_3 < n ) then
         v_x(j_3+1:n) = zero ! in case
      end if

      return
      end subroutine ! eigen_trd_2update_v

