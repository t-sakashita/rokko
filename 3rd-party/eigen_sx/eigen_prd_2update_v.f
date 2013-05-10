      subroutine eigen_prd_2update_v(
     &              ux, vx, nv,
     &              u_t, v_t,
     &              i, i_base, m0)
!----
      use communication_sx, only : reduce_dbl
     &                           , translate_g2l
     &                           , get_loop_end
!----
      implicit none
!
      integer, intent(in)    ::  nv, i, i_base, m0
      real(8), intent(inout) ::  ux(1:nv,*)
      real(8), intent(inout) ::  vx(1:nv,*)
      real(8), intent(inout) ::  u_t(4,*)
      real(8), intent(inout) ::  v_t(*)
!----
      integer                ::  j, l, n, ll
      integer                ::  k_1, k_2
!----
      integer                ::  j_1, j_2, j_3
      integer                ::  l_1, l_2, l_4
      integer                ::  jj_1, jj_2, jj_3
      integer                ::  lx
!----
      include 'trd.h'
      include 'param.h'
      include 'cstab.h'
!----
      real(8)                ::  w0, w1
      real(8)                ::  u0_0, v0_0
      real(8)                ::  u1_0, v1_0
      real(8)                ::  u0_1, v0_1
      real(8)                ::  u1_1, v1_1
      real(8)                ::  ux0, vx0
      real(8)                ::  ux1, vx1


      k_1 = i - i_base
      k_2 = m0

      if ( k_2 <= k_1 ) return

      l = i - mband
      n  = translate_g2l(l, size_of_col,my_col)

! for attention to unexpected overflow or nan
      j_3 = get_loop_end(l, size_of_col,my_col)
      if ( j_3 < n ) then
         vx(j_3+1:n,k_1-0) = zero ! in case
         vx(j_3+1:n,k_1-1) = zero ! in case
      end if
!----
! v=v-(uv+vu)u
!----
      l_2 = k_2-k_1
      do j=1,l_2*(2*mband)
         u_t(j,1) = zero
      end do

      l_4 = mod(k_2-k_1,2)+k_1+1
      lx  = l1_lsize*l1_way/16

      ll = (n-1)/size_of_row+1
      ll = ((ll-1)/2+1)*2

      jj_2 = 1+ll*(my_row-1)
      jj_3 = min(n,ll*my_row)
      do jj_1=jj_2,jj_3,lx
         j_2 = jj_1; j_3 = min(jj_1+lx-1,jj_3)

         do l_1=k_1+1,l_4-1                        ! 0

            j = l_1-k_1
            u0_0 = u_t(1,j+0)
            v0_0 = u_t(2,j+0)
            u1_0 = u_t(3,j+0)
            v1_0 = u_t(4,j+0)
!dir$ ivdep
!dir$ vector always
            do j_1=j_2,j_3
               w0 = ux(j_1,k_1-0)
               w1 = ux(j_1,k_1-1)
               u0_0 = u0_0 + vx(j_1,l_1+0) * w0
               v0_0 = v0_0 + ux(j_1,l_1+0) * w0
               u1_0 = u1_0 + vx(j_1,l_1+0) * w1
               v1_0 = v1_0 + ux(j_1,l_1+0) * w1
            end do! j_1
            u_t(1,j+0) = u0_0
            u_t(2,j+0) = v0_0
            u_t(3,j+0) = u1_0
            u_t(4,j+0) = v1_0
         end do! l_1
         do l_1=l_4,k_2,2                      ! 1

            j = l_1-k_1
            u0_0 = u_t(1,j+0)
            v0_0 = u_t(2,j+0)
            u1_0 = u_t(3,j+0)
            v1_0 = u_t(4,j+0)
            u0_1 = u_t(1,j+1)
            v0_1 = u_t(2,j+1)
            u1_1 = u_t(3,j+1)
            v1_1 = u_t(4,j+1)
!dir$ ivdep
!dir$ vector always
            do j_1=j_2,j_3
               w0 = ux(j_1,k_1-0)
               w1 = ux(j_1,k_1-1)
               u0_0 = u0_0 + vx(j_1,l_1+0) * w0
               v0_0 = v0_0 + ux(j_1,l_1+0) * w0
               u1_0 = u1_0 + vx(j_1,l_1+0) * w1
               v1_0 = v1_0 + ux(j_1,l_1+0) * w1
               u0_1 = u0_1 + vx(j_1,l_1+1) * w0
               v0_1 = v0_1 + ux(j_1,l_1+1) * w0
               u1_1 = u1_1 + vx(j_1,l_1+1) * w1
               v1_1 = v1_1 + ux(j_1,l_1+1) * w1
            end do! j_1
            u_t(1,j+0) = u0_0
            u_t(2,j+0) = v0_0
            u_t(3,j+0) = u1_0
            u_t(4,j+0) = v1_0
            u_t(1,j+1) = u0_1
            u_t(2,j+1) = v0_1
            u_t(3,j+1) = u1_1
            u_t(4,j+1) = v1_1
         end do! l_1

      end do! jj_1

      call reduce_dbl(u_t, v_t, (k_2-k_1)*(2*mband), 1, mpi_comm_row)
      call reduce_dbl(u_t, v_t, (k_2-k_1)*(2*mband), 1, mpi_comm_col)

      if ( n > 100000 ) then
         jj_2 = 1+ll*(my_row-1)
         jj_3 = min(n,ll*my_row)
      else
         jj_2 = 1
         jj_3 = n
      end if
      do jj_1=jj_2,jj_3,lx
         j_2 = jj_1; j_3 = min(jj_1+lx-1,jj_3)

         do l_1=k_1+1,l_4-1                        ! 0

            j = l_1-k_1

            u0_0 = u_t(1,j+0)
            v0_0 = u_t(2,j+0)
            u1_0 = u_t(3,j+0)
            v1_0 = u_t(4,j+0)
!dir$ ivdep
!dir$ vector always
            do j_1=j_2,j_3
               w0 = vx(j_1,k_1-0)
               w1 = vx(j_1,k_1-1)
               ux0 = ux(j_1,l_1+0)
               vx0 = vx(j_1,l_1+0)
               w0 = w0
     &            - ux0 * u0_0
     &            - vx0 * v0_0
               w1 = w1
     &            - ux0 * u1_0
     &            - vx0 * v1_0
               vx(j_1,k_1-0) = w0
               vx(j_1,k_1-1) = w1
            end do! j_1
         end do! l_1
         do l_1=l_4,k_2,2                      ! 1

            j = l_1-k_1

            u0_0 = u_t(1,j+0)
            v0_0 = u_t(2,j+0)
            u1_0 = u_t(3,j+0)
            v1_0 = u_t(4,j+0)
            u0_1 = u_t(1,j+1)
            v0_1 = u_t(2,j+1)
            u1_1 = u_t(3,j+1)
            v1_1 = u_t(4,j+1)
!dir$ ivdep
!dir$ vector always
            do j_1=j_2,j_3
               w0 = vx(j_1,k_1-0)
               w1 = vx(j_1,k_1-1)
               ux0 = ux(j_1,l_1+0)
               vx0 = vx(j_1,l_1+0)
               w0 = w0
     &            - ux0 * u0_0
     &            - vx0 * v0_0
               w1 = w1
     &            - ux0 * u1_0
     &            - vx0 * v1_0
               ux1 = ux(j_1,l_1+1)
               vx1 = vx(j_1,l_1+1)
               w0 = w0
     &            - ux1 * u0_1
     &            - vx1 * v0_1
               w1 = w1
     &            - ux1 * u1_1
     &            - vx1 * v1_1
               vx(j_1,k_1-0) = w0
               vx(j_1,k_1-1) = w1
            end do! j_1
         end do! l_1

      end do! jj_1

      if ( n > 100000 ) then
         do j_1=1,jj_2-1
            vx(j_1,k_1-0) = zero
            vx(j_1,k_1-1) = zero
         end do! j_1
         do j_1=jj_3+1,n
            vx(j_1,k_1-0) = zero
            vx(j_1,k_1-1) = zero
         end do! j_1

         do j_1=1,n
            v_t(j_1)   = vx(1,k_1-0)
            v_t(j_1+n) = vx(1,k_1-1)
         end do! j_1
         call reduce_dbl(v_t, u_t, mband*n, 1, mpi_comm_row)
         do j_1=1,n
            vx(1,k_1-0) = v_t(j_1)
            vx(1,k_1-1) = v_t(j_1+n)
         end do! j_1
      end if

!  for attention to unexpected overflow or nan
      j_3 = get_loop_end(l, size_of_col,my_col)
      if ( j_3 < n ) then
         vx(j_3+1:n,k_1-0) = zero ! in case
         vx(j_3+1:n,k_1-1) = zero ! in case
      end if


      return
      end subroutine ! eigen_prd_2update_v

