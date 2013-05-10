      subroutine eigen_prd_2update_local(
     &              w, nm,
     &              ux, uy, vx, vy, nv,
     &              i_base, i)
!----
      use communication_sx, only : get_loop_start
     &                           , get_loop_end
     &                           , translate_l2g
!----
      implicit none
!----
      integer, intent(in)    ::  nm, nv
      real(8), intent(inout) ::  w(1:nm,*)
      real(8), intent(in)    ::  ux(1:nv,*),uy(1:nv,*)
      real(8), intent(in)    ::  vx(1:nv,*),vy(1:nv,*)
      integer, intent(in)    ::  i_base, i
!----
      integer                ::  k_1
      integer                ::  j, l
!----
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3
      integer                ::  l_1
      integer                ::  jj_1, jj_2, jj_3
      integer                ::  kk_1, kk_2, kk_3
      integer                ::  lx
!----
      include 'trd.h'
      include 'cstab.h'
!----
      real(8)                :: u_x0, v_x0
      real(8)                :: u_x1, v_x1
      real(8)                :: uy0_0, vy0_0
      real(8)                :: uy1_0, vy1_0
      real(8)                :: uy0_1, vy0_1
      real(8)                :: uy1_1, vy1_1
      real(8)                :: uy0_2, vy0_2
      real(8)                :: uy1_2, vy1_2
      real(8)                :: uy0_3, vy0_3
      real(8)                :: uy1_3, vy1_3
      real(8)                :: uy0_4, vy0_4
      real(8)                :: uy1_4, vy1_4
      real(8)                :: uy0_5, vy0_5
      real(8)                :: uy1_5, vy1_5
      real(8)                :: w0
      real(8)                :: w1
      real(8)                :: w2
      real(8)                :: w3
      real(8)                :: w4
      real(8)                :: w5

      integer, parameter     :: vlen  = 336

      k_1 = i - i_base
      if ( k_1 <= 1 ) return
!----
      lx  = l1_lsize*l1_way/16

      i_2 = get_loop_start(i_base+1, size_of_row,my_row)
      i_3 = get_loop_end  (i-mband,  size_of_row,my_row)
      i_4 = mod(i_3-i_2+1,6)+i_2
      if ( i_2 > i_3 ) return
!----
      l = i - mband
      jj_2 = get_loop_start(1, size_of_col,my_col)
      jj_3 = get_loop_end  (l, size_of_col,my_col)

      do jj_1=jj_2,jj_3,lx
         j_2 = jj_1; j_3 = min(jj_1+lx-1,jj_3)
         do i_1=i_2,i_4-1                         ! 0
            j   = translate_l2g(i_1, size_of_row,my_row)
            l_1 = j-i_base
            uy0_0 = uy(i_1+0,k_1-0)
            vy0_0 = vy(i_1+0,k_1-0)
            uy1_0 = uy(i_1+0,k_1-1)
            vy1_0 = vy(i_1+0,k_1-1)



!dir$ ivdep
!dir$ vector always
            do j_1=j_2,j_3
               u_x0 = ux(j_1,k_1-0)
               v_x0 = vx(j_1,k_1-0)
               u_x1 = ux(j_1,k_1-1)
               v_x1 = vx(j_1,k_1-1)

               w0 = w(j_1,l_1+0*size_of_row)
               w0 = w0
     &               - (u_x0 * vy0_0)
     &               - (v_x0 * uy0_0)
               w0 = w0
     &               - (u_x1 * vy1_0)
     &               - (v_x1 * uy1_0)
               w(j_1,l_1+0*size_of_row) = w0
            end do! j_1


         end do! l_1
         do i_1=i_4,i_3,6                     ! 5
            j   = translate_l2g(i_1, size_of_row,my_row)
            l_1 = j-i_base
            uy0_0 = uy(i_1+0,k_1-0)
            vy0_0 = vy(i_1+0,k_1-0)
            uy1_0 = uy(i_1+0,k_1-1)
            vy1_0 = vy(i_1+0,k_1-1)
            uy0_1 = uy(i_1+1,k_1-0)
            vy0_1 = vy(i_1+1,k_1-0)
            uy1_1 = uy(i_1+1,k_1-1)
            vy1_1 = vy(i_1+1,k_1-1)
            uy0_2 = uy(i_1+2,k_1-0)
            vy0_2 = vy(i_1+2,k_1-0)
            uy1_2 = uy(i_1+2,k_1-1)
            vy1_2 = vy(i_1+2,k_1-1)
            uy0_3 = uy(i_1+3,k_1-0)
            vy0_3 = vy(i_1+3,k_1-0)
            uy1_3 = uy(i_1+3,k_1-1)
            vy1_3 = vy(i_1+3,k_1-1)
            uy0_4 = uy(i_1+4,k_1-0)
            vy0_4 = vy(i_1+4,k_1-0)
            uy1_4 = uy(i_1+4,k_1-1)
            vy1_4 = vy(i_1+4,k_1-1)
            uy0_5 = uy(i_1+5,k_1-0)
            vy0_5 = vy(i_1+5,k_1-0)
            uy1_5 = uy(i_1+5,k_1-1)
            vy1_5 = vy(i_1+5,k_1-1)

            do kk_1=j_2,j_3,vlen
               kk_2=kk_1; kk_3=min(kk_1+vlen-1,j_3)


!dir$ ivdep
!dir$ vector always
               do j_1=kk_2,kk_3
                  u_x0 = ux(j_1,k_1-0)
                  v_x0 = vx(j_1,k_1-0)
                  u_x1 = ux(j_1,k_1-1)
                  v_x1 = vx(j_1,k_1-1)

                  w0 = w(j_1,l_1+0*size_of_row)
                  w1 = w(j_1,l_1+1*size_of_row)
                  w0 = w0
     &               - (u_x0 * vy0_0)
     &               - (v_x0 * uy0_0)
                  w1 = w1
     &               - (u_x0 * vy0_1)
     &               - (v_x0 * uy0_1)
                  w0 = w0
     &               - (u_x1 * vy1_0)
     &               - (v_x1 * uy1_0)
                  w1 = w1
     &               - (u_x1 * vy1_1)
     &               - (v_x1 * uy1_1)
                  w(j_1,l_1+0*size_of_row) = w0
                  w(j_1,l_1+1*size_of_row) = w1
               end do! j_1

!dir$ ivdep
!dir$ vector always
               do j_1=kk_2,kk_3
                  u_x0 = ux(j_1,k_1-0)
                  v_x0 = vx(j_1,k_1-0)
                  u_x1 = ux(j_1,k_1-1)
                  v_x1 = vx(j_1,k_1-1)

                  w2 = w(j_1,l_1+2*size_of_row)
                  w3 = w(j_1,l_1+3*size_of_row)
                  w2 = w2
     &               - (u_x0 * vy0_2)
     &               - (v_x0 * uy0_2)
                  w3 = w3
     &               - (u_x0 * vy0_3)
     &               - (v_x0 * uy0_3)
                  w2 = w2
     &               - (u_x1 * vy1_2)
     &               - (v_x1 * uy1_2)
                  w3 = w3
     &               - (u_x1 * vy1_3)
     &               - (v_x1 * uy1_3)
                  w(j_1,l_1+2*size_of_row) = w2
                  w(j_1,l_1+3*size_of_row) = w3


                  w4 = w(j_1,l_1+4*size_of_row)
                  w5 = w(j_1,l_1+5*size_of_row)
                  w4 = w4
     &               - (u_x0 * vy0_4)
     &               - (v_x0 * uy0_4)
                  w5 = w5
     &               - (u_x0 * vy0_5)
     &               - (v_x0 * uy0_5)
                  w4 = w4
     &               - (u_x1 * vy1_4)
     &               - (v_x1 * uy1_4)
                  w5 = w5
     &               - (u_x1 * vy1_5)
     &               - (v_x1 * uy1_5)
                  w(j_1,l_1+4*size_of_row) = w4
                  w(j_1,l_1+5*size_of_row) = w5
               end do! j_1

            end do! kk_1

         end do! l_1
      end do! jj_1
!----
      return
      end subroutine ! eigen_prd_2update_local

