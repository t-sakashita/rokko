      subroutine eigen_trd_d_preserve(a, w, nm,
     &           d_t,
     &           u_x, u_y, v_x, v_y, nv,
     &           m0, i_base, i_block)
!----
      use communication_s, only : get_loop_start
     &                          , get_loop_end
     &                          , translate_l2g
     &                          , translate_g2l
!----
      implicit none

      integer, intent(in)    ::  nm, nv, m0, i_base, i_block
      real(8), intent(inout) ::  a(1:nm,*)
      real(8), intent(out)   ::  w(1:nm,*)
      real(8), intent(out)   ::  u_x(1:nv,*), u_y(1:nv,*)
      real(8), intent(out)   ::  v_x(1:nv,*), v_y(1:nv,*)
      real(8), intent(out)   ::  d_t(*)

      real(8), parameter     ::  zero = 0d0, one = 1d0

      integer                ::  col_local
      integer                ::  row_local
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3
      integer                ::  k_1, k_2, k_3
      integer                ::  i, j, l

      include 'trd.h'


      i_2 = get_loop_start(i_base+1,  size_of_row,my_row)
      i_3 = get_loop_end  (i_base+m0, size_of_row,my_row)
      j_2 = get_loop_start(1,         size_of_col,my_col)
      j_3 = get_loop_end  (i_base+m0, size_of_col,my_col)

      do i_1=i_2,i_3
         j = translate_l2g(i_1, size_of_row,my_row)
         do j_1=j_2,j_3
            w(j_1,j-i_base) = a(j_1,i_1)
         end do! j_1
      end do! i_1

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = get_loop_end  (i_base+m0-1, size_of_row,my_row)
         i_4 = size_of_col/n_common
         if ( size_of_row == size_of_col ) then
!dir$ ivdep
!dir$ vector always
            do i_1=i_2,i_3,i_4
               d_t(i_1)   = a(i_1,i_1)
               a(i_1,i_1) = zero
            end do! i_1
         else
            do i_1=i_2,i_3,i_4
               j   = translate_l2g(i_1, size_of_row,my_row)
               j_1 = translate_g2l(j,   size_of_col,my_col)
               d_t(i_1)   = a(j_1,i_1)
               a(j_1,i_1) = zero
            end do! i_1
         end if
      end if

      i = i_base+m0
      l = i - 2
      row_local  = translate_g2l(i, size_of_row,my_row)
      col_local  = translate_g2l(i, size_of_col,my_col)

      k_2 = m0
      k_3 = max(1, 3*(2-i_block))
      do k_1=1,k_2
         do j_1=1,col_local
            u_x(j_1,k_1) = zero
            v_x(j_1,k_1) = zero
         end do! j_1
         do j_1=1,row_local
            u_y(j_1,k_1) = zero
            v_y(j_1,k_1) = zero
         end do! j_1
      end do! k_1

      return
      end subroutine eigen_trd_d_preserve


      subroutine eigen_trd_d_restore(a, w, nm,
     &           d_t,
     &           m0, i_base)
!----
      use communication_s, only : get_loop_start
     &                          , get_loop_end
     &                          , translate_l2g
!----
      integer, intent(in)    ::  nm, m0, i_base
      real(8), intent(out)   ::  a(1:nm,*)
      real(8), intent(in)    ::  w(1:nm,*), d_t(*)

      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4

      include 'trd.h'


      i_2 = get_loop_start(i_base+1,  size_of_row,my_row)
      i_3 = get_loop_end  (i_base+m0, size_of_row,my_row)
      do i_1=i_2,i_3
         j   = translate_l2g(i_1, size_of_row,my_row)
         j_3 = get_loop_end  (j, size_of_col,my_col)
         do j_1=1,j_3
            a(j_1,i_1) = w(j_1,j-i_base)
         end do! j_1
      end do! i_1

      if ( diag_0 > 0 ) then
         i_2 = diag_0
         i_3 = get_loop_end  (i_base, size_of_row,my_row)
         i_4 = size_of_col/n_common
         if ( size_of_row == size_of_col ) then
            do i_1=i_2,i_3,i_4
               j_1 = i_1
               a(j_1,i_1) = d_t(i_1)
            end do! i_1
         else
            do i_1=i_2,i_3,i_4
!----              j   = translate_l2g(i_1, size_of_row,my_row)
!----              j_1 = translate_g2l(j,   size_of_col,my_col)
               j   = (i_1-1)*size_of_row+my_row
               j_1 = (j-1)/size_of_col+1
               a(j_1,i_1) = d_t(i_1)
            end do! i_1
         end if
      end if

      return
      end subroutine eigen_trd_d_restore

