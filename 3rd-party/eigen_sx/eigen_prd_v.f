      subroutine eigen_prd_v(u_x, v_x, v_y, nv, u_t, v_t, c, i)
!----
      use communication_sx, only : reduce_dbl
     &                           , datacast_dbl
     &                           , translate_g2l
     &                           , get_loop_start
     &                           , get_loop_end
     &                           , get_loop_node
!----
      implicit none

      integer, intent(in)    :: nv, i
      real(8), intent(in)    :: u_x(1:nv,*)
      real(8), intent(inout) :: v_x(1:nv,*), v_y(1:nv,*)
      real(8), intent(out)   :: u_t(*), v_t(*)
      real(8), intent(inout) :: c(2,2)

      real(8)                ::  tt(4), ss(4)
      real(8)                ::  s11,s12,s21,s22
      real(8)                ::  c11,c12,c21,c22
      real(8)                ::  u12
      real(8)                ::  sx(2,2), tx(2,2)
      real(8)                ::  u1, u2
      real(8)                ::  v1, v2

      integer                ::  col_local, owner_node_col
      integer                ::  row_local, owner_node_row
      integer                ::  j_1, j_2, j_3
      integer                ::  l
      integer                ::  n, ll, jj_2, jj_3

!---- integer, parameter     ::  vtol = 512
      integer, parameter     ::  vtol = 51200

      include 'trd.h'
      include 'param.h'


      l = i-2

      owner_node_col = get_loop_node (l, size_of_col,my_col)
      owner_node_row = get_loop_node (i, size_of_row,my_row)

      col_local       = translate_g2l(l, size_of_col,my_col)
      row_local       = translate_g2l(i, size_of_row,my_row)

      j_2         = get_loop_start(1, size_of_col,my_col)
      j_3         = get_loop_end  (l, size_of_col,my_col)

      n = j_3-j_2+1
      ll = (n-1)/size_of_row+1
      ll = ((ll-1)/2+1)*2

      if ( n > vtol ) then
         jj_2 = j_2+(1+ll*(my_row-1))-1
         jj_3 = j_2+(min(n,ll*my_row))-1
      else
         jj_2 = j_2
         jj_3 = j_3
      end if

! s:=v^tu=u^tau
      s11 = zero
      s12 = zero
      s22 = zero
      u12 = zero
!dir$ unroll(2)
!dir$ ivdep
!dir$ vector always
      do j_1=jj_2,jj_3
         v1 = v_x(j_1,1)
         v2 = v_x(j_1,2)
         u1 = u_x(j_1,1)
         u2 = u_x(j_1,2)
         s11  = s11 + v1 * u1
         s12  = s12 + v1 * u2
         s22  = s22 + v2 * u2
!----    u12  = u12 + u1 * u2
      end do! j_1
      ss(1) = s11
      ss(2) = s12
      ss(3) = s22
!---- ss(4) = u12


      if ( n > vtol ) then
         call reduce_dbl(ss(1), tt(1), 4, 1, mpi_comm_row)
      end if
      call reduce_dbl(ss(1), tt(1), 4, 1, mpi_comm_col)

!---- u12 = ss(4)
      u12 = c(1,2)
      c(2,1) = -c(2,2)*(c(1,1)*u12)

!---- print*,"v2'=",v_x(1:l,2)

! sx:=s
      sx(1,1) = ss(1)
      sx(2,1) = ss(2)
      sx(1,2) = ss(2)
      sx(2,2) = ss(3)

! tx:=sx*c
      tx(1,1) = sx(1,1)*c(1,1) + sx(1,2)*c(2,1)
      tx(2,1) = sx(2,1)*c(1,1) + sx(2,2)*c(2,1)
      tx(1,2) =                  sx(1,2)*c(2,2)
      tx(2,2) =                  sx(2,2)*c(2,2)
            
! sx:=c^t*tx
      sx(1,1) = c(1,1)*tx(1,1) + c(2,1)*tx(2,1)
      sx(2,1) =                  c(2,2)*tx(2,1)
      sx(1,2) = c(1,1)*tx(1,2) + c(2,1)*tx(2,2)
      sx(2,2) =                  c(2,2)*tx(2,2)

      tx(1,1) = sx(1,1)
      tx(2,2) = sx(2,2)
      tx(2,1) = sx(1,2)+sx(2,1)
      tx(1,2) = zero

      s11 = -tx(1,1)/2
      s21 = -tx(2,1)/2
      s12 = -tx(1,2)/2
      s22 = -tx(2,2)/2

      c11 = c(1,1)
      c21 = c(2,1)
      c12 = zero
      c22 = c(2,2)

!---- print*,"s=[",s11,s12
!---- print*,"  [",s21,s22

! v:=vc^t
!dir$ unroll(2)
!dir$ ivdep
!dir$ vector always
      do j_1=jj_2,jj_3
         v1 = v_x(j_1,1)
         v2 = v_x(j_1,2)
         v1         =      v1 * c(1,1) + v2 * c(2,1)
         v2         =                    v2 * c(2,2)
!---- end do! j_1
!
!---- print*,"v2''=",v_x(1:l,2)
!
! v:=v-us
!---- do j_1=jj_2,jj_3
         u1 = u_x(j_1,1)
         u2 = u_x(j_1,2)
         v_x(j_1,1) = v1 + s11 * u1 + s21 * u2
         v_x(j_1,2) = v2            + s22 * u2
      end do! j_1

!---- print*,"u2'''=",u_x(1:l,2)
!---- print*,"u1'''=",u_x(1:l,1)
!---- print*,"v2'''=",v_x(1:l,2)
!---- print*,"v1'''=",v_x(1:l,1)

      if ( n > vtol ) then
         do j_1=j_2,jj_2-1
            v_x(j_1,1)=zero
            v_x(j_1,2)=zero
         end do! j_1
         do j_1=jj_3+1,j_3
            v_x(j_1,1)=zero
            v_x(j_1,2)=zero
         end do! j_1

         if ( size_of_row > 1 ) then
            do j_1=j_2,j_3
               v_t(j_1)  =v_x(j_1,1)
               v_t(j_1+n)=v_x(j_1,2)
            end do
            call reduce_dbl(v_t(1), u_t(1), 2*n, 1, mpi_comm_row)
            do j_1=j_2,j_3
               v_x(j_1,1) = v_t(j_1)
               v_x(j_1,2) = v_t(j_1+n)
            end do
         end if
      end if

      call datacast_dbl(v_y(1,1), v_x(1,1), u_t(1), v_t(1), col_local)
      call datacast_dbl(v_y(1,2), v_x(1,2), u_t(1), v_t(1), col_local)

! for attention to unexpected overflow or nan
      j_3 = get_loop_end      (l, size_of_col, my_col)
      n   = translate_g2l(l, size_of_col, my_col)
      if ( j_3 < n ) then
         v_x(j_3+1:n,1) = zero ! in case
         v_x(j_3+1:n,2) = zero ! in case
      end if

      j_3 = get_loop_end      (l, size_of_row, my_row)
      n   = translate_g2l(l, size_of_row, my_row)
      if ( j_3 < n ) then
         v_y(j_3+1:n,1) = zero ! in case
         v_y(j_3+1:n,2) = zero ! in case
      end if


      return
      end subroutine eigen_prd_v

