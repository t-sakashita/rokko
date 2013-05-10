      subroutine eigen_trd_v(u_x, v_x, v_y, nv, u_t, v_t, beta, i)
!----
      use communication_s, only : reduce_dbl
     &                          , allgather_dbl
     &                          , get_loop_start
     &                          , get_loop_end
     &                          , get_owner_node
     &                          , translate_g2l
     &                          , datacast_dbl
!----
      implicit none

      integer, intent(in)    ::  nv, i
      real(8), intent(in)    ::  u_x(1:nv)
      real(8), intent(inout) ::  v_x(1:nv), v_y(1:nv)
      real(8), intent(out)   ::  u_t(*), v_t(*)
      real(8), intent(in)    ::  beta

!---- real(8)                ::  alpha, temp
      real(8)                ::  alpha(1:1), temp(1:1)
!---- save                       alpha

      real(8), parameter     ::  zero = 0d0, one = 1d0

      integer                ::  col_local
      integer                ::  row_local
      integer                ::  j_1, j_2, j_3
      integer                ::  l
      integer                ::  n, ll, jj_1, jj_2, jj_3

!---- integer, parameter     ::  vtol = 512
      integer, parameter     ::  vtol = 51200

      include 'trd.h'


      l = i-1

      col_local   = translate_g2l(l, size_of_col,my_col)
      row_local   = translate_g2l(i, size_of_row,my_row)

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
      endif
!----
!---- v':= v-((u,v)/2|u|^2)u
!----
      if ( beta /= zero ) then

         alpha(1) = zero
!---- dir$ ivdep
         do j_1=jj_2,jj_3
            alpha(1)    = alpha(1)+v_x(j_1)*u_x(j_1)
         end do! j_1

         if ( n > vtol ) then
            call reduce_dbl(alpha, temp, 1, 1, mpi_comm_row)
         end if
         call reduce_dbl(alpha, temp, 1, 1, mpi_comm_col)

         alpha(1) = alpha(1)/(2*beta)

!---- dir$ ivdep
         do j_1=jj_2,jj_3
            v_x(j_1) = (v_x(j_1)-alpha(1)*u_x(j_1))/beta
         end do! j_1

         if ( n > vtol ) then
            call allgather_dbl(v_x(jj_2), v_t(j_2), ll, mpi_comm_row)
            v_x(j_2:j_3) = v_t(j_2:j_3)
         end if

         call datacast_dbl(v_y(1), v_x(1), u_t(1), v_t(1), col_local)

      else

         col_local   = translate_g2l(i, size_of_col,my_col)
         row_local   = translate_g2l(i, size_of_row,my_row)

         v_x(1:col_local) = zero
         v_y(1:row_local) = zero

      end if

!---- for attention to unexpected overflow or nan
      j_3 = get_loop_end      (l, size_of_col, my_col)
      n   = translate_g2l(l, size_of_col, my_col)
      if ( j_3 < n ) then
         v_x(j_3+1:n) = zero ! in case
      end if


      return
      end subroutine eigen_trd_v

