      subroutine  eigen_trd_u(
     &            a, nm,
     &            u_x, u_y, nv,
     &            u_t, v_t, i, beta, e)
!---- $     use omp_lib
!----
      use communication_s, only : bcast_dbl
     &                          , reduce_dbl
     &                          , get_loop_start
     &                          , get_loop_end
     &                          , get_owner_node
     &                          , translate_g2l
     &                          , datacast_dbl
!----
      implicit none

      integer, intent(in)    ::  nm, nv, i
      real(8), intent(inout) ::  a(1:nm)
      real(8), intent(inout) ::  u_x(1:nv), u_y(1:nv)
      real(8), intent(inout) ::  u_t(*), v_t(*)
      real(8), intent(inout) ::  beta
      real(8), intent(out)   ::  e(*)

      include 'trd.h'

      real(8)                ::  anorm2, a_n, g_n
      real(8)                ::  tt(4), ss(4), t
!---- save anorm2

      integer                ::  owner_node_col, col_local
      integer                ::  owner_node_row, row_local
      integer                ::  j_1, j_2, j_3
      integer                ::  l

      include 'param.h'
      logical, parameter     :: use_ux_buffer = .true.
!---- logical, parameter     :: use_ux_buffer = .false.


      l = i-1

      owner_node_col = get_owner_node (l, size_of_col, my_col)
      owner_node_row = get_owner_node (i, size_of_row, my_row)

      col_local   = translate_g2l(l, size_of_col, my_col)
      row_local   = translate_g2l(i, size_of_row, my_row)

      j_2         = get_loop_start(1, size_of_col, my_col)
      j_3         = get_loop_end  (l, size_of_col, my_col)

!----
!---- u=...
!----
      if ( use_ux_buffer ) then

         if ( owner_node_row == my_row ) then
            if ( beta < zero ) then
               anorm2 = zero
               do  j_1=j_2,j_3
                  t = a(j_1)
                  anorm2   =  anorm2 + t**2
                  u_x(j_1) =  t
               end do ! j_1
            else
               anorm2 = beta
            end if
               u_x(col_local+1) =  anorm2
         end if

         call bcast_dbl(u_x(1), col_local+1, owner_node_row,
     &                  mpi_comm_row)

         anorm2 = u_x(col_local+1)
         u_x(col_local+1) = zero

      else

         call bcast_dbl(a(1), col_local, owner_node_row, mpi_comm_row)

         anorm2 = zero
         do  j_1 = j_2, j_3
            anorm2   =  anorm2 + a(j_1)**2
         end do ! j_1

      end if

      if ( owner_node_col == my_col ) then
         if ( use_ux_buffer ) then
            a_n =  u_x(col_local)
         else
            a_n =  a(col_local)
         end if
      else
         a_n =  zero
      end if

      tt(1) =  anorm2
      tt(2) =  a_n
      call reduce_dbl(tt(1), ss(1), 2, 1, mpi_comm_col)
      anorm2 =  tt(1)
      a_n    =  tt(2)

      if ( anorm2 /= zero ) then

      g_n   = -sign(sqrt(anorm2), a_n)
      beta  =  anorm2 - a_n * g_n
      e (i) =  g_n

         if ( use_ux_buffer ) then
            continue
         else
            u_x(j_2:j_3) =  a(j_2:j_3)
         end if

         if ( owner_node_col == my_col ) then
            u_x(col_local) =  a_n - g_n
            if ( owner_node_row == my_row ) then
               a(col_local) =  u_x(col_local)
            end if
         else
            if ( j_3 < col_local ) then
               u_x(j_3+1:col_local) = zero
            end if
         end if

         call datacast_dbl(u_y(1), u_x(1), u_t(1), v_t(1), col_local)

      else

         beta  =  zero
         e (i) =  zero

      end if


      return
      end subroutine  eigen_trd_u

