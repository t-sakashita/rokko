      subroutine eigen_prd_u(
     &              a, nm,
     &              u_x, u_y, nv,
     &              u_t, v_t, i, c, e, ne)
!----
      use communication_sx, only : bcast_dbl
     &                           , reduce_dbl   
     &                           , datacast_dbl
     &                           , datacast_dbl2
     &                           , translate_g2l
     &                           , get_loop_start
     &                           , get_loop_end
     &                           , get_loop_node
!----
      implicit none

      integer, intent(in)    :: nm, nv, ne, i
      real(8), intent(inout) :: a(1:nm,*)
      real(8), intent(out)   :: u_x(1:nv,*), u_y(1:nv,*)
      real(8), intent(out)   :: u_t(*), v_t(*)

      include 'trd.h'

      real(8), intent(out)   :: c(mband,mband), e(1:ne,*)

      real(8)                :: g_n(mband), a_n
      real(8)                :: t, tt(2*mband**2), ss(2*mband**2)
      real(8)                :: alpha(mband), beta(mband), scale(mband)

      real(8), parameter     :: zero = 0.0d+00, one = 1.0d+00

      integer                :: owner_node_col, col_local
      integer                :: owner_node_row
      integer                :: i_1
      integer                :: j_1, j_2, j_3
      integer                :: k_1
      integer                :: l, ll

      real(8) :: t11, t12, t22, tol
      real(8) :: a1, a2
      real(8) :: al(mband,mband), ax(mband, mband)
      real(8), parameter :: eps = 2d-16


      l  = i-mband
      ll = l-mband

      alpha(1:mband) = one
      beta (1:mband) = one
      g_n  (1:mband) = zero
      scale(1:mband) = zero

      c(1:mband, 1:mband)  = zero


      do i_1=mband,1,-1
         if ( ll + i_1 >= 1  ) then

            owner_node_row = get_loop_node(l+i_1, size_of_row, my_row)
            col_local       = translate_g2l(i-1,   size_of_col, my_col)
            call bcast_dbl(a(1,i_1), col_local, owner_node_row,
     &                     mpi_comm_row)

         end if
      end do


      al(1:mband, 1:mband) = zero

      do i_1=mband,1,-1
         if ( ll + i_1 >= 1  ) then

            owner_node_col = get_loop_node(ll+i_1, size_of_col, my_col)
            col_local      = translate_g2l(ll+i_1, size_of_col, my_col)

            do j_1=mband,1,-1
               if ( my_col == owner_node_col ) then
                  al(i_1,j_1) = a(col_local, j_1)
               endif
            end do

         end if
      end do

      j_2 = get_loop_start(1, size_of_col, my_col)
      j_3 = get_loop_end  (l, size_of_col, my_col)

      t11 = zero; t12 = zero; t22 = zero
      do j_1=j_2,j_3
         a2  =a(j_1, 2)
         a1  =a(j_1, 1)
         t22 = t22 + a2 * a2
         t12 = t12 + a1 * a2
         t11 = t11 + a1 * a1
      end do

      tt(1) =  t22
      tt(2) =  t12
      tt(3) =  t11
      tt(4) =  al(2,2)
      tt(5) =  al(1,2)
      tt(6) =  al(2,1)
      tt(7) =  al(1,1)
      call reduce_dbl(tt, ss, 7, 1, mpi_comm_col)
      ax(2,2) = tt(1)
      ax(1,2) = tt(2)
      ax(2,1) = ax(1,2)
      ax(1,1) = tt(3)
      al(2,2) = tt(4)
      al(1,2) = tt(5)
      al(2,1) = tt(6)
      al(1,1) = tt(7)

      do i_1=mband,1,-1
         if ( ll + i_1 >= 1  ) then

            j_2 = get_loop_start(1,      size_of_col, my_col)
            j_3 = get_loop_end  (ll+i_1, size_of_col, my_col)

            owner_node_col = get_loop_node(ll+i_1, size_of_col, my_col)
            col_local      = translate_g2l(ll+i_1, size_of_col, my_col)

            if ( i_1 == 2 ) then

               alpha(i_1) = ax(i_1, i_1)
               a_n        = al(i_1, i_1)

               do j_1=j_2,j_3
                  u_x(j_1, i_1) = a(j_1, i_1)
               end do
            else
               tol = 16 * ax(1, 1) * eps
               ax(1, 1) = max(zero, ax(1, 1) - al(2, 1)**2)
               ax(1, 2) =           ax(1, 2) - al(2, 1)*al(2, 2)

               alpha(i_1) = ax(i_1, i_1)
               a_n        = al(i_1, i_1)

! singularity found, so re-calculate norm and inner products
               if ( alpha(i_1) < tol ) then
                  do j_1=j_2,j_3
                     u_x(j_1, 1) = a(j_1, 1)
                  end do
                  do k_1=1,mband
                     t = zero
                     do j_1=j_2,j_3
                        t = t + u_x(j_1, i_1)*u_x(j_1, k_1)
                     enddo
                     tt(1+k_1)=t
                  end do
                  if ( my_col == owner_node_col ) then
                     a_n = a(col_local, i_1)
                  else
                     a_n = zero
                  end if
                  tt(1) =  a_n
                  call reduce_dbl(tt, ss, 1+mband, 1, mpi_comm_col)
                  a_n        = tt(1)
                  do k_1=1,mband
                     ax(i_1,k_1) = tt(1+k_1)
                  end do
                  do k_1=1,mband
                     if ( k_1 /= i_1 ) then
                        ax(k_1,i_1) = ax(i_1,k_1)
                     end if
                  end do
                  alpha(i_1) = ax(i_1,i_1)
               else
                  do j_1=j_2,j_3
                     u_x(j_1, i_1) = a(j_1, i_1)
                  end do
               end if

            endif

            if ( alpha(i_1) > zero ) then

               scale(i_1) = sqrt(alpha(i_1))
!----
! scaling for evading numerical instability
!----
               t = one / scale(i_1)
               do j_1=j_2,j_3
                  u_x(j_1, i_1) = u_x(j_1, i_1) * t
               end do
               do k_1 = 1, mband
                  al (k_1, i_1) = al (k_1, i_1) * t
                  ax (k_1, i_1) = ax (k_1, i_1) * t
                  ax (i_1, k_1) = ax (i_1, k_1) * t
               end do

               a_n          =  a_n * t
               g_n(i_1)     = -sign(one, a_n)
               beta(i_1)    =  one - a_n * g_n(i_1)
               al(i_1, i_1) =  a_n - g_n(i_1)

               if ( my_col == owner_node_col ) then
                  u_x(col_local, i_1) =  al(i_1, i_1)
               end if

               do k_1=1,i_1-1
                  ax(k_1, i_1) = ax(k_1, i_1) - g_n(i_1) * al(i_1, k_1)
                  ax(i_1, k_1) = ax(k_1, i_1)
               end do
               ax(i_1, i_1) = 2 * beta(i_1)
               do k_1=i_1+1, mband
                  ax(i_1, k_1) = ax(i_1, k_1) - g_n(i_1) * al(i_1, k_1)
                  ax(k_1, i_1) = ax(i_1, k_1)
               end do

               do k_1=1,i_1-1
                  t            = ax(k_1, i_1) / beta(i_1)
                  j_2 = get_loop_start(1,      size_of_col, my_col)
                  j_3 = get_loop_end  (ll+i_1, size_of_col, my_col)
                  do j_1=j_2,j_3
                     a (j_1, k_1) = a (j_1, k_1) - t * u_x(j_1, i_1)
                  end do! j_1
                  do j_1=1,i_1
                     al(j_1, k_1) = al(j_1, k_1) - t * al (j_1, i_1)
                  end do
                  ax(k_1, i_1) = - ax(k_1, i_1)
                  ax(i_1, k_1) =   ax(k_1, i_1)
               end do

            end if

         end if
      end do
!----
      do i_1=1,mband
         if ( ll + i_1 >= 1  ) then
!----
! scaling back
!----
            if ( alpha(i_1) > zero ) then
               g_n(i_1)   = g_n(i_1)  * scale(i_1)
               beta(i_1)  = beta(i_1) * alpha(i_1)
            else
               g_n(i_1)   = zero
               beta(i_1)  = one
            endif
            j_2 = get_loop_start(1,      size_of_col, my_col)
            j_3 = get_loop_end  (ll+i_1, size_of_col, my_col)
            t = scale(i_1)
            do j_1=j_2,j_3
               u_x(j_1, i_1) = u_x(j_1, i_1) * t
            end do
         end if
      end do
!----
      owner_node_row = get_loop_node(i, size_of_row, my_row)
      if ( owner_node_row == my_row ) then
         if ( i > 2 ) then
            owner_node_col = get_loop_node(i-2, size_of_col, my_col)
            if ( owner_node_col == my_col ) then
               col_local = translate_g2l(i-2, size_of_col, my_col)
               e(i-1, 1) = a(col_local, 1)
            end if
         end if
         if ( i > 1 ) then
            owner_node_col = get_loop_node(i-1, size_of_col, my_col)
            if ( owner_node_col == my_col ) then
               col_local = translate_g2l(i-1, size_of_col, my_col)
               e(i, 1) = a(col_local, 2)
            end if
         end if
         owner_node_col = get_loop_node(i,   size_of_col, my_col)
         if ( owner_node_col == my_col ) then
            e(i-1, 2) = g_n(1)
            e(i  , 2) = g_n(2)
         end if
      end if

      do i_1=1,mband
         if ( ll + i_1 >= 1  ) then
            j_2 = get_loop_start(ll+i_1+1, size_of_col, my_col)
            j_3 = get_loop_end  (i-1,      size_of_col, my_col)
            u_x(j_2:j_3,i_1) = zero

            j_3 = get_loop_end  (ll+i_1+mband-1, size_of_col, my_col)
            a(1:j_3,i_1) = u_x(1:j_3,i_1)
         end if
      end do
      col_local = translate_g2l(i-1, size_of_col, my_col)
      if ( mod(mband,2)==1 ) then; i_1 = 1
         call datacast_dbl(u_y(1, i_1),u_x(1, i_1),u_t,v_t,col_local)
      end if
      do i_1=mod(mband,2)+1,mband,2
         call datacast_dbl2(u_y(1, i_1),u_y(1, i_1+1),
     &                      u_x(1, i_1),u_x(1, i_1+1),u_t,v_t,col_local)
      end do

      c(1:mband, 1:mband) = zero
      do i_1=1,mband
         c(i_1, i_1) = one / beta(i_1)
      end do
      do i_1=1,mband
         do j_1=1,i_1-1
            c(j_1,i_1) = ax(j_1,i_1)*scale(j_1)*scale(i_1)
         end do
      end do

      return
      end subroutine eigen_prd_u

