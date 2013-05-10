      subroutine  eigen_prd_Au(
     &               a, nm,
     &               u_x, u_y, v_x, nv,
     &               u_t, v_t, d_t,
     &               i, i_base, m)
!----
      use communication_sx, only : reduce_dbl
     &                           , get_loop_start
     &                           , get_loop_end
     &                           , translate_g2l
!----
!$    use omp_lib
      implicit none

      integer, intent(in)    :: nm, nv, i, i_base, m
      real(8), intent(inout) :: a(1:nm,*)
      real(8), intent(inout) :: u_x(1:nv,*), u_y(1:nv,*)
      real(8), intent(inout) :: v_x(1:nv,*)
      real(8), intent(inout) :: u_t(*), v_t(1:nv,*)
      real(8), intent(in)    :: d_t(*)

      integer                :: n, blk_0
      integer                :: col_local, row_local
      integer                :: n1, n2, n3, n4
      integer                :: i_1, i_2, i_3
      integer                :: j_1, j_2, j_3
      integer                :: k_1, k_2, k_3
      integer                :: l

      integer                :: thread_size, thread_rank

      real(8), pointer       :: u0_z(:),v0_z(:)
      real(8), pointer       :: u1_z(:),v1_z(:)
      integer                :: offset1, offset2
      integer                :: offset3, offset4
      common/tred_au_common/    u0_z, v0_z, u1_z, v1_z,
     &                          offset1, offset2, offset3, offset4

      include 'trd.h'
      include 'param.h'


      thread_rank = 0
      thread_size = 1
!$    thread_rank = omp_get_thread_num()
!$    thread_size = omp_get_num_threads()

      l = i-2

      i_2 = get_loop_start(1, size_of_row,my_row)
      i_3 = get_loop_end  (l, size_of_row,my_row)
      j_2 = get_loop_start(1, size_of_col,my_col)
      j_3 = get_loop_end  (l, size_of_col,my_col)

      row_local  = translate_g2l(i, size_of_row,my_row)
      col_local  = translate_g2l(l, size_of_col,my_col)

      k_1   = i-i_base
      k_2   = m
      n     = i_3
      blk_0 = k_1-k_2


      n1 = offset1+nv*thread_rank
      n2 = offset2+nv*thread_rank
      n3 = offset3+nv*thread_rank
      n4 = offset4+nv*thread_rank

      do j_1=1,col_local
         u0_z(j_1+n1) = zero
         u1_z(j_1+n3) = zero
      end do! j_1
      do j_1=1,row_local
         v0_z(j_1+n2) = zero
         v1_z(j_1+n4) = zero
      end do! j_1

      call  eigen_prd_Au_step1(
     &      a, nm,
     &      u_x(1,1), u_y(1,1),
     &      u_x(1,2), u_y(1,2),
     &      u0_z(1+n1), v0_z(1+n2),
     &      u1_z(1+n3), v1_z(1+n4),
     &      1, n, nv, blk_0
     &      ,thread_rank, thread_size
     &      )

!$omp barrier

      call  eigen_prd_Au_step2(
     &      v_x(1,1), u0_z(1+offset1), u1_z(1+offset3),
     &      v_t(1,1), v0_z(1+offset2), v1_z(1+offset4),
     &      nv, col_local, row_local
     &      ,thread_rank, thread_size
     &      )

!$omp barrier

!$omp master

      if ( nprocs > 1 ) then
!----    call reduce_dbl(v_t(1,1), u_t, row_local, 1, mpi_comm_col)
!----    call reduce_dbl(v_t(1,2), u_t, row_local, 1, mpi_comm_col)
         u_t(1:row_local) = v_t(1:row_local,1)
         u_t(1+row_local:2*row_local) = v_t(1:row_local,2)
         call reduce_dbl(u_t, v_t, 2*row_local, 1, mpi_comm_col)
         v_t(1:row_local,1) = u_t(1:row_local)
         v_t(1:row_local,2) = u_t(1+row_local:2*row_local)
      end if

      call  eigen_prd_Au_step3(
     &      u_y, v_x,
     &      v_t, d_t,
     &      1, n, nv
     &      )

      if ( nprocs > 1 ) then
         v_t(1:col_local,1) = v_x(1:col_local,1)
         v_t(1+col_local:2*col_local,1) = v_x(1:col_local,2)
         call reduce_dbl(v_t,u_t,2*col_local,size_of_col,mpi_comm_row)
         v_x(1:col_local,1) = v_t(1:col_local,1)
         v_x(1:col_local,2) = v_t(1+col_local:2*col_local,1)
      end if

!$omp end master


      return
      end subroutine ! eigen_prd_Au

      subroutine  eigen_prd_Au_step1(
     &            a, nm,
     &            u0_x, u0_y, u1_x, u1_y,
     &            u0_t, v0_t, u1_t, v1_t,
     &            n1, n2, nv, blk_0
     &            ,thread_rank, thread_size
     &            )
!----
      use communication_sx, only : get_loop_start
     &                           , get_loop_end
     &                           , translate_l2g
!----
      implicit none

      integer, intent(in)    :: nm, nv, n1, n2
      real(8), intent(inout) :: a(1:nm,*)
      real(8), intent(inout) :: u0_x(1:nv), u0_y(1:nv)
      real(8), intent(inout) :: u1_x(1:nv), u1_y(1:nv)
      real(8), intent(inout) :: u0_t(1:nv), v0_t(1:nv)
      real(8), intent(inout) :: u1_t(1:nv), v1_t(1:nv)
      integer, intent(in)    :: blk_0
      integer                :: thread_size, thread_rank

      integer                :: i_0, i_1, i_2, i_3, i_4
      integer                :: j_1, j_2, j_3, j_4
      integer                :: k_1, k_2, k_3, k_4
      integer                :: l_1, l_2, l_3, l_4
      integer                :: i, j, k

      real(8)                :: v0, u0
      real(8)                :: v1, u1
      real(8)                :: a0_0
      real(8)                :: a0_1
      real(8)                :: a0_2
      real(8)                :: a0_3
      real(8)                :: a0_4
      real(8)                :: a0_5
      real(8)                :: u0_y0
      real(8)                :: u1_y0
      real(8)                :: u0_y1
      real(8)                :: u1_y1
      real(8)                :: u0_y2
      real(8)                :: u1_y2
      real(8)                :: u0_y3
      real(8)                :: u1_y3
      real(8)                :: u0_y4
      real(8)                :: u1_y4
      real(8)                :: u0_y5
      real(8)                :: u1_y5

      integer                ::   lx, ly
      integer                ::   ii_2, ii_3, ii_4
      integer                ::   jj_2, jj_3
      integer                ::   kk_1, kk_4

      include 'trd.h'
      include 'param.h'


      i_2 = n1
      i_3 = n2

!----
! v:= au
!----
      if ( blk_0 == 0 ) then
         do i_1=i_2+thread_rank,i_3,thread_size
            j   = translate_l2g(i_1, size_of_row,my_row)
            j_3 = get_loop_end  (j, size_of_col,my_col)
            j   = j+size_of_row*(6*2-1)
            j_4 = get_loop_end  (j, size_of_col,my_col)
            do j_1=j_3+1,min(j_4,nm)
               a(j_1,i_1) = zero
            end do! j_1
         end do! i_1
!$omp barrier
      end if
!----
!----
      lx = 32*1000; ly = 32*1
!----
      l_2 = i_2; l_3 = i_3
      l_1 = l_2; l_4 = l_3

      k_2 = 1
      k   = translate_l2g(l_4, size_of_row,my_row)
      k_3 = get_loop_end  (k-1, size_of_col,my_col)
      do k_1=k_2,k_3,lx; k_4 = min(k_3,k_1+lx-1)

         j    = translate_l2g(k_1, size_of_col,my_col)
         ii_2 = get_loop_start(j,size_of_row,my_row)
         ii_2 = max(l_1,ii_2)
         ii_3 = l_4
         if ( ii_2 > ii_3 ) cycle
         ii_4 = mod(ii_3-ii_2+1,6*2)+ii_2

         do i_1=ii_2+thread_rank,ii_4-1,thread_size

            j    = translate_l2g(i_1, size_of_row,my_row)
            j    = j+(1-1)*size_of_row
            jj_2 = k_1
            jj_3 = get_loop_end  (j-1, size_of_col,my_col)
            jj_3 = min(k_4,jj_3)
            if ( jj_2 > jj_3 ) cycle

            do kk_1=jj_2,jj_3,64*64; kk_4=min(kk_1+64*64-1,jj_3)

!dir$ ivdep
!dir$ vector aligned
               do j_1=kk_1,kk_4
!----          do j_1=jj_2,jj_3

                  u0 = u0_x(j_1+0)
                  u1 = u1_x(j_1+0)

                  a0_0 = a(j_1+0,i_1+0)
                  v0_t(i_1+0) = v0_t(i_1+0)
     &                        + (a0_0*u0)
                  v1_t(i_1+0) = v1_t(i_1+0)
     &                        + (a0_0*u1)

               end do! j_1

               u0_y0 = u0_y(i_1+0)
               u1_y0 = u1_y(i_1+0)

!dir$ ivdep
!dir$ vector aligned
               do j_1=kk_1,kk_4
!----          do j_1=jj_2,jj_3

                  v0 = u0_t(j_1+0)
                  v1 = u1_t(j_1+0)

                  a0_0 = a(j_1+0,i_1+0)
                  v0 = v0
     &                        + (a0_0*u0_y0)
                  v1 = v1
     &                        + (a0_0*u1_y0)
                  u0_t(j_1+0) = v0
                  u1_t(j_1+0) = v1

               end do! j_1

            enddo ! kk_1

         end do! i_1

         do i_0=ii_4+thread_rank*6*2,ii_3, thread_size*6*2

            j    = translate_l2g(i_0, size_of_row,my_row)
            j    = j+(6*2-1)*size_of_row
            jj_2 = k_1
            jj_3 = get_loop_end  (j-1, size_of_col,my_col)
            jj_3 = min(k_4,jj_3)
            if ( jj_2 > jj_3 ) cycle

            do kk_1=jj_2,jj_3,336; kk_4=min(kk_1+336-1,jj_3)
               do i_1 = i_0,i_0+6*2-1,6

!dir$ ivdep
!dir$ vector aligned
                  do j_1=kk_1,kk_4
!----             do j_1=jj_2,jj_3

                     u0 = u0_x(j_1+0)
                     u1 = u1_x(j_1+0)

                     a0_0 = a(j_1+0,i_1+0)
                     a0_1 = a(j_1+0,i_1+1)

                     v0_t(i_1+0) = v0_t(i_1+0)
     &                           + (a0_0*u0)
                     v1_t(i_1+0) = v1_t(i_1+0)
     &                           + (a0_0*u1)
                     v0_t(i_1+1) = v0_t(i_1+1)
     &                           + (a0_1*u0)
                     v1_t(i_1+1) = v1_t(i_1+1)
     &                           + (a0_1*u1)

                  end do! j_1

                  u0_y0 = u0_y(i_1+0)
                  u0_y1 = u0_y(i_1+1)
                  u0_y2 = u0_y(i_1+2)
                  u0_y3 = u0_y(i_1+3)
                  u0_y4 = u0_y(i_1+4)
                  u0_y5 = u0_y(i_1+5)
                  u1_y0 = u1_y(i_1+0)
                  u1_y1 = u1_y(i_1+1)
                  u1_y2 = u1_y(i_1+2)
                  u1_y3 = u1_y(i_1+3)
                  u1_y4 = u1_y(i_1+4)
                  u1_y5 = u1_y(i_1+5)


!dir$ ivdep
!dir$ vector aligned
                  do j_1=kk_1,kk_4
!----             do j_1=jj_2,jj_3

                     u0 = u0_x(j_1+0)
                     u1 = u1_x(j_1+0)

                     v0 = u0_t(j_1+0)
                     v1 = u1_t(j_1+0)

                     a0_0 = a(j_1+0,i_1+0)
                     a0_1 = a(j_1+0,i_1+1)

                     v0 = v0
     &                  + (a0_0*u0_y0)
     &                  + (a0_1*u0_y1)
                     v1 = v1
     &                  + (a0_0*u1_y0)
     &                  + (a0_1*u1_y1)

                     a0_2 = a(j_1+0,i_1+2)
                     a0_3 = a(j_1+0,i_1+3)

                     v0_t(i_1+2) = v0_t(i_1+2)
     &                           + (a0_2*u0)
                     v1_t(i_1+2) = v1_t(i_1+2)
     &                           + (a0_2*u1)
                     v0_t(i_1+3) = v0_t(i_1+3)
     &                           + (a0_3*u0)
                     v1_t(i_1+3) = v1_t(i_1+3)
     &                           + (a0_3*u1)

                     v0 = v0
     &                  + (a0_2*u0_y2)
     &                  + (a0_3*u0_y3)
                     v1 = v1
     &                  + (a0_2*u1_y2)
     &                  + (a0_3*u1_y3)

                     a0_4 = a(j_1+0,i_1+4)
                     a0_5 = a(j_1+0,i_1+5)

                     v0_t(i_1+4) = v0_t(i_1+4)
     &                           + (a0_4*u0)
                     v1_t(i_1+4) = v1_t(i_1+4)
     &                           + (a0_4*u1)
                     v0_t(i_1+5) = v0_t(i_1+5)
     &                           + (a0_5*u0)
                     v1_t(i_1+5) = v1_t(i_1+5)
     &                           + (a0_5*u1)

                     v0 = v0
     &                  + (a0_4*u0_y4)
     &                  + (a0_5*u0_y5)
                     v1 = v1
     &                  + (a0_4*u1_y4)
     &                  + (a0_5*u1_y5)

                     u0_t(j_1+0) = v0
                     u1_t(j_1+0) = v1

                  end do! j_1

               end do! i_1
            end do! kk_1

         end do! i_0

      end do! k_1


      return
      end subroutine ! eigen_prd_Au_step1

      subroutine  eigen_prd_Au_step2(
     &            u_t, u0_z, u1_z, v_t, v0_z, v1_z, nv,
     &            col_local, row_local
     &            ,thread_rank, thread_size
     &            )
      implicit none

      integer, intent(in)    :: nv, col_local, row_local
      real(8), intent(out)   :: u_t(1:nv,*),  v_t(1:nv,*)
      real(8), intent(in)    :: u0_z(1:nv,*), v0_z(1:nv,*)
      real(8), intent(in)    :: u1_z(1:nv,*), v1_z(1:nv,*)

      integer                :: thread_rank, thread_size

      integer                :: j_1, j_3, j_4
      integer                :: jj_1, jj_2, jj_3, jj_4
      integer                :: j

      integer, parameter     :: lx = 1024


      jj_1 = col_local
      jj_2 = (jj_1-1)/thread_size+1
      jj_3 =    (jj_2*(thread_rank+0)     )+1
      jj_4 = min(jj_2*(thread_rank+1),jj_1)

      do jj_1=jj_3,jj_4,lx
         j_3=jj_1; j_4=min(jj_1+lx-1,jj_4)
         if ( mod(thread_size,2) == 1 ) then
            do j_1=j_3,j_4
               u_t(j_1,1) = u0_z(j_1,1)
               u_t(j_1,2) = u1_z(j_1,1)
            end do
         else
            do j_1=j_3,j_4
               u_t(j_1,1) = 0.0d+00
               u_t(j_1,2) = 0.0d+00
            end do
         end if
         do j=mod(thread_size,2)+1,thread_size,2
            do j_1=j_3,j_4
               u_t(j_1,1) = u_t(j_1,1)+u0_z(j_1,j+0)+u0_z(j_1,j+1)
               u_t(j_1,2) = u_t(j_1,2)+u1_z(j_1,j+0)+u1_z(j_1,j+1)
            end do
         end do
      end do

      jj_1 = row_local
      jj_2 = (jj_1-1)/thread_size+1
      jj_3 =    (jj_2*(thread_rank+0)     )+1
      jj_4 = min(jj_2*(thread_rank+1),jj_1)

      do jj_1=jj_3,jj_4,lx
         j_3=jj_1; j_4=min(jj_1+lx-1,jj_4)
         if ( mod(thread_size,2) == 1 ) then
            do j_1=j_3,j_4
               v_t(j_1,1) = v0_z(j_1,1)
               v_t(j_1,2) = v1_z(j_1,1)
            end do
         else
            do j_1=j_3,j_4
               v_t(j_1,1) = 0.0d+00
               v_t(j_1,2) = 0.0d+00
            end do
         end if
         do j=mod(thread_size,2)+1,thread_size,2
            do j_1=j_3,j_4
               v_t(j_1,1) = v_t(j_1,1)+v0_z(j_1,j+0)+v0_z(j_1,j+1)
               v_t(j_1,2) = v_t(j_1,2)+v1_z(j_1,j+0)+v1_z(j_1,j+1)
            end do
         end do
      end do

      return
      end subroutine ! eigen_prd_Au_step2

      subroutine  eigen_prd_Au_step3(
     &            u_y, v_x,
     &            v_t, d_t,
     &            n1, n2, nv
     &            )
!----
      use communication_sx, only : translate_g2l
     &                           , translate_l2g
!----
      implicit none

      integer, intent(in)    :: nv, n1, n2
      real(8), intent(in)    :: u_y(1:nv,*)
      real(8), intent(inout) :: v_x(1:nv,*)
      real(8), intent(in)    :: v_t(1:nv,*)
      real(8), intent(in)    :: d_t(1:nv)

      integer                :: i_2, i_3
      integer                :: i, j
      integer                :: nm1, nm2

      include 'trd.h'


      i_2 = n1
      i_3 = n2

      if ( diag_0 > 0 ) then

         j = translate_l2g(diag_0, size_of_row,my_row)
         j = translate_g2l(j, size_of_col,my_col)
         if ( j > nv ) return

         nm1 = size_of_row/n_common
         nm2 = size_of_col/n_common

         do i=1,mband

            if ( nm2 == 1 ) then
               if ( nm1 == 1 ) then
                  call eigen_prd_Au_step3_sub3(
     &               v_x(j,i), v_t(diag_0,i),
     &               d_t(diag_0), u_y(diag_0,i),
     &               (i_3-diag_0)/(size_of_col/n_common)+1)
               else
                  call eigen_prd_Au_step3_sub2(
     &               v_x(j,i), v_t(diag_0,i),
     &               d_t(diag_0), u_y(diag_0,i),
     &               (i_3-diag_0)/(size_of_col/n_common)+1, nm1)
               end if
            else
               call eigen_prd_Au_step3_sub1(
     &            v_x(j,i), v_t(diag_0,i), d_t(diag_0), u_y(diag_0,i),
     &            (i_3-diag_0)/(size_of_col/n_common)+1, nm1, nm2)
            end if

         end do! i

      end if


      return
      end subroutine ! eigen_prd_Au_step3

      subroutine eigen_prd_Au_step3_sub2(v_x,v_t,d_t,u_y, n,nm1)
      implicit none

      integer, intent(in)    :: n, nm1
      real(8), intent(inout) :: v_x(nm1,*)
      real(8), intent(in)    :: v_t(*)
      real(8), intent(in)    :: d_t(*)
      real(8), intent(in)    :: u_y(*)

      integer                :: i


*soption unroll(4)
!dir$ vector always
      do i=1,n
         v_x(1,i) = v_x(1,i)+v_t(i)+d_t(i)*u_y(i)
      end do! i


      return
      end subroutine ! eigen_prd_Au_step3_sub2

      subroutine eigen_prd_Au_step3_sub3(v_x,v_t,d_t,u_y, n)
      implicit none

      integer, intent(in)    :: n
      real(8), intent(inout) :: v_x(1:n)
      real(8), intent(in)    :: v_t(1:n)
      real(8), intent(in)    :: d_t(1:n)
      real(8), intent(in)    :: u_y(1:n)

      integer                :: i


*soption unroll(4)
!dir$ vector always
      do i=1,n
         v_x(i) = v_x(i)+v_t(i)+d_t(i)*u_y(i)
      end do! i


      return
      end subroutine ! eigen_prd_Au_step3_sub3

      subroutine eigen_prd_Au_step3_sub1(v_x,v_t,d_t,u_y, n,nm1,nm2)
      implicit none

      integer, intent(in)    :: n, nm1, nm2
      real(8), intent(inout) :: v_x(nm1,*)
      real(8), intent(in)    :: v_t(nm2,*)
      real(8), intent(in)    :: d_t(nm2,*)
      real(8), intent(in)    :: u_y(nm2,*)

      integer                :: i


*soption unroll(4)
!dir$ vector always
      do i=1,n
         v_x(1,i) = v_x(1,i)+v_t(1,i)+d_t(1,i)*u_y(1,i)
      end do! i


      return
      end subroutine ! eigen_prd_Au_step3_sub1
