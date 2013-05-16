      subroutine  eigen_trd_Au(
     &            a, nm,
     &            u_x, u_y, v_x, nv,
     &            u_t, v_t, d_t,
     &            i, i_base, m)
!----
      use communication_s, only : reduce_dbl
     &                          , get_loop_start
     &                          , get_loop_end
     &                          , get_owner_node
     &                          , translate_g2l
!----
      implicit none

      integer, intent(in)    ::  nm, nv, i, i_base, m
      real(8), intent(inout) ::  a(1:nm,*)
      real(8), intent(inout) ::  u_x(1:nv), u_y(1:nv)
      real(8), intent(inout) ::  v_x(1:nv)
      real(8), intent(inout) ::  u_t(*), v_t(*)
      real(8), intent(in)    ::  d_t(*)

      integer                ::  n, blk_0
      integer                ::  col_local, row_local
      integer                ::  n1, n2
      integer                ::  i_1, i_2, i_3
      integer                ::  j_1, j_2, j_3
      integer                ::  k_1, k_2
      integer                ::  l
      integer                ::  ierr
      integer                ::  thread_size, thread_rank

      include 'trd_au.h'

      real(8)      d1, d2
      real(8)      t2_reduce1, t2_reduce2
      common /t2/  t2_reduce1, t2_reduce2

      include 'mpif.h'
!$    include 'omp_lib.h'
      include 'trd.h'
      include 'CSTAB.h'

      thread_rank = 0
      thread_size = 1
!$    thread_rank = omp_get_thread_num()
!$    thread_size = omp_get_num_threads()

      l = i-1

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

      if ( thread_size > 1 ) then

         n1 = offset1+nv*thread_rank
         n2 = offset2+nv*thread_rank

         do j_1=1,col_local
            u0_z(j_1+n1) = 0.0d+00
         enddo
         do j_1=1,row_local
            v0_z(j_1+n2) = 0.0d+00
         enddo

         call  eigen_trd_Au_step1(
     &           a, nm,
     &           u_x, u_y, u0_z(1+n1), v0_z(1+n2),
     &           1, n, col_local, row_local, nv, blk_0
     &           ,thread_rank, thread_size
     &         )

!$omp barrier

         call  eigen_trd_Au_step2(
     &           v_x, u0_z(1+offset1),
     &           v_t, v0_z(1+offset2),
     &           nv, col_local, row_local,
     &           thread_rank, thread_size
     &         )

!$omp barrier

      else

         do j_1=1,col_local
            v_x(j_1) = 0.0d+00
         enddo
         do j_1=1,row_local
            v_t(j_1) = 0.0d+00
         enddo

         call  eigen_trd_Au_step1(
     &           a, nm,
     &           u_x, u_y, v_x, v_t,
     &           1, n, col_local, row_local, nv, blk_0,
     &           thread_rank, thread_size
     &         )

      endif

!----          if ( thread_rank == 0 ) then
!$omp master

      if ( nprocs > 1 ) then
c #ifdef TIMER
c          call mpi_barrier(mpi_comm_col,ierr)
c          d1=mpi_wtime()
c #endif

         call reduce_dbl(v_t, u_t, row_local, 1, mpi_comm_col)

c #ifdef TIMER
c          d2=mpi_wtime()
c          t2_reduce1=t2_reduce1+(d2-d1)
c #endif
      end if

      call  eigen_trd_Au_step3(
     &        u_x, u_y, v_x,
     &        u_t, v_t, d_t,
     &        1, n, col_local, row_local, nv
     &      )

      if ( nprocs > 1 ) then
c #ifdef TIMER
c          call mpi_barrier(mpi_comm_col,ierr)
c          d1=mpi_wtime()
c #endif
         call reduce_dbl(v_x, v_t, col_local, size_of_col, mpi_comm_row)

c #ifdef TIMER
c          d2=mpi_wtime()
c          t2_reduce2=t2_reduce2+(d2-d1)
c #endif
      end if

!$omp end master
!----           end if

      return
      end subroutine ! eigen_trd_Au

      subroutine  eigen_trd_Au_step1(
     &            a, nm,
     &            u_x, u_y, u_t, v_t,
     &            n1, n2, col_local, row_local, nv, blk_0,
     &            thread_rank, thread_size
     &            )
!----
      use communication_s, only : get_loop_start
     &                          , get_loop_end
     &                          , translate_l2g
!----
      implicit none

      integer, intent(in)    :: nm, nv, n1, n2
      real(8), intent(inout) :: a(1:nm,*)
      real(8), intent(inout) :: u_x(1:nv), u_y(1:nv)
      real(8), intent(inout) :: u_t(1:nv), v_t(1:nv)
      integer, intent(in)    :: blk_0
      integer                :: col_local, row_local
      integer                :: thread_size, thread_rank

      integer                :: i_0
      integer                :: i_1, i_2, i_3, i_4
      integer                :: j_1, j_2, j_3, j_4
      integer                :: k_1, k_2, k_3, k_4
      integer                :: l_1, l_2, l_3, l_4
      integer                :: i, j, k

      real(8)                :: v0, u0
      real(8)                :: a0_0
      real(8)                :: a0_1
      real(8)                :: a0_2
      real(8)                :: a0_3
      real(8)                :: a0_4
      real(8)                :: a0_5
      real(8)                :: u_y0
      real(8)                :: u_y1
      real(8)                :: u_y2
      real(8)                :: u_y3
      real(8)                :: u_y4
      real(8)                :: u_y5

      integer                :: lx, ly
      integer                :: ii_1, ii_2, ii_3, ii_4, ii_5
      integer                :: jj_1, jj_2, jj_3, jj_4, jj_5
      integer                :: kk_1, kk_2, kk_3, kk_4, kk_5

      include 'mpif.h'
      include 'trd.h'

          i_2 = n1
          i_3 = n2
!----
!---- v:= au
!----
      if ( blk_0 == 0 ) then
         do i_1=i_2+thread_rank,i_3,thread_size
            j   = translate_l2g(i_1, size_of_row,my_row)
            j_3 = get_loop_end  (j, size_of_col,my_col)
            j   = j+size_of_row*(6*2-1)
            j_4 = get_loop_end  (j, size_of_col,my_col)
            do j_1=j_3+1,min(j_4,nm)
               a(j_1,i_1) = 0.0d+00
            end do! j_1
         end do! i_1
!$omp barrier
      end if
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

               u_y0 = u_y(i_1+0)

!dir$ ivdep
!dir$ vector aligned
               do j_1=kk_1,kk_4
!              do j_1=jj_2,jj_3

                  u0 = u_x(j_1+0)
                  a0_0 = a(j_1+0,i_1+0)
                  v0 = u_t(j_1+0)

                  v0 = v0
     &                            + (a0_0*u_y0)
                  v_t(i_1+0) = v_t(i_1+0)
     &                            + (a0_0*u0)
                  u_t(j_1+0) = v0

               end do! j_1

            end do! kk_1

         end do! i_1

         do i_0=ii_4+thread_rank*6*2,ii_3, thread_size*6*2

            j    = translate_l2g(i_0, size_of_row,my_row)
            j    = j+(6*2-1)*size_of_row
            jj_2 = k_1
            jj_3 = get_loop_end  (j-1, size_of_col,my_col)
            jj_3 = min(k_4,jj_3)
            if ( jj_2 > jj_3 ) cycle

            do kk_1=jj_2,jj_3,16*31; kk_4=min(kk_1+16*31-1,jj_3)
               do i_1 = i_0,i_0+6*2-1,6

                  u_y0 = u_y(i_1+0)
                  u_y1 = u_y(i_1+1)
                  u_y2 = u_y(i_1+2)
                  u_y3 = u_y(i_1+3)
                  u_y4 = u_y(i_1+4)
                  u_y5 = u_y(i_1+5)

!dir$ ivdep
!dir$ vector aligned
!dir$ unroll(8)
                  do j_1=kk_1,kk_4

                     u0 = u_x(j_1+0)

                     a0_0 = a(j_1+0,i_1+0)
                     a0_1 = a(j_1+0,i_1+1)

                     v_t(i_1+0) = v_t(i_1+0)
     &                            + (a0_0*u0)
                     v_t(i_1+1) = v_t(i_1+1)
     &                            + (a0_1*u0)

                  enddo

!---- ccc            call mm_prefetch(a(kk_1,i_1),1)

!dir$ ivdep
!dir$ vector aligned
!dir$ unroll(4)
                  do j_1=kk_1,kk_4

                     u0 = u_x(j_1+0)
                     v0 = 0.0d+00

                     a0_0 = a(j_1+0,i_1+0)
                     a0_1 = a(j_1+0,i_1+1)

                     v0 = v0
     &                            + (a0_0*u_y0)
     &                            + (a0_1*u_y1)

                     a0_2 = a(j_1+0,i_1+2)
                     a0_3 = a(j_1+0,i_1+3)

                     v_t(i_1+2) = v_t(i_1+2)
     &                            + (a0_2*u0)
                     v_t(i_1+3) = v_t(i_1+3)
     &                            + (a0_3*u0)

                     v0 = v0
     &                            + (a0_2*u_y2)
     &                            + (a0_3*u_y3)

                     a0_4 = a(j_1+0,i_1+4)
                     a0_5 = a(j_1+0,i_1+5)

                     v_t(i_1+4) = v_t(i_1+4)
     &                            + (a0_4*u0)
                     v_t(i_1+5) = v_t(i_1+5)
     &                            + (a0_5*u0)

                     v0 = v0
     &                            + (a0_4*u_y4)
     &                            + (a0_5*u_y5)

                     u_t(j_1+0) = v0 + u_t(j_1+0)

                  end do! j_1

               enddo! i_1
            enddo! kk_1

         end do! i_1
      end do! k_1

      return
      end subroutine ! eigen_trd_Au_step1

      subroutine  eigen_trd_Au_step2(
     &            u_t, u_z, v_t, v_z, nv,
     &            col_local, row_local,
     &            thread_rank, thread_size
     &            )
      implicit none

      integer, intent(in)    :: nv, col_local, row_local
      real(8), intent(out)   :: u_t(1:nv),   v_t(1:nv)
      real(8), intent(in)    :: u_z(1:nv,*), v_z(1:nv,*)

      integer                :: thread_rank, thread_size

      integer                :: j_1, j_2, j_3, j_4
      integer                :: jj_1, jj_2, jj_3, jj_4
      integer                :: i, j, k

      integer, parameter     :: lx = 1024

      jj_1 = col_local
      jj_2 = max(512,(jj_1-1)/thread_size+1)
      jj_3 =    (jj_2*(thread_rank+0)     )+1
      jj_4 = min(jj_2*(thread_rank+1),jj_1)

      do jj_1=jj_3,jj_4,lx
         j_3=jj_1; j_4=min(jj_1+lx-1,jj_4)
         if ( thread_size == 4 ) then
            do j_1=j_3,j_4
               u_t(j_1) = u_z(j_1,1)+u_z(j_1,2)+u_z(j_1,3)+u_z(j_1,4)
            end do
         else
            if ( mod(thread_size,2) == 1 ) then
               do j_1=j_3,j_4
                  u_t(j_1) = u_z(j_1,1)
               end do
            else
               do j_1=j_3,j_4
                  u_t(j_1) = 0.0d+00
               end do
            end if
            do j=mod(thread_size,2)+1,thread_size,2
               do j_1=j_3,j_4
                  u_t(j_1) = u_t(j_1)+u_z(j_1,j+0)+u_z(j_1,j+1)
               end do
            end do
         end if
      end do

      jj_1 = row_local
      jj_2 = max((jj_1-1)/thread_size+1,512)
      jj_3 =    (jj_2*(thread_rank+0)     )+1
      jj_4 = min(jj_2*(thread_rank+1),jj_1)

      do jj_1=jj_3,jj_4,lx
         j_3=jj_1; j_4=min(jj_1+lx-1,jj_4)
         if ( thread_size == 4 ) then
            do j_1=j_3,j_4
               v_t(j_1) = v_z(j_1,1)+v_z(j_1,2)+v_z(j_1,3)+v_z(j_1,4)
            end do
         else
            if ( mod(thread_size,2) == 1 ) then
               do j_1=j_3,j_4
                  v_t(j_1) = v_z(j_1,1)
               end do
            else
               do j_1=j_3,j_4
                  v_t(j_1) = 0.0d0
               end do
            end if
            do j=mod(thread_size,2)+1,thread_size,2
               do j_1=j_3,j_4
                  v_t(j_1) = v_t(j_1)+v_z(j_1,j+0)+v_z(j_1,j+1)
               end do
            end do
         end if
      end do

      return
      end subroutine ! eigen_trd_Au_step2

      subroutine  eigen_trd_Au_step3(
     &            u_x, u_y, v_x,
     &            u_t,v_t,d_t,
     &            n1, n2, col_local, row_local, nv
     &            )
!----
      use communication_s, only : translate_l2g
     &                          , translate_g2l
!----
      implicit none

      integer, intent(in)    :: nv, n1, n2
      real(8), intent(in)    :: u_x(1:nv),u_y(1:nv)
      real(8), intent(inout) :: v_x(1:nv)
      real(8), intent(in)    :: u_t(*),v_t(*)
      real(8), intent(in)    :: d_t(1:nv)
      integer                :: col_local, row_local

      integer                :: i_1, i_2, i_3, i_4
      integer                :: i, j, k
      integer                :: nm1, nm2

      include 'mpif.h'
      include 'trd.h'

      i_2 = n1
      i_3 = n2

      if ( diag_0 > 0 ) then

         j = translate_l2g(diag_0, size_of_row,my_row)
         j = translate_g2l(j, size_of_col,my_col)
         if ( j > nv ) return

            nm1 = size_of_row/n_common
            nm2 = size_of_col/n_common

            if ( nm2 == 1 ) then
               if ( nm1 == 1 ) then
                  call eigen_trd_Au_step3_sub3(
     &                   v_x(j), v_t(diag_0), d_t(diag_0), u_y(diag_0),
     &                   (i_3-diag_0)/(size_of_col/n_common)+1,nm1,nm2
     &                 )
               else
                  call eigen_trd_Au_step3_sub2(
     &                   v_x(j), v_t(diag_0), d_t(diag_0), u_y(diag_0),
     &                   (i_3-diag_0)/(size_of_col/n_common)+1,nm1,nm2
     &                 )
               end if
            else
               call eigen_trd_Au_step3_sub1(
     &                v_x(j), v_t(diag_0), d_t(diag_0), u_y(diag_0),
     &                (i_3-diag_0)/(size_of_col/n_common)+1,nm1,nm2
     &              )

            end if

         end if

      return
      end subroutine ! eigen_trd_Au_step3

      subroutine eigen_trd_Au_step3_sub1(v_x,v_t,d_t,u_y, n,nm1,nm2)
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
      end subroutine ! eigen_trd_Au_step3_sub1

      subroutine eigen_trd_Au_step3_sub2(v_x,v_t,d_t,u_y, n,nm1,nm2)
      implicit none

      integer, intent(in)    :: n, nm1, nm2
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
      end subroutine ! eigen_trd_Au_step3_sub2

      subroutine eigen_trd_Au_step3_sub3(v_x,v_t,d_t,u_y, n,nm1,nm2)
      implicit none

      integer, intent(in)    :: n, nm1, nm2
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
      end subroutine ! eigen_trd_Au_step3_sub3

