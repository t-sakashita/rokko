      module communication_sx

      implicit none

      include 'mpif.h'

      private :: isend_dbl
      private :: recv_dbl
      private :: irecv_dbl
      private :: wait_dbl
      private :: waitall_dbl
      public  :: bcast_dbl
      public  :: reduce_dbl
      public  :: get_loop_start
      public  :: get_loop_end
      public  :: translate_l2g
      public  :: translate_g2l
      public  :: get_loop_node
      public  :: eigen_init
      public  :: eigen_free
      public  :: datacast_dbl
      public  :: datacast_dbl2
      private :: ierr

      integer :: ierr

      contains
!----
      subroutine isend_dbl(buf, n, idest, ireq, icom)
      implicit none

      integer, intent(in)    :: n, idest, icom
      integer, intent(inout) :: ireq
      real(8), intent(inout) :: buf(1:n)

      call mpi_isend(buf,n,mpi_double_precision,
     &               idest-1, 1, icom, ireq, ierr)


      return
      end subroutine ! isend_dbl
!----
      subroutine recv_dbl(buf, n, isrc, icom)
      implicit none

      integer, intent(in)    :: n, isrc, icom
      real(8), intent(inout) :: buf(1:n)

      integer                :: status(mpi_status_size)

      call mpi_recv(buf, n, mpi_double_precision,
     &              isrc-1, 1, icom, status, ierr)

      return
      end subroutine ! recv_dbl
!----
      subroutine irecv_dbl(buf, n, isrc, ireq, icom)
      implicit none

      integer, intent(in)    :: n, isrc, icom
      integer, intent(inout) :: ireq
      real(8), intent(inout) :: buf(1:n)

      call mpi_irecv(buf, n, mpi_double_precision,
     &               isrc-1, 1, icom, ireq, ierr)

      return
      end subroutine ! irecv_dbl
!----
      subroutine wait_dbl(ireq)
      implicit none

      integer, intent(inout) :: ireq

      integer                :: status(mpi_status_size)

      call mpi_wait(ireq, status, ierr)

      return
      end subroutine ! wait_dbl
!----
      subroutine waitall_dbl(n, ireq)
      implicit none

      integer, intent(in   ) :: n
      integer, intent(inout) :: ireq(n)

      integer, pointer       :: status(:)

      allocate(status(mpi_status_size*n))
      call mpi_waitall(n, ireq, status(1), ierr)
      deallocate(status)

      return
      end subroutine ! waitall_dbl
!----
      subroutine bcast_dbl(buf, n, iroot, icom)
      implicit none

      integer, intent(in)    :: n, iroot, icom
      real(8), intent(inout) :: buf(1:n)

      integer :: my_rank, i, j

      real(8) :: d1,d2
      real(8) ::    time_bcast, time_reduce, time_redist, time_gather
      common /stat/ time_bcast, time_reduce, time_redist, time_gather

      d1=mpi_wtime()

      call mpi_bcast(buf, n, mpi_double_precision,
     &               iroot-1, icom, ierr)

!     call d_bbcast(buf, n, iroot-1, icom, ierr)

      d2=mpi_wtime()
      time_bcast=time_bcast+(d2-d1)

      return
      end subroutine ! bcast_dbl
!----
      subroutine reduce_dbl(buf, wrk, n, dist, icom)
      implicit none

      integer, intent(in)    :: n, dist, icom
      real(8), intent(inout) :: buf(1:n), wrk(1:n)

      real(8) :: d1,d2
      real(8) ::    time_bcast, time_reduce, time_redist, time_gather
      common /stat/ time_bcast, time_reduce, time_redist, time_gather

      d1=mpi_wtime()

      call mpi_allreduce(buf, wrk, n, mpi_double_precision,
     &                   mpi_sum, icom, ierr)
      buf(1:n) = wrk(1:n)
!---- call mpi_allreduce(mpi_in_place, buf, n, mpi_double_precision,
!----&                   mpi_sum, icom, ierr)

      d2=mpi_wtime()
      time_reduce=time_reduce+(d2-d1)

      return
      end subroutine ! reduce_dbl
!----
      integer function get_loop_start(istart, nnod, inod)
      implicit none

      integer, intent(in)    :: istart, nnod, inod

      get_loop_start = (istart+nnod-1-inod)/nnod+1

      return
      end function ! get_loop_start
!----
      integer function get_loop_end(iend, nnod, inod)
      implicit none

      integer, intent(in)    :: iend, nnod, inod

      get_loop_end   = (iend  +nnod  -inod)/nnod+0

      return
      end function ! get_loop_end
!----
      integer function translate_l2g(ictr, nnod, inod)
      implicit none

      integer, intent(in)    :: ictr, nnod, inod

      translate_l2g   = (ictr-1)*nnod+inod

      return
      end function ! get_loop_end
!----
      integer function get_loop_node(ictr, nnod, inod)
      implicit none

      integer, intent(in)    :: ictr, nnod, inod

      get_loop_node   = mod(ictr-1,nnod)+1

      return
      end function ! get_loop_end
!----
      integer function translate_g2l(ictr, nnod, inod)
      implicit none

      integer, intent(in)    :: ictr, nnod, inod

      translate_g2l   = (ictr-1)/nnod+1

      return
      end function ! get_loop_end
!----
      subroutine eigen_init(ndim)
      implicit none

      include 'trd.h'

      integer ::  ndim
      integer ::  n1, n2, n3, i0, i1, j0, i, k, ierr

      real(8) ::    time_bcast, time_reduce, time_redist, time_gather
      common /stat/ time_bcast, time_reduce, time_redist, time_gather

      time_bcast = 0d0
      time_reduce= 0d0
      time_redist= 0d0
      time_gather= 0d0

      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      call mpi_comm_rank(mpi_comm_world, myrank, ierr)
      myrank = myrank+1

      if ( ndim == 1 ) then
         size_of_col = 1 !! fixed
      else
         size_of_col = int(sqrt(dble(nprocs)))
         i = 1
         if ( mod(nprocs,i) == 0 ) then
            k = i
         else
            k = 1
         end if
         do
            if ( size_of_col <= k ) exit
            if ( mod(size_of_col,k) == 0 .and.
     &           mod(nprocs,size_of_col) == 0 ) exit
            size_of_col = size_of_col-1
         end do!!
      endif
      size_of_row = nprocs/size_of_col
!----
!---- my_col =    (myrank-1)/size_of_row +1
!---- my_row = mod(myrank-1, size_of_row)+1
      my_col = mod(myrank-1, size_of_col)+1
      my_row =    (myrank-1)/size_of_col +1

      call mpi_comm_split(mpi_comm_world,my_row,my_col,
     &                    mpi_comm_col,ierr)
      call mpi_comm_split(mpi_comm_world,my_col,my_row,
     &                    mpi_comm_row,ierr)


      n1 = max(size_of_col,size_of_row)
      n2 = min(size_of_col,size_of_row)

      do
         if ( n1 == n2 ) then
            n_common = n1
            exit
         end if
         n3 = n1-n2
         n1 = max(n2,n3)
         n2 = min(n2,n3)
      end do!!

      p0_ = -1
      q0_ = -1
      do i0=1,size_of_col
         if ( mod(i0-1,n_common) == mod(my_row-1,n_common) ) then
            n1 = my_row-i0
            if ( n1 >= 0 ) then
               do i1=1,size_of_col
                  k = +n1+(i1-1)*size_of_row
                  if ( mod(k,size_of_col) == 0 ) then
                     p0_(i0) = k/size_of_col
                     q0_(i0) = (i1-1)
                     exit
                  end if
               end do! i1
            else
               do i1=1,size_of_row
                  k = -n1+(i1-1)*size_of_col
                  if ( mod(k,size_of_row) == 0 ) then
                     q0_(i0) = k/size_of_row
                     p0_(i0) = (i1-1)
                     exit
                  end if
               end do! i1
            end if
         end if
      end do! i0
      p0_ = p0_+1
      q0_ = q0_+1


      diag_0 = 0
      diag_1 = 0
      do i0=1,size_of_row/n_common
         j0 = (i0-1)*size_of_row+my_row
         k = mod(j0-1,size_of_col)+1
         if ( k == my_col ) then
            diag_0 = i0
            diag_1 = (j0-1)/size_of_col+1
            exit
         end if
      end do! i0

!---- call mpi_barrier(mpi_comm_world, ierr)

      return
      end subroutine ! eigen_init
!----
      subroutine eigen_free(flag)
      implicit none

      include 'trd.h'

      integer                :: flag, ierr

      real(8) ::    time_bcast, time_reduce, time_redist, time_gather
      common /stat/ time_bcast, time_reduce, time_redist, time_gather

      call mpi_comm_free(mpi_comm_col,ierr)
      call mpi_comm_free(mpi_comm_row,ierr)

#ifdef DETAIL
      if ( flag == 1 .and. myrank == 1 ) then
         print*, "  communication time in \"eigen_prd\""
      endif
      if ( flag == 2 .and. myrank == 1 ) then
         print*, "  communication time in \"eigen_dcx\""
      endif
      if ( flag == 3 .and. myrank == 1 ) then
         print*, " "
         print*, "detail of exectime in \"eigen_pbk\""
         print*, "  communication time in \"eigen_pbk\""
      endif
      if ( flag >= 1 .and. flag <=3 .and. myrank == 1 ) then
         print*, "   bcast  :: ", time_bcast,"(sec)"
         print*, "   reduce :: ", time_reduce,"(sec)"
         print*, "   redist :: ", time_redist,"(sec)"
         print*, "   gather :: ", time_gather,"(sec)"
      endif
#endif

      return
      end subroutine ! eigen_free
!----
      subroutine datacast_dbl(u_y, u_x, u_t, u_s, n)
      implicit none

      integer, intent(in)    :: n
      real(8), intent(inout) :: u_y(1:n), u_x(1:n), u_t(1:n), u_s(1:n)

      include 'trd.h'

      integer :: nx, ny, ic, i0, i1, j0, j1, k0, i, j
      integer :: req(1024), reqr(2), x_snod, y_snod
      integer :: his_rank, her_rank

      real(8) :: d1,d2

      real(8) ::    time_bcast, time_reduce, time_redist, time_gather
      common /stat/ time_bcast, time_reduce, time_redist, time_gather

      d1=mpi_wtime()

      if ( size_of_col == 1 ) then
         if ( size_of_row == 1 ) then
            u_y(1:n) = u_x(1:n)
         else
            ny = (n-1)/size_of_row+1
            do i0=1,ny
               j0 = my_row+size_of_row*(i0-1)
               u_y(i0) = u_x(j0)
            end do! i0
         end if
         return
      end if

      if ( size_of_col == size_of_row ) then
         if ( my_col == my_row ) then
            u_y(1:n) = u_x(1:n)
         end if
         call bcast_dbl(u_y, n, my_row, mpi_comm_col)
         return
      end if

      x_snod = size_of_col/n_common
      y_snod = size_of_row/n_common

      if ( p0_(my_col) > 0 ) then

         nx = (n-p0_(my_col))/y_snod+1
         do i0=1,nx
            j0 = p0_(my_col)+y_snod*(i0-1)
            k0 = q0_(my_col)+x_snod*(i0-1)
            u_t(i0) = u_x(j0)
            u_y(k0) = u_x(j0)
         end do! i0

         do i0=1,x_snod-1
!----
!---- receiving message : length    = ifloor(n/y_snod)
!----                   : #sender   = size_of_col-1
!----
            his_rank = mod(my_col-1 +size_of_col +i0*n_common,
     &                     size_of_col)+1
            ny = (n-p0_(his_rank))/y_snod+1
            call irecv_dbl(u_s, ny, his_rank, reqr(1), mpi_comm_col)

!----
!---- sending message   : length    = ifloor(n/y_snod)
!----                   : #receiver = size_of_col-1
!----
            her_rank = mod(my_col-1 +size_of_col -i0*n_common,
     &                     size_of_col)+1
            call isend_dbl(u_t, nx, her_rank, req(i0), mpi_comm_col)

            call wait_dbl(reqr(1))

            do i1=1,ny
               j1 = q0_(his_rank)+x_snod*(i1-1)
               u_y(j1) = u_s(i1)
            end do! i1

         end do! i0
         call waitall_dbl(x_snod-1, req)

         do i0=1,n_common-1
            her_rank = mod(my_col-1 +size_of_col +i0,size_of_col)+1
            call isend_dbl(u_y, n, her_rank, req(i0), mpi_comm_col)
         end do! i0
         call waitall_dbl(n_common-1, req)

      else

         i = mod(my_row-1,n_common)
         j = mod(my_col-1,n_common)
         ic = mod(j-i+n_common,n_common)
         his_rank = mod(my_col-1 +size_of_col -ic,size_of_col)+1
         call recv_dbl(u_y, n, his_rank, mpi_comm_col)

      end if

9999  continue

      d2=mpi_wtime()
      time_redist=time_redist+(d2-d1)

      return
      end subroutine ! datacast_dbl
!----
      subroutine datacast_dbl2(ur_y,ui_y, ur_x,ui_x, u_t,u_s, n)
      implicit none

      integer, intent(in)    :: n
      real(8), intent(inout) :: ur_y(1:n), ui_y(1:n)
      real(8), intent(inout) :: ur_x(1:n), ui_x(1:n)
      real(8), intent(inout) :: u_t(1:2*n), u_s(1:2*n)

      include 'trd.h'

      integer :: nx, ny, ic, i0, i1, j0, k0, k1, i, j
      integer :: req(1024), reqr(2), x_snod, y_snod
      integer :: his_rank, her_rank

      real(8) :: d1,d2
      real(8) ::    time_bcast, time_reduce, time_redist, time_gather
      common /stat/ time_bcast, time_reduce, time_redist, time_gather


      d1=mpi_wtime()

      if ( size_of_col == 1 ) then
         if ( size_of_row == 1 ) then
            ur_y(1:n)=ur_x(1:n)
            ui_y(1:n)=ui_x(1:n)
         else
            ny=(n-1)/size_of_row+1
            do i0=1,ny
               j0=my_row+size_of_row*(i0-1)
               ur_y(i0)=ur_x(j0)
               ui_y(i0)=ui_x(j0)
            end do! i0
         end if
         return
      end if

      if ( size_of_col == size_of_row ) then
         if ( my_col == my_row ) then
            u_t(  1:  n) = ur_x(1:n)
            u_t(n+1:n+n) = ui_x(1:n)
         end if
         call bcast_dbl(u_t, 2*n, my_row, mpi_comm_col)
         ur_y(1:n) = u_t(  1:  n)
         ui_y(1:n) = u_t(n+1:n+n)
         return
      end if

      x_snod = size_of_col/n_common
      y_snod = size_of_row/n_common

      if ( p0_(my_col) > 0 ) then

         nx = (n-p0_(my_col))/y_snod+1
         do i0=1,nx
            j0 = p0_(my_col)+y_snod*(i0-1)
            k0 = q0_(my_col)+x_snod*(i0-1)
            u_t(   i0) = ur_x(j0)
            u_t(nx+i0) = ui_x(j0)
            ur_y(k0) = ur_x(j0)
            ui_y(k0) = ui_x(j0)
         end do! i0

         do i0=1,x_snod-1
!----
!---- receiving message : length    = ifloor(n/y_snod)
!----                   : #sender   = size_of_col-1
!----
            his_rank = mod(my_col-1 +size_of_col +i0*n_common,
     &                     size_of_col)+1
            ny = (n-p0_(his_rank))/y_snod+1
            call irecv_dbl(u_s, 2*ny,
     &                     his_rank, reqr(1), mpi_comm_col)

!----
!---- sending message   : length    = ifloor(n/y_snod)
!----                   : #receiver = size_of_col-1
!----
            her_rank = mod(my_col-1 +size_of_col -i0*n_common,
     &                     size_of_col)+1
            call isend_dbl(u_t, 2*nx,
     &                     her_rank, req(i0), mpi_comm_col)

            call wait_dbl(reqr(1))

            do i1=1,ny
               k1 = q0_(his_rank)+x_snod*(i1-1)
               ur_y(k1) = u_s(   i1)
               ui_y(k1) = u_s(ny+i1)
            end do! i

         end do! i0
         call waitall_dbl(x_snod-1, req)

         do i0=1,n_common-1
            her_rank = mod(my_col-1 +size_of_col +i0,size_of_col)+1
            u_t(  1:  n) = ur_y(1:n)
            u_t(n+1:n+n) = ui_y(1:n)
            call isend_dbl(u_t, 2*n,
     &                     her_rank, req(i0), mpi_comm_col)
         end do! i0
         call waitall_dbl(n_common-1, req)

      else

         i = mod(my_row-1,n_common)
         j = mod(my_col-1,n_common)
         ic = mod(j-i+n_common,n_common)
         his_rank = mod(my_col-1 +size_of_col -ic,size_of_col)+1
         call recv_dbl(u_s, 2*n,
     &                 his_rank, mpi_comm_col)
         ur_y(1:n) = u_s(  1:  n)
         ui_y(1:n) = u_s(n+1:n+n)

      end if

9999  continue

      d2=mpi_wtime()
      time_redist=time_redist+(d2-d1)

      return
      end subroutine ! datacast_dbl2

      end module communication_sx
