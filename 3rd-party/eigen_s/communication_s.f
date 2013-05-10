      module communication_s

      implicit none

      include 'mpif.h'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  get_loop_start  : calculate a local loop start position from the
!                    global loop start position.
!  get_loop_end    : calculate a local loop end position from the
!                    global loop end position.
!  get_owner_node  : calculate number of the specified rank.
!  translate_g2l   : calculate a local array index from a global array
!                    index.
!  translate_l2g   : calculate a global array index from a local array
!                    index.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      public  :: get_loop_start
      public  :: get_loop_end
      public  :: get_owner_node
      public  :: translate_l2g
      public  :: translate_g2l

      public  :: bcast_dbl
      public  :: reduce_dbl
      public  :: allgather_dbl
      public  :: eigen_init
      public  :: eigen_free
      public  :: datacast_dbl

      private :: isend_dbl
      private :: recv_dbl
      private :: irecv_dbl
      private :: wait_dbl
      private :: waitall_dbl
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

      real(8) ::    d1, d2
      real(8) ::    time_bcast, time_reduce, time_redist, time_gather
      common /stat/ time_bcast, time_reduce, time_redist, time_gather

#ifdef TIMER
      call mpi_barrier(icom,ierr)
      d1=mpi_wtime()
#endif

      call mpi_bcast(buf, n, mpi_double_precision,
     &               iroot-1, icom, ierr)

!---- call d_bbcast(buf, n, iroot-1, icom, ierr)

#ifdef TIMER
      d2=mpi_wtime()
      time_bcast=time_bcast+(d2-d1)
#endif

      return
      end subroutine ! bcast_dbl
!----
      subroutine reduce_dbl(buf, wrk, n, dist, icom)
      implicit none

      integer, intent(in)    :: n, dist, icom
      real(8), intent(inout) :: buf(1:n), wrk(1:n)

      real(8) ::    d1, d2
      real(8) ::    time_bcast, time_reduce, time_redist, time_gather
      common /stat/ time_bcast, time_reduce, time_redist, time_gather

#ifdef TIMER
      call mpi_barrier(icom,ierr)
      d1=mpi_wtime()
#endif

      call mpi_allreduce(buf, wrk, n, mpi_double_precision,
     &                   mpi_sum, icom, ierr)
      buf(1:n) = wrk(1:n)

#ifdef TIMER
      d2=mpi_wtime()
      time_reduce=time_reduce+(d2-d1)
#endif

      return
      end subroutine ! reduce_dbl
!----
      subroutine allgather_dbl(buf, wrk, n, icom)
      implicit none

      integer, intent(in)    :: n, icom
      real(8), intent(inout) :: buf(1:n), wrk(1:n)

      real(8) ::    d1, d2
      real(8) ::    time_bcast, time_reduce, time_redist, time_gather
      common /stat/ time_bcast, time_reduce, time_redist, time_gather

#ifdef TIMER
      call mpi_barrier(icom,ierr)
      d1=mpi_wtime()
#endif

      call mpi_allgather(buf, n, mpi_double_precision,
     &                   wrk, n, mpi_double_precision,
     &                   icom, ierr)

#ifdef TIMER
      d2=mpi_wtime()
      time_gather=time_gather+(d2-d1)
#endif

      return
      end subroutine ! allgather_dbl
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
      end function ! translate_l2g
!----
      integer function get_owner_node(ictr, nnod, inod)
      implicit none

      integer, intent(in)    :: ictr, nnod, inod


      get_owner_node   = mod(ictr-1,nnod)+1


      return
      end function ! get_owner_node
!----
      integer function translate_g2l(ictr, nnod, inod)
      implicit none

      integer, intent(in)    :: ictr, nnod, inod


      translate_g2l   = (ictr-1)/nnod+1


      return
      end function ! translate_g2l
!----
      subroutine eigen_init(ndim)
      implicit none

      include 'trd.h'

      integer ::  ndim
      integer ::  n1, n2, n3, i, j, k

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
!---- print*,myrank,"pe partition = ",size_of_col,size_of_row
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
      do i=1,size_of_col
         if ( mod(i-1,n_common) == mod(my_row-1,n_common) ) then
            n1 = my_row-i
            if ( n1 >= 0 ) then
               do j=1,size_of_col
                  k = +n1+(j-1)*size_of_row
                  if ( mod(k,size_of_col) == 0 ) then
                     p0_(i) = k/size_of_col
                     q0_(i) = (j-1)
                     exit
                  end if
               end do! j
            else
               do j=1,size_of_row
                  k = -n1+(j-1)*size_of_col
                  if ( mod(k,size_of_row) == 0 ) then
                     q0_(i) = k/size_of_row
                     p0_(i) = (j-1)
                     exit
                  end if
               end do! j
            end if
         end if
      end do! i
      p0_ = p0_+1
      q0_ = q0_+1


      diag_0 = 0
      diag_1 = 0
      do i=1,size_of_row/n_common
         j = (i-1)*size_of_row+my_row
         k = mod(j-1,size_of_col)+1
         if ( k == my_col ) then
            diag_0 = i
            diag_1 = (j-1)/size_of_col+1
            exit
         end if
      end do! i_1

!---- call mpi_barrier(mpi_comm_world, ierr)

      return
      end subroutine ! eigen_init
!----
      subroutine eigen_free(flag)
      implicit none

      include 'trd.h'

      integer :: flag

      real(8) ::    time_bcast, time_reduce, time_redist, time_gather
      common /stat/ time_bcast, time_reduce, time_redist, time_gather


      call mpi_comm_free(mpi_comm_col,ierr)
      call mpi_comm_free(mpi_comm_row,ierr)

#ifdef DETAIL
      if ( flag == 1 .and. myrank == 1 ) then
         print*, "  communication time in \"eigen_trd\""
      endif
      if ( flag == 2 .and. myrank == 1 ) then
         print*, " "
         print*, "detail of exectime in \"eigen_dc\""
         print*, "  communication time in \"eigen_dc\""
      endif
      if ( flag == 3 .and. myrank == 1 ) then
         print*, "  communication time in \"eigen_tbk\""
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

      integer :: nx, ny, ic, i, j, k
      integer :: req(1024), reqr(2), x_snod, y_snod
      integer :: his_rank, her_rank

      real(8) ::    d1, d2
      real(8) ::    time_bcast, time_reduce, time_redist, time_gather
      common /stat/ time_bcast, time_reduce, time_redist, time_gather

#ifdef TIMER
      call mpi_barrier(mpi_comm_col,ierr)
      d1=mpi_wtime()
#endif

      if ( size_of_col == 1 ) then
         if ( size_of_row == 1 ) then
            u_y(1:n) = u_x(1:n)
         else
            ny = (n-1)/size_of_row+1
            do i=1,ny
               j = my_row+size_of_row*(i-1)
               u_y(i) = u_x(j)
            end do! i
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
         do i=1,nx
            j = p0_(my_col)+y_snod*(i-1)
            k = q0_(my_col)+x_snod*(i-1)
            u_t(i) = u_x(j)
            u_y(k) = u_x(j)
!----       print*,"i am",myrank,"sdrv ",myrank
         end do! i

         do ic=1,x_snod-1

!----
!----       receiving message : length    = ifloor(n/y_snod)
!----                         : #sender   = size_of_col-1
!----
            his_rank =
     &            mod(my_col-1 +size_of_col +ic*n_common,size_of_col)+1
            ny = (n-p0_(his_rank))/y_snod+1
!----       print*,"i am",myrank,"recv ",(his_rank-1)*size_of_row+my_row
            call irecv_dbl(u_s, ny, his_rank, reqr(1), mpi_comm_col)

!----
!----       sending message   : length    = ifloor(n/y_snod)
!----                         : #receiver = size_of_col-1
!----
            her_rank =
     &            mod(my_col-1 +size_of_col -ic*n_common,size_of_col)+1
            call isend_dbl(u_t, nx, her_rank, req(ic), mpi_comm_col)
!----       print*,"i am",myrank,"send ",(her_rank-1)*size_of_row+my_row

            call wait_dbl(reqr(1))

            do i=1,ny
               k = q0_(his_rank)+x_snod*(i-1)
               u_y(k) = u_s(i)
            end do! i

         end do! ic
         call waitall_dbl(x_snod-1, req)

         do ic=1,n_common-1
            her_rank = mod(my_col-1 +size_of_col +ic,size_of_col)+1
!----       print*,"i am",myrank,"send ",(her_rank-1)*size_of_row+my_row
            call isend_dbl(u_y, n, her_rank, req(ic), mpi_comm_col)
!----       call wait_dbl(req(ic))
         end do! ic
         call waitall_dbl(n_common-1, req)

      else

         i = mod(my_row-1,n_common)
         j = mod(my_col-1,n_common)
         ic = mod(j-i+n_common,n_common)
         his_rank = mod(my_col-1 +size_of_col -ic,size_of_col)+1
!----    print*,"i am",myrank,"recv ",(his_rank-1)*size_of_row+my_row
         call recv_dbl(u_y, n, his_rank, mpi_comm_col)

      end if

#ifdef TIMER
      d2=mpi_wtime()
      time_redist=time_redist+(d2-d1)
#endif

      return
      end subroutine ! datacast_dbl

      end module communication_s
