!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  n    :  dimension of matrix
!  d    :  array to store a main diagonal element
!  e    :  array to store a sub diagonal element
!  nme  :  size of "e" array
!  z    :  array for eigen vector
!  nmz  :  size of "z" array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine eigen_dcx(n, d, e, nme, z, nmz)
!----
      use communication_sx, only : eigen_init
     &                           , eigen_free
     &                           , bcast_dbl   
!----
      integer, intent(in)           :: n, nmz
      real(8), intent(inout)        :: d(1:n), e(1:2*nme)
      real(8), intent(out)          :: z(1:nmz,*)

! parameters blacs array descritor(the position of entry tags), etc
      integer, parameter            :: block_cyclic_2d = 1
      integer, parameter            :: dlen_  = 9
      integer, parameter            :: dtype_ = 1
      integer, parameter            :: ctxt_  = 2
      integer, parameter            :: m_     = 3
      integer, parameter            :: n_     = 4
      integer, parameter            :: mb_    = 5
      integer, parameter            :: nb_    = 6
      integer, parameter            :: rsrc_  = 7
      integer, parameter            :: csrc_  = 8
      integer, parameter            :: lld_   = 9

      integer                       :: descz( dlen_ )
      integer                       :: descw( dlen_ )
 
      integer                       :: my_rank
      integer                       :: lwork, liwork
 
      real(8), pointer              :: work(:)
      integer, pointer              :: iwork(:)

      real(8)  :: d1,d2,d3,d4,d5,d6

      include 'mpif.h'
      include 'trd.h'
 
 
      d1 = mpi_wtime()
! blacs/pblas/scalapack initialization
      call blacs_pinfo( my_rank, nprocs )
      if ( nprocs < 1 ) then
!  mpi group setup
         call mpi_comm_size( mpi_comm_world, nprocs, ierr )
         call mpi_comm_rank( mpi_comm_world, my_rank, ierr )
         call blacs_setup( my_rank, nprocs )
      end if


      call blacs_get( -1, 0, ictxt )
   
      call eigen_init(2)

      nprow = size_of_col
      npcol = size_of_row

      call blacs_gridinit( ictxt, 'column-major', nprow, npcol )
      call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )

      d2 = mpi_wtime()

! blacs array registration
      nb = 64+32                    ! size of block width
      np = numroc( n, nb, myrow, 0, nprow )
      nq = numroc( n, nb, mycol, 0, npcol )
      lddz = (n-1)/nprow+1
      lddz = ((lddz-1)/nb+1)*nb     ! first dimension size of array "z"
      lddw = (n-1)/npcol+1
      lddw = ((lddw-1)/nb+1)*nb     ! second dimension size of array "z"
      call descinit( descz, n, n, nb, nb, 0, 0, ictxt, lddz, info )
      if(info/=0) then
         print*,"Error in 'descinit' routine."
         print*,"  error code =",info
         stop 
      endif
!----
! preparing working arrays
      nx     = (n-1)/npcol+1
      lwork  = max(1+6*n+2*np*(nq+max(nq,nb)), lddz*lddw, nmz*nx)
      liwork = 2+7*n+8*npcol
      allocate(work(lwork), iwork(liwork), stat=istat)
      if(istat /= 0) then
         print*,"memory exhausted"
         call flush(6)
         call mpi_abort( mpi_comm_world, 1, ierr )
      end if

! somehow, z must be nullified (originally next loop is not required.)
      z(1:lddz*lddw,1) = 0.0d+00

c #ifdef DETAIL
c       d3 = mpi_wtime()
c       if(my_rank==0) then
c          print*," "
c          print*,"detail of exectime in \"eigen_dcx\""
c          print*,"   before pdstedc =",d3-d1,"(sec)"
c       endif
c #endif

      d4 = mpi_wtime ()
      call eigen_pdsxedc('i',2, n, d(1), e(1), nme, z(1,1), 1, 1, descz,
     &               work(1), lwork, iwork(1), liwork, info)
      d5 = mpi_wtime()
      if(info/=0) then
         print*,"Error in 'eigen_pdsxedc' routine."
         print*,"  error code =",info
         stop 
      endif
c #ifdef DETAIL
c       if(my_rank==0) print*,"   pdstedc =",d5-d4,"(sec)"
c #endif

      if(nb==1)then
         do i0=nx,1,-1
            work(1:lddz)=z(1+(i0-1)*lddz:lddz+(i0-1)*lddz,1)
            z(1:lddz,i0)=work(1:lddz)
         enddo
         do i0=1,nx
            z(lddz+1:nmz,i0)=0.0d0
         enddo
      else
         call descinit( descw, n, n, 1, 1, 0, 0, ictxt, nmz, info )
         call pdgemr2d( n, n, z, 1, 1, descz, work, 1, 1, descw, ictxt )
         z(1:nmz*nx,1)=work(1:nmz*nx)
      endif

      call bcast_dbl( d(1), n, 1, mpi_comm_world )
      d6 = mpi_wtime()
      call eigen_free(2)

c #ifdef TIMER
c       if(my_rank==0) then
c          print*,"Exectime of \"eigen_dcx\" routine  =",d6-d1,"(sec)"
c       endif
c #endif

! freeing working arrays
      deallocate(work)
      deallocate(iwork)

! blacs/pblas/scalapack finalize
      call blacs_gridexit( ictxt )

      return
      end subroutine

