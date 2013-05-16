      subroutine eigen_dc(n, d, e, z, ldz, info)
!----
      use communication_s, only : bcast_dbl
     &                          , eigen_init
     &                          , eigen_free
!----
      integer, intent(in)       :: n, ldz
      real(8), intent(inout)    :: d(1:n), e(1:n)
      real(8), intent(out)      :: z(1:ldz,*)

!---- parameters blacs array descritor(the position of entry tags), etc
      integer, parameter        :: block_cyclic_2d = 1
      integer, parameter        :: dlen_  = 9
      integer, parameter        :: dtype_ = 1
      integer, parameter        :: ctxt_  = 2
      integer, parameter        :: m_     = 3
      integer, parameter        :: n_     = 4
      integer, parameter        :: mb_    = 5
      integer, parameter        :: nb_    = 6
      integer, parameter        :: rsrc_  = 7
      integer, parameter        :: csrc_  = 8
      integer, parameter        :: lld_   = 9
      integer                   :: descz( dlen_ )
      integer                   :: descw( dlen_ )
      integer                   :: my_rank, ierr
      integer                   :: lwork, liwork
      real(8), pointer          :: work(:)
      integer, pointer          :: iwork(:)
      real(8)                   :: hs0, hs1

      include 'mpif.h'
      include 'trd.h'
 
 
!---- blacs/pblas/scalapack initialization
      call blacs_pinfo( my_rank, nprocs )
      if ( nprocs < 1 ) then
!---- mpi group setup
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


!---- blacs array registration
      nb = 32                             ! blocking width
      np = numroc( n, nb, myrow, 0, nprow )
      nq = numroc( n, nb, mycol, 0, npcol )
      lddz = (n-1)/nprow+1
      lddz = ((lddz-1)/nb+1)*nb           ! size of array "z"
      lddw = (n-1)/npcol+1
      lddw = ((lddw-1)/nb+1)*nb           ! size of second dimension of
                                          !  array "z".
!----
      call descinit( descz, n, n, nb, nb, 0, 0, ictxt, lddz, ierr )
 
!---- preparing working arrays
      nx     = (n-1)/npcol+1
      lwork  = max(1+6*n+2*np*(nq+max(nq,nb)), lddz*lddw, ldz*nx)
      liwork = 2+7*n+8*npcol
      allocate(work(lwork), iwork(liwork), stat=istat)
      if(istat /= 0) then
         print*,"memory exhausted"
         call flush(6)
         call mpi_abort( mpi_comm_world, 1, ierr )
      end if


!---- somehow, z must be nullified (originally next loop is not
!----& required.)
      z(1:lddz*lddw,1) = 0.0d+00

       hs0=mpi_wtime()
!---- mkl_mode = mkl_get_dynamic()
!---- call mkl_set_dynamic(0)
      call pdstedc( 'i', n, d(1), e(1), z(1,1), 1, 1, descz,
     &              work(1), lwork, iwork(1), liwork, info )
!---- call mkl_set_dynamic(mkl_mode)

      if(nb==1)then
         do i=nx,1,-1
            work(1:lddz)=z(1+(i-1)*lddz:lddz+(i-1)*lddz,1)
            z(1:lddz,i)=work(1:lddz)
         enddo
         do i=1,nx
            z(lddz+1:ldz,i)=0.0d0
         enddo
      else
         call descinit( descw, n, n, 1, 1, 0, 0, ictxt, ldz, ierr )
         call pdgemr2d( n, n, z, 1, 1, descz, work, 1, 1, descw, ictxt )
         z(1:ldz*nx,1)=work(1:ldz*nx)
      endif
       hs1=mpi_wtime()

      call bcast_dbl( d(1), n, 1, mpi_comm_world )
      call eigen_free(2)
c #ifdef TIMER
c       if ( myrank == 1 ) then
c          print *,"Exectime of 'eigen_dc'  routine =",hs1-hs0,"(sec)"
c       endif
c #endif

!---- freeing working arrays
      deallocate(work)
      deallocate(iwork)

 
!---- blacs/pblas/scalapack finalize
      call blacs_gridexit( ictxt )

      return
      end subroutine

