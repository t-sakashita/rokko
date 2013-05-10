      use communication_s, only : eigen_init
     &                          , eigen_free
      implicit double precision (a-h,o-z)

      real(8),allocatable :: a(:,:)
      real(8), pointer :: b(:),z(:),w(:)
!
      include 'mpif.h'
      include 'trd.h'
!
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,i$inod,ierr)
      call mpi_comm_size(mpi_comm_world,i$nnod,ierr)
*
      n=1000
      m=32
      allocate(a(n,n))
      call eigen_init(2)
      NPROW = size_of_col
      NPCOL = size_of_row
      nx = ((n-1)/NPROW+1)
      call CSTAB_get_optdim(nx, 6, 16*4, 16*4*2, nm)
      call eigen_free(0)

      NB  = 64+32
      nmz = ((n-1)/NPROW+1)
      nmz = ((nmz-1)/NB+1)*NB+1
      nmw = ((n-1)/NPCOL+1)
      nmw = ((nmw-1)/NB+1)*NB+1

      larray = MAX(nmz,nm)*nmw
      allocate( b(larray), z(larray), w(n), stat=istat)
      if(istat.ne.0) then
         print*,"Memory exhausted"
         call flush(6)
         stop
      endif
*
      call matrix_set(n, a)
      call matrix_adjust_s(n, a, b, nm)
      call eigen_s(n, b, nm, w(1), z(1), nm, m, 0)
*-
      call matrix_set(n, a)
      call matrix_adjust_s(n, a, b, nm)
      call ev_test_2D(n, b, nm, w(1), z(1), nm)

      deallocate(a)
      deallocate(b)
      deallocate(z)
      deallocate(w)
*
      call MPI_Finalize(ierr)
      end

