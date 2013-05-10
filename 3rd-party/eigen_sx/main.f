!----
      use communication_sx, only : eigen_init
     &                           , eigen_free
     &                           , translate_l2g
!----
      implicit double precision (a-h,o-z)

      real(8),allocatable :: a(:,:)
      real(8), pointer ::
     &           b(:),z(:),d(:),e(:),w(:)
      integer, parameter :: NDIM = 2
!
      include 'mpif.h'
      include 'trd.h'
!
!
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,i$inod,ierr)
      call mpi_comm_size(mpi_comm_world,i$nnod,ierr)
*
      n=1000
      m=32
      allocate(a(n,n))
      call eigen_init(NDIM)
      NPROW = size_of_col
      NPCOL = size_of_row
      nx = ((n-1)/NPROW+1)
      call CSTAB_get_optdim(nx, 6, 16*2, 16*4, nm)
      call eigen_free(0)

      NB  = 64+32
      nmz = ((n-1)/NPROW+1)
      nmz = ((nmz-1)/NB+1)*NB+1
      nmw = ((n-1)/NPCOL+1)
      nmw = ((nmw-1)/NB+1)*NB+1
      nme = ((n-1)/2+1)*2

      nh = (n-1)/4+1
      call CSTAB_get_optdim(nh, 4, 16*2, 16*4, nnh)
      nnh = 4 * nnh

      n1x = ((n-1)/nprocs+1)
      larray = MAX(MAX(nmz,nm,nnh)*nmw, n*n1x)
      allocate(
     &         b(larray), z(larray),
     &         d(nme), e(2*nme), w(nme),
     &         stat=istat)
      if(istat.ne.0) then
         print*,"Memory exhausted"
         call flush(6)
         stop
      endif
*
      call matrix_set(n,a)
      call matrix_adjust_sx(n,a,b,nm)
*-
      call eigen_sx(n,b,nm,w,z,d,e,nme,m,0)
*-
      call matrix_set(n,a)
      call matrix_adjust_sx(n,a,b,nm)
      call ev_test_2D(n, b, nm, w, z, nm, t)
*
      deallocate(a)
      deallocate(z)
      deallocate(d)
      deallocate(e)
      deallocate(w)
*
      call MPI_Finalize(ierr)
      end

