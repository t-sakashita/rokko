program EigenExa
  use eigen_libs
  use mpi
  implicit none
  integer :: n, m, nm, ny, mtype, i, i_inod, ierr
  real(8), pointer :: a(:,:),z(:,:),w(:)
  call mpi_init_thread( MPI_THREAD_MULTIPLE, i, ierr )
  call mpi_comm_rank( mpi_comm_world, i_inod, ierr )
  call eigen_init( )

  n = 100
  call eigen_get_matdims( n, nm, ny )
  allocate( a(nm, ny), z(nm, ny), w(n))
  call mat_set( n, a(1,1), nm, 0 ) ! Frank matrix

  call eigen_sx( n, n, a, nm, w, z, nm, 48, 128, 'A')

  if ( i_inod == 0 ) then
     print*, "Matrix dimension = ", N
     do i = 1, n
        print *, i, w(i)
     end do
  end if

  call eigen_free( )
  call MPI_Finalize( ierr )
end program EigenExa
