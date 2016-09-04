program EigenExa
  use MPI
  use eigen_libs
  implicit none
  double precision   init_tick, gen_tick, diag_tick, end_tick
  integer :: n, m, nm, ny, mtype, i, i_inod, ierr
  real(8), pointer :: a(:,:),z(:,:),w(:)
  character(len=10) :: tmp_str
  integer arg_len, status

  call mpi_init_thread( MPI_THREAD_MULTIPLE, i, ierr )
  if (command_argument_count().eq.1) then
     call get_command_argument(1, tmp_str, arg_len, status)
     read(tmp_str, *) n
  else
     write(*,'(A)') "Error: eigen_exa dimension"
     stop
  endif
  
  init_tick = mpi_wtime()
  call mpi_comm_rank( mpi_comm_world, i_inod, ierr )
  call eigen_init()
  gen_tick = mpi_wtime()
  call eigen_get_matdims( n, nm, ny )
  allocate( a(nm, ny), z(nm, ny), w(n))
  call mat_set( n, a(1,1), nm, 3 ) ! minij matrix

  diag_tick = mpi_wtime()
  call eigen_sx( n, n, a, nm, w, z, nm, 48, 128, 'A')
  end_tick = mpi_wtime()
  if ( i_inod == 0 ) then
     print *, "Matrix dimension = ", n
     do i = 1, min(100, n)
        print *, i, w(i)
     end do
     print *, "init_time = ", gen_tick - init_tick
     print *, "gen_time = ", diag_tick - gen_tick
     print *, "diag_time = ", end_tick - diag_tick
  end if

  call eigen_free( )
  call MPI_Finalize( ierr )
end program EigenExa
