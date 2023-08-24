      subroutine mat_set(n, a, nm)

      use eigen_libs_mod
      use mpi
!$    use omp_lib

      implicit none

      integer, intent(in)    :: n, nm
      real(8), intent(out)   :: a(1:nm,*)

      real(8), parameter     :: ZERO= 0.0D0
      real(8), parameter     :: ONE = 1.0D0

      integer                :: COMM, x_COMM, y_COMM
      integer                :: nnod, x_nnod, y_nnod
      integer                :: inod, x_inod, y_inod
      integer                :: iloop_sta, iloop_end
      integer                :: jloop_sta, jloop_end
      integer                :: i, i_1
      integer                :: j, j_1

      call eigen_get_comm (COMM, x_COMM, y_COMM)
      call eigen_get_procs(nnod, x_nnod, y_nnod)
      call eigen_get_id   (inod, x_inod, y_inod)

      jloop_sta = eigen_loop_start(1, x_nnod, x_inod)
      jloop_end = eigen_loop_end  (n, x_nnod, x_inod)
      iloop_sta = eigen_loop_start(1, y_nnod, y_inod)
      iloop_end = eigen_loop_end  (n, y_nnod, y_inod)

      if (iloop_sta <= iloop_end) then
!$OMP PARALLEL DO &
!$OMP& PRIVATE(i,j,i_1,j_1)
         do i_1 = iloop_sta, iloop_end
            i = eigen_translate_l2g(i_1, y_nnod, y_inod)
            do j_1 = jloop_sta, jloop_end
               j = eigen_translate_l2g(j_1, x_nnod, x_inod)
               a(j_1, i_1) = min(i,j)*ONE
            end do
         end do
!$OMP END PARALLEL DO
      end if

      end subroutine  mat_set


program EigenExa
  use eigen_libs_mod
  use mpi
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

  call mpi_barrier(mpi_comm_world, ierr)
  init_tick = mpi_wtime()
  call mpi_comm_rank( mpi_comm_world, i_inod, ierr )
  call eigen_init()
  call mpi_barrier(mpi_comm_world, ierr)
  gen_tick = mpi_wtime()
  call eigen_get_matdims( n, nm, ny )
  allocate( a(nm, ny), z(nm, ny), w(n))
  call mat_set( n, a, nm )

  call mpi_barrier(mpi_comm_world, ierr)
  diag_tick = mpi_wtime()
  call eigen_sx( n, n, a, nm, w, z, nm, 48, 128, 'A')
  call mpi_barrier(mpi_comm_world, ierr)
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
