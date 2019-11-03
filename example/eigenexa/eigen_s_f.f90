program main
  use eigen_libs_mod
  use mpi
  implicit none
  integer, parameter :: n = 8
  integer :: provided, ierr, myrank, nprocs, nprow, npcol
  integer :: inod, x_inod, y_inod
  integer :: nm, ny, i, j
  integer :: il, jl, iloop_sta, iloop_end, jloop_sta, jloop_end
  double precision, allocatable :: a(:,:), z(:,:), w(:)

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call eigen_init(MPI_COMM_WORLD)
  call eigen_get_procs(nprocs, npcol, nprow)
  if (myrank == 0) then
     print *, "n =", n
     print *, "nprocs =", nprocs
     print *, "nprow =", nprow
     print *, "npcol =", npcol
  end if
  call eigen_get_id(inod, x_inod, y_inod)

  call eigen_get_matdims(n, nm, ny)
  allocate(a(nm, ny), z(nm, ny), w(n))

  jloop_sta = eigen_loop_start(1, nprow, y_inod)
  jloop_end = eigen_loop_end(n, nprow, y_inod)
  iloop_sta = eigen_loop_start(1, npcol, x_inod)
  iloop_end = eigen_loop_end(n, npcol, x_inod)
  if (jloop_sta <= jloop_end) then
     do jl = jloop_sta, jloop_end
        j = eigen_translate_l2g(jl, nprow, y_inod)
        do il = iloop_sta, iloop_end
           i = eigen_translate_l2g(il, npcol, x_inod)
           a(il, jl) = dble(min(i, j))
        end do
     end do
  end if

  call eigen_s(n, n, a, nm, w, z, nm, 48, 128, 'A')
  if (myrank == 0) then
     print *, "eigenvalues:", w(:)
  end if

  deallocate(a, z, w)
  call eigen_free()
  call MPI_Finalize(ierr)
end program main
