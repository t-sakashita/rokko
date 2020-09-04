program main
  use mpi
  implicit none
  integer, parameter :: n = 8, nb = 1
  integer :: provided, ierr, myrank, nprocs
  integer :: icontxt, npcol, nprow
  integer :: info
  double precision, allocatable :: d(:), e(:), w(:)
  integer :: lwork, liwork
  double precision, allocatable :: work(:)
  integer, allocatable :: iwork(:)
  integer :: il, iu, m, nsplit
  double precision :: vl, vu, abstol
  integer, allocatable :: isplit(:), iblock(:)

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  nprow = int(sqrt(1. * nprocs + 0.5))
  npcol = int(1. * nprocs / nprow)
  if (myrank == 0) then
     print *, "n =", n
     print *, "nprocs =", nprocs
     print *, "nprow =", nprow
     print *, "npcol =", npcol
     if ((nprocs /= nprow * npcol) .or. (mod(n, nprow) /= 0) .or. &
          (mod(n, npcol) /= 0)) then
        print *, "incompatible matrix size and number of processes"
        call MPI_Abort(MPI_COMM_WORLD, 127, ierr)
     end if
  end if
  call BLACS_get(0, 0, icontxt)
  call BLACS_gridinit(icontxt, 'R', nprow, npcol)  ! not used explicitly, but needed for pdstebz. If you remove it, segmentation fault occurs.

  allocate(d(n), e(n))
  call generate_laplacian(n, d, e)

  allocate(w(n))

  abstol = 0d0
  lwork = -1
  liwork = -1
  allocate(work(1), iwork(1))
  allocate(isplit(n))
  allocate(iblock(n))

  call pdstebz(icontxt, 'A', 'E', n, &
       & vl, vu, il, iu, &
       & abstol, d, e, m, nsplit, &
       & w, iblock, isplit, &
       & work, lwork, iwork, liwork, &
       & info)
  lwork = int(work(1))
  liwork = iwork(1)
  deallocate(work, iwork)
  allocate(work(lwork), iwork(liwork))

  call pdstebz(icontxt, 'A', 'E', n, &
       & vl, vu, il, iu, &
       & abstol, d, e, m, nsplit, &
       & w, iblock, isplit, &
       & work, lwork, iwork, liwork, &
       & info)

  call mpi_barrier(mpi_comm_world, ierr)
  if (myrank == 0) then
     print *, "eigenvalues:", w(:)
  end if

  deallocate(isplit, iblock)
  deallocate(d, e, w, work, iwork)
  call BLACS_gridexit(icontxt)
  call MPI_finalize(ierr)

contains

  subroutine generate_laplacian(n, d, e)
    implicit none
    integer, intent(in) :: n
    double precision, intent(out) :: d(:), e(:)
    integer :: i

    d(1) = 1d0;
    do i = 2, n
       d(i) = 2d0
    enddo
    do i = 1, n-1
       e(i) = -1d0
    enddo
  end subroutine generate_laplacian

end program main
