MODULE laplacian
  use MPI
  use rokko
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  integer :: comm
  integer :: nprocs, myrank
  integer(c_int), private :: dim, local_offset, num_local_rows
  integer(c_int), private :: start_row, end_row
  integer(c_int), private :: start_k, end_k
  logical, private :: is_first_proc, is_last_proc
  double precision, private :: buf_p = 0, buf_m = 0
  INTEGER, private :: status_m(MPI_STATUS_SIZE), status_p(MPI_STATUS_SIZE)
CONTAINS
  SUBROUTINE initialize (mat, dim_in)
    type(rokko_distributed_mfree), intent(inout) :: mat
    integer(c_int), intent(in) :: dim_in
    integer :: ierr
    integer :: tmp, rem
    comm = MPI_COMM_WORLD
    call MPI_comm_rank(comm, myrank, ierr)
    call MPI_comm_size(comm, nprocs, ierr)
    dim = dim_in
    tmp = dim / nprocs
    rem = mod(dim, nprocs)
    num_local_rows = (dim + nprocs - myrank - 1) / nprocs
    start_row = tmp * myrank + min(rem, myrank) + 1  ! extra plus 1
    end_row = start_row + num_local_rows - 1 ! extra plus 1
    
    if (start_row == 1) then
       is_first_proc = .true.
    else
       is_first_proc = .false.
    endif
 
    if (end_row == dim) then
       is_last_proc = .true.
    else
       is_last_proc = .false.
    endif

    end_k = num_local_rows
    print*, "myrank=", myrank, "start_row=", start_row, "end_row=", end_row,&
         "is_first_proc", is_first_proc, " is_last_proc", is_last_proc
    print*, "myrank=", myrank, "num_local_rows=", num_local_rows
    call rokko_distributed_mfree_construct(mat, multiply, dim, num_local_rows)
  END SUBROUTINE initialize
  ! subroutine passed to C function.
  ! It must be interoperable!
  SUBROUTINE multiply (n, x, y) BIND(C)
    INTEGER(C_INT), INTENT(IN), VALUE :: n
    DOUBLE PRECISION, INTENT(IN) :: x(n)
    DOUBLE PRECISION, INTENT(OUT) :: y(n)
    integer :: ierr
    integer :: k

    if (num_local_rows == 0) then
       return
    endif
    
    if (.not.(is_first_proc) .and. (nprocs /= 1)) then
       print*, "recv myrank=", myrank
       call MPI_Send(x(1), 1, MPI_DOUBLE_PRECISION, myrank-1, 0, comm, ierr)
       call MPI_Recv(buf_m, 1, MPI_DOUBLE_PRECISION, myrank-1, 0, comm, status_m, ierr)
       !std::cout << "buffff=" << buf << std::endl
    endif
    
    if (.not.(is_last_proc) .and. (nprocs /= 1)) then
       print*, "send myrank=", myrank
       call MPI_Recv(buf_p, 1, MPI_DOUBLE_PRECISION, myrank+1, 0, comm, status_p, ierr)
       call MPI_Send(x(end_k), 1, MPI_DOUBLE_PRECISION, myrank+1, 0, comm, ierr)
       !std::cout << "buffff=" << buf2 << std::endl
    endif
    
    if (is_first_proc) then
       if (num_local_rows /= 1) then
          y(2) = x(2) - x(1)
          if (nprocs /= 1) then
             y(end_k) = - x(end_k - 1) + 2 * x(end_k) - buf_p
          endif
       else 
          y(2) = x(2) - buf_p
       endif
    endif
    
    if (is_last_proc) then
       if (num_local_rows /= 1) then
          if (nprocs /= 1) then
             y(1) = - buf_m + 2 * x(1) - x(2)
          endif
          y(end_k) = 2 * x(end_k) - x(end_k - 1)      
       else
          y(end_k) = 2 * x(end_k) - buf_m
       endif
    endif
    if (.not.(is_first_proc .or. is_last_proc)) then ! neither first or last process
       if (num_local_rows /= 1) then
          y(1) = - buf_m + 2 * x(1) - x(2)
          y(end_k) = - x(end_k - 1) + 2 * x(end_k) - buf_p
       endif
    else
       y(1) = - buf_m + 2 * x(2) - buf_p
    endif
    
    ! from 2 to end
    do k=2, end_k-1
       y(k) = - x(k-1) + 2 * x(k) - x(k+1)
    enddo
  END SUBROUTINE multiply
END MODULE laplacian

program main
  use MPI
  use rokko
  use laplacian
  implicit none
  integer :: provided, ierr
  
  double precision :: eig_val
  double precision, allocatable, dimension(:) :: eig_vec

  type(rokko_parallel_sparse_solver) :: solver
  character(len=100) :: solver_name, tmp_str
  integer arg_len, status
  type(rokko_distributed_mfree) :: mat
  integer :: dim
  integer :: num_evals, block_size, max_iters
  integer :: num_local_rows, num_conv
  double precision :: tol

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(comm, myrank, ierr)
  call MPI_comm_size(comm, nprocs, ierr)

  solver_name = "anasazi"  ! default
  dim = 10
  if (command_argument_count().ge.1) then
     call get_command_argument(1, tmp_str, arg_len, status)
     solver_name = trim(tmp_str)
  endif

  if (command_argument_count().eq.2) then
     call get_command_argument(2, tmp_str, arg_len, status)
     read(tmp_str, *) dim
  endif

  if (myrank == 0) then
     write(*,*) "solver name = ", trim(solver_name)
     write(*,*) "matrix dimension = ", dim
  endif

  call rokko_parallel_sparse_solver_construct(solver, solver_name)
  call initialize(mat, dim)

  num_evals = 2
  block_size = 5
  max_iters = 500
  tol = 1.0e-8
  call rokko_parallel_sparse_solver_diagonalize_distributed_mfree(solver, mat, num_evals, block_size, max_iters, tol)

  num_conv = rokko_parallel_sparse_solver_num_conv(solver)
  eig_val = rokko_parallel_sparse_solver_eigenvalue(solver, 0)
  num_local_rows = rokko_distributed_mfree_num_local_rows(mat)
  print*, "eigval=", eig_val
  call rokko_distributed_mfree_destruct(mat)
  call rokko_parallel_sparse_solver_destruct(solver)

  call MPI_finalize(ierr)
end program main
