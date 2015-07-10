MODULE mfree
  USE, INTRINSIC :: ISO_C_BINDING
  use rokko!, only : type(rokko_distributed_mfree)
  IMPLICIT NONE
  TYPE(C_FUNPTR) :: cproc

  ! interface for C function.
!  INTERFACE
!     SUBROUTINE rokko_distributed_mfree_construct2(mat, func_multiply, dim, num_local_rows)
!       type(rokko_distributed_mfree), intent(in) :: mat
!       TYPE(C_FUNPTR), INTENT(IN), VALUE :: func_multiply
!       INTEGER(C_INT), INTENT(IN) :: dim, num_local_rows
!     END SUBROUTINE rokko_distributed_mfree_construct2
!  END INTERFACE
  
CONTAINS
  SUBROUTINE simple2() bind(c)
    USE, INTRINSIC :: ISO_C_BINDING
  END SUBROUTINE simple2
  ! Wrapper for C function.
  SUBROUTINE rokko_distributed_mfree_construct_f(mat, simple_in, dim, num_local_rows)
    USE, INTRINSIC :: ISO_C_BINDING
    type(rokko_distributed_mfree), intent(inout) :: mat
    INTEGER(C_INT), INTENT(IN) :: dim, num_local_rows
    INTERFACE
       SUBROUTINE simple_in()
         USE, INTRINSIC :: ISO_C_BINDING
       END SUBROUTINE simple_in
    END INTERFACE
    ! Get C procedure pointer.
    cproc = C_FUNLOC(simple2)
    ! call wrapper in C.
    call rokko_distributed_mfree_construct2(mat, cproc, dim, num_local_rows)
  END SUBROUTINE rokko_distributed_mfree_construct_f
  
END MODULE mfree
 
MODULE laplacian
  use MPI
  use rokko
  use mfree!, only : rokko_distributed_mfree_construct_f
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  INTEGER(C_INT) myrank, nprocs
  TYPE type_vars
     SEQUENCE
     integer :: comm
     integer(c_int) :: nprocs, myrank
     integer(c_int) :: dim, local_offset, num_local_rows
     integer(c_int) :: start_row, end_row
     integer(c_int) :: start_k, end_k
     logical :: is_first_proc, is_last_proc
     double precision :: buf_m, buf_p
     INTEGER :: status_m(MPI_STATUS_SIZE), status_p(MPI_STATUS_SIZE)
  END TYPE type_vars
  TYPE(type_vars) :: vars
CONTAINS
  SUBROUTINE initialize (mat, dim_in)
    type(rokko_distributed_mfree), intent(inout) :: mat
    integer(c_int), intent(in) :: dim_in
    integer :: ierr
    integer :: tmp, rem
    vars%comm = MPI_COMM_WORLD
    call MPI_comm_rank(MPI_COMM_WORLD, vars%myrank, ierr)
    call MPI_comm_size(MPI_COMM_WORLD, vars%nprocs, ierr)
    vars%dim = dim_in
    tmp = vars%dim / vars%nprocs;
    rem = mod(vars%dim, vars%nprocs);
    vars%num_local_rows = (vars%dim + vars%nprocs - vars%myrank - 1) / vars%nprocs
    vars%start_row = tmp * vars%myrank + min(rem, vars%myrank)
    vars%end_row = vars%start_row + vars%num_local_rows - 1
    
    if (vars%start_row == 1) then
       vars%is_first_proc = .true.
    else
       vars%is_first_proc = .false.
    endif
 
    if (vars%end_row == (vars%dim-1)) then
       vars%is_last_proc = .true.
    else
       vars%is_last_proc = .false.
    endif

    vars%end_k = vars%num_local_rows - 1;
    print*, "myrank=", vars%myrank, "start_row=", vars%start_row, "end_row=", vars%end_row
    print*, "myrank=", vars%myrank, "num_local_rows=", vars%num_local_rows
    call rokko_distributed_mfree_construct_f(mat, simple, vars%dim, vars%num_local_rows)
  END SUBROUTINE initialize
  ! subroutine passed to C function.
  ! It must be interoperable!
  SUBROUTINE multiply (n, x, y) BIND(C)
    INTEGER(C_INT), INTENT(IN), VALUE :: n
    REAL(C_DOUBLE), INTENT(IN) :: x(n)
    REAL(C_DOUBLE), INTENT(OUT) :: y(n)
    integer :: ierr
    integer :: k
    print*, "calleddddddd"
    
    if (vars%num_local_rows == 0) then
       return
    endif
    
    if (.not.(vars%is_first_proc) .and. (vars%nprocs /= 1)) then
       !std::cout << "recv myrank=" << vars%myrank << std::endl;
       call MPI_Send(x(1), 1, MPI_DOUBLE, vars%myrank-1, 0, vars%comm);
       call MPI_Recv(vars%buf_m, 1, MPI_DOUBLE, vars%myrank-1, 0, vars%comm, vars%status_m, ierr);
       !std::cout << "buffff=" << buf << std::endl;
    endif
    
    if (.not.(vars%is_last_proc) .and. (vars%nprocs /= 1)) then
       !std::cout << "send myrank=" << vars%myrank << std::endl;
       call MPI_Recv(vars%buf_p, 1, MPI_DOUBLE, vars%myrank+1, 0, vars%comm, vars%status_p);
       call MPI_Send(x(vars%end_k), 1, MPI_DOUBLE, vars%myrank+1, 0, vars%comm);
       !std::cout << "buffff=" << buf2 << std::endl;
    endif
    
    if (vars%is_first_proc) then
       if (vars%num_local_rows /= 1) then
          y(2) = x(2) - x(1);
          if (vars%nprocs /= 1) then
             y(vars%end_k) = - x(vars%end_k - 1) + 2 * x(vars%end_k) - vars%buf_p
          endif
       else 
          y(2) = x(2) - vars%buf_p
       endif
    endif
    
    if (vars%is_last_proc) then
       if (vars%num_local_rows /= 1) then
          if (vars%nprocs /= 1) then
             y(1) = - vars%buf_m + 2 * x(1) - x(2)
          endif
          y(vars%end_k) = 2 * x(vars%end_k) - x(vars%end_k - 1)      
       else
          y(vars%end_k) = 2 * x(vars%end_k) - vars%buf_m
       endif
    endif
    if (.not.(vars%is_first_proc .or. vars%is_last_proc)) then ! neither first or last process
       if (vars%num_local_rows /= 1) then
          y(1) = - vars%buf_m + 2 * x(1) - x(2);
          y(vars%end_k) = - x(vars%end_k - 1) + 2 * x(vars%end_k) - vars%buf_p;
       endif
    else
       y(1) = - vars%buf_m + 2 * x(2) - vars%buf_p
    endif
    
    ! from 2 to end
    do k=2, vars%end_k-1
       y(k) = - x(k-1) + 2 * x(k) - x(k+1)
    enddo
  END SUBROUTINE multiply
  SUBROUTINE simple () BIND(C)
  END SUBROUTINE simple
END MODULE laplacian

program main
  use MPI
  use rokko
  use mfree
  use laplacian
  implicit none
  integer :: provided, ierr
  
  double precision :: eig_val
  double precision, allocatable, dimension(:) :: eig_vec

  type(rokko_parallel_sparse_solver) :: solver
  character(len=100) :: solver_name, tmp_str
  type(rokko_distributed_mfree) :: mat
  integer :: dim
  integer :: num_evals, block_size, max_iters
  integer :: num_local_rows, num_conv
  double precision :: tol

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)

  solver_name = "anasazi"
  dim = 10

  if (myrank == 0) then
     write(*,*) "solver name = ", trim(solver_name)
     write(*,*) "matrix dimension = ", dim
  endif

  call rokko_parallel_sparse_solver_construct(solver, solver_name)
  call initialize(mat, dim)

  num_evals = 10
  block_size = 5
  max_iters = 500
  tol = 1.0e-8
  call rokko_parallel_sparse_solver_diagonalize_distributed_mfree(solver, mat, num_evals, block_size, max_iters, tol)

  !num_conv = rokko_parallel_sparse_solver_num_conv(solver);
  !eig_val = rokko_parallel_sparse_solver_eigenvalue(solver, 0);
  !num_local_rows = rokko_distributed_crs_matrix_num_local_rows(mat);

  call rokko_distributed_mfree_destruct(mat)
  call rokko_parallel_sparse_solver_destruct(solver)

  call MPI_finalize(ierr)
end program main


