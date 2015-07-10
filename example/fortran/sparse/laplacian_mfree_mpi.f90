MODULE mfree
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE
  
  ! interface for C function.
  INTERFACE
     SUBROUTINE construct_c(func_multiply, vars) BIND(C)
       USE, INTRINSIC :: ISO_C_BINDING
       TYPE(C_FUNPTR), INTENT(IN), VALUE :: func_multiply
       TYPE(C_PTR), INTENT(INOUT) :: vars
     END SUBROUTINE construct_c
  END INTERFACE
  
CONTAINS 
  ! Wrapper for C function.
  SUBROUTINE construct(func_in, vars)
!    USE, INTRINSIC :: ISO_C_BINDING
    TYPE(C_FUNPTR) :: cproc
    TYPE(C_PTR), INTENT(inout) :: vars
    INTERFACE
       SUBROUTINE func_in(x, y, vars)
         USE, INTRINSIC :: ISO_C_BINDING
         REAL(KIND=C_DOUBLE), INTENT(IN), dimension(5) :: x
         REAL(KIND=C_DOUBLE), INTENT(OUT), dimension(5) :: y
         TYPE(C_PTR), INTENT(INOUT) :: vars
       END SUBROUTINE func_in
    END INTERFACE
    ! Get C procedure pointer.
    cproc = C_FUNLOC (func_in)

    ! call wrapper in C.
    call construct_c(cproc, vars)
  END SUBROUTINE construct
  
END MODULE mfree
 
MODULE laplacian
  use MPI
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
  TYPE(type_vars) vars
CONTAINS
  SUBROUTINE initialize (dim_in)
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
    print*, "myrank=", vars%myrank, "num_local_rows=%d\n", vars%num_local_rows
  END SUBROUTINE initialize
  ! subroutine passed to C function.
  ! It must be interoperable!
  SUBROUTINE multiply (x, y, vars_in) BIND(C)
    REAL(KIND=C_DOUBLE), INTENT(IN), dimension(5) :: x
    REAL(KIND=C_DOUBLE), INTENT(OUT), dimension(5) :: y
    TYPE(C_PTR), INTENT(INOUT) :: vars_in
    TYPE(type_vars), pointer :: vars
    integer :: ierr
    integer :: k

    call C_F_POINTER(vars_in, vars)
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
END MODULE laplacian

program main
  use MPI
  use mfree
  use laplacian
  implicit none
  integer :: provided, ierr
  TYPE(C_PTR) :: vars_main

  call MPI_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
  call initialize(10)
  call construct(multiply, vars_main)

  call MPI_finalize(ierr)
end program main


