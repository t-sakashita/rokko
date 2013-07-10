program frank_matrix
  use MPI
  use rokko
  implicit none
  integer::dim
  type(distributed_matrix)::Mat,Z !defined in rokko
  type(grid)::g !defined in rokko
  type(solver)::solver !defined in rokko
  
  real(8),allocatable::w(:),vec(:) !localized_vector
  char(len=100)::solver_name

  !--MPI variables--
  integer::ierr,irank,nprocs,comm

  !---loop variables---
  integer::i

  call MPI_init(ierr) 
  call MPI_comm_rank(MPI_COMM_WORLD, irank, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  comm = MPI_COMM_WORLD

  open(10, file="input.data")
  read(10,*) solver_name 
  read(10,*) dim
  close(10)

  call set_solver(solver, solver_name)
  call set_grid(g, comm)
  call initialize_distributed_matrix(mat, dim, dim, g)
  call initialize_distributed_matrix(Z, dim, dim, g)
  allocate(w(dim));
  
  call generate_franc_matrix(mat)

  call Diagnalize(solver, mat, w, Z)
  
  write(*,"(a,100e25.15)") "Eigen values=", w


  allocate(vec(dim))
  do i=1, dim
    call get_column_from_distributed_matrix(vec, Z, i)
    write(*,"(a,i,a,100e25.15)") "Eigen vector of ",i, "=", vec
  end do
    
  deallocate(w)
  deallocate(vec)
  call MPI_finalize(ierr)
end program frank_matrix
