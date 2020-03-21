#include <slepceps.h>
#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>
#include <rokko/utility/lattice.hpp>

struct model {
  MPI_Comm comm;
  int L;
  std::vector<std::pair<int, int>> lattice;
  double* buffer;
};

/*
   User-defined routines
*/
PetscErrorCode MatMult(Mat A, Vec x, Vec y);
PetscErrorCode MatGetDiagonal(Mat A, Vec diag);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  PetscMPIInt    rank, nproc;
  PetscInt       nev;
  PetscErrorCode ierr;
  double init_tick, initend_tick, gen_tick, diag_tick, end_tick;

  int provided;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  init_tick = MPI_Wtime();
  SlepcInitialize(&argc, &argv, (char*)0, 0);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&nproc); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
  MPI_Barrier(MPI_COMM_WORLD);
  initend_tick = MPI_Wtime();

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  model m;
  std::string lattice_file("xyz.dat");
  if (argc >= 2) lattice_file = argv[1];
  int L;
  rokko::read_lattice_file(lattice_file, L, m.lattice);
  
  MPI_Barrier(MPI_COMM_WORLD);
  gen_tick = MPI_Wtime();
  int dim = 1 << L;
  PetscInt N_global = dim;
  PetscInt N_local = N_global / nproc;
  m.L = L;
  m.comm = PETSC_COMM_WORLD;
  m.buffer = new double[N_local];
  if (m.buffer == 0) {
    MPI_Abort(MPI_COMM_WORLD, 4);
  }
  ierr = MatCreateShell(PETSC_COMM_WORLD, N_local, N_local, N_global, N_global, &m, &A); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)())MatMult); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(void(*)())MatMult); CHKERRQ(ierr);
  //ierr = MatShellSetOperation(A,MATOP_GET_DIAGONAL,(void(*)())MatGetDiagonal); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  MPI_Barrier(MPI_COMM_WORLD);
  diag_tick = MPI_Wtime();
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps); CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,NULL); CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP); CHKERRQ(ierr);
  //ierr = EPSSetDimensions(eps, 5, 100, 100); CHKERRQ(ierr);
  //ierr = EPSSetDimensions(eps, 1, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);
  //ierr = EPSSetTolerances(eps, (PetscScalar) 1., (PetscInt) 2000);   CHKERRQ(ierr);
  /*  Vec v0;
  MatCreateVecs(A, &v0, NULL);
  VecSet(v0,1.0);
  EPSSetInitialSpace(eps,1,&v0);*/

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  MPI_Barrier(MPI_COMM_WORLD);
  end_tick = MPI_Wtime();
  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetType(eps,&type); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"solution method: %s\n\n",type); CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"number of requested eigenvalues: %D\n",nev); CHKERRQ(ierr);
  int num_conv;
  EPSGetConverged(eps, &num_conv);
  PetscScalar eval_r, eval_i;
  if (rank == 0) {
    std::cout << "number of converged eigenpairs = " << num_conv << std::endl;
    for (int i=0; i<num_conv; ++i) {
      EPSGetEigenvalue(eps, i, &eval_r, &eval_i);
      std::cout << ' ' << eval_r;
    }
    std::cout << std::endl;
    std::cout << "init_time = " << initend_tick - init_tick << std::endl
	      << "gen_time = " << diag_tick - gen_tick << std::endl
	      << "diag_time = " << end_tick - diag_tick << std::endl;
  }
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  //ierr = EPSValuesView(eps,NULL); CHKERRQ(ierr);
  ierr = EPSDestroy(&eps); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = SlepcFinalize();

  MPI_Finalize();
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "MatMult"
PetscErrorCode MatMult(Mat A,Vec x,Vec y)
{
  PetscFunctionBeginUser;
  PetscErrorCode ierr;
  model *m;
  ierr = MatShellGetContext(A, &m); CHKERRQ(ierr);

  const PetscScalar *px;
  PetscScalar       *py;

  ierr = VecGetArrayRead(x, &px); CHKERRQ(ierr);
  ierr = VecGetArray(y, &py); CHKERRQ(ierr);
  PetscInt len_x, len_y;
  ierr = VecGetLocalSize(x, &len_x);  CHKERRQ(ierr);
  ierr = VecGetLocalSize(y, &len_y); CHKERRQ(ierr);
  //for(int j = 0; j < len_y; ++j) {
  //  py[j] = 0.;
  //}
  rokko::heisenberg_hamiltonian::multiply(m->comm, m->L, m->lattice, px, py, m->buffer);
  ierr = VecRestoreArrayRead(x,&px); CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&py); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetDiagonal"
PetscErrorCode MatGetDiagonal(Mat A, Vec diag)
{
  PetscFunctionBeginUser;
  PetscErrorCode ierr;
  model *m;
  ierr = MatShellGetContext(A, &m); CHKERRQ(ierr);

  PetscScalar       *pd;

  ierr = VecGetArray(diag, &pd); CHKERRQ(ierr);
  rokko::heisenberg_hamiltonian::fill_diagonal(m->comm, m->L, m->lattice, pd);
  ierr = VecRestoreArray(diag ,&pd); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
