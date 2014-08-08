#include <slepceps.h>
#include <petscblaslapack.h>
#include <rokko/localized_matrix.hpp>
#include <rokko/utility/heisenberg_hamiltonian_mpi.hpp>

struct model {
  MPI_Comm comm;
  int L;
  std::vector<std::pair<int, int> > lattice;
  double* buffer;
};

/*
   User-defined routines
*/
PetscErrorCode MatMult_myMat(Mat A, Vec x, Vec y);
PetscErrorCode MatGetDiagonal_myMat(Mat A, Vec diag);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  PetscMPIInt    nproc;
  PetscInt       nev;
  PetscErrorCode ierr;

  SlepcInitialize(&argc, &argv, (char*)0, 0);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&nproc); CHKERRQ(ierr);

  int L = 8;
  ierr = PetscOptionsGetInt(NULL,"-L", &L, NULL); CHKERRQ(ierr);
  PetscInt N_global = 1 << L;
  PetscInt N_local = N_global / nproc;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  model m;
  m.comm = PETSC_COMM_WORLD;
  m.L = L;
  m.buffer = new double[N_local];
  if (m.buffer == 0) {
    MPI_Abort(MPI_COMM_WORLD, 4);
  }

  for (int i=0; i<L; ++i) {
    m.lattice.push_back(std::make_pair(i, (i+1)%L));
  }

  ierr = MatCreateShell(PETSC_COMM_WORLD, N_local, N_local, N_global, N_global, &m, &A); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)())MatMult_myMat); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(void(*)())MatMult_myMat); CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_GET_DIAGONAL,(void(*)())MatGetDiagonal_myMat); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps); CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,NULL); CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP); CHKERRQ(ierr);
  ierr = EPSSetDimensions(eps, 5, 100, 100); CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps, (PetscScalar) 1., (PetscInt) 200);   CHKERRQ(ierr);
  Vec v0;
  MatGetVecs(A, &v0, NULL);
  VecSet(v0,1.0);
  EPSSetInitialSpace(eps,1,&v0);
  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSSolve(eps);CHKERRQ(ierr);

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetType(eps,&type); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type); CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSPrintSolution(eps,NULL); CHKERRQ(ierr);
  ierr = EPSDestroy(&eps); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "MatMult_myMat"
PetscErrorCode MatMult_myMat(Mat A,Vec x,Vec y)
{
  PetscFunctionBeginUser;
  PetscErrorCode ierr;
  model *m;
  ierr = MatShellGetContext(A, &m); CHKERRQ(ierr);

  const PetscScalar *px;
  PetscScalar       *py;

  ierr = VecGetArrayRead(x, &px); CHKERRQ(ierr);
  ierr = VecGetArray(y, &py); CHKERRQ(ierr);
  rokko::heisenberg_hamiltonian::multiply(m->comm, m->L, m->lattice, px, py, m->buffer);
  ierr = VecRestoreArrayRead(x,&px); CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&py); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MatGetDiagonal_myMat"
PetscErrorCode MatGetDiagonal_myMat(Mat A, Vec diag)
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

