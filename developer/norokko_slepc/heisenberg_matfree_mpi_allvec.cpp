#include <slepceps.h>
#include <petscblaslapack.h>
#include <rokko/localized_matrix.hpp>
#include <rokko/utility/heisenberg_hamiltonian.hpp>

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
  PetscMPIInt    nproc, myproc;
  PetscInt       nev;
  PetscErrorCode ierr;

  SlepcInitialize(&argc, &argv, (char*)0, 0);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&nproc); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myproc); CHKERRQ(ierr);

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
  ierr = EPSSetTolerances(eps, (PetscScalar) 1e-1, (PetscInt) 100);   CHKERRQ(ierr);
  //Vec v0;
  //MatCreateVecs(A, &v0, NULL);
  //VecSet(v0,1.0);
  //EPSSetInitialSpace(eps,1,&v0);
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

  PetscInt len_x, len_y;
  ierr = VecGetLocalSize(x, &len_x);  CHKERRQ(ierr);
  ierr = VecGetLocalSize(y, &len_y); CHKERRQ(ierr);
  //std::cout << "len_x=" << len_x << "len_y=" << len_y << std::endl;

  const PetscScalar *px;
  PetscScalar       *py, *pgx;

  VecScatter ctx;
  Vec vout;
  VecScatterCreateToAll(x,&ctx,&vout);
  // scatter as many times as you need
  VecScatterBegin(ctx,x,vout,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,x,vout,INSERT_VALUES,SCATTER_FORWARD);

  ierr = VecGetArray(vout, &pgx); CHKERRQ(ierr);
  ierr = VecGetArrayRead(x, &px); CHKERRQ(ierr);
  ierr = VecGetArray(y, &py); CHKERRQ(ierr);

  PetscInt Istart, Iend;
  ierr = MatGetOwnershipRange(A, &Istart, &Iend); CHKERRQ(ierr);
  PetscInt A_m, n;
  ierr = MatGetSize(A, &A_m, &n); CHKERRQ(ierr);
  PetscScalar* pgy = new PetscScalar[n];
  rokko::heisenberg_hamiltonian::multiply(m->L, m->lattice, pgx, pgy);

  for(int j = 0; j < len_y; ++j) {
    py[j] = 0.;
  }

  int k=0;
  for(int i = Istart; i < Iend; ++i) {
    py[k] = pgy[i];
    ++k;
  }

  // destroy scatter context and local vector when no longer needed
  VecScatterDestroy(&ctx);
  //VecView(vout, PETSC_VIEWER_STDOUT_WORLD);
  VecDestroy(&vout);

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

  PetscInt Istart, Iend;
  ierr = MatGetOwnershipRange(A, &Istart, &Iend); CHKERRQ(ierr);
  PetscInt A_m, n;
  ierr = MatGetSize(A, &A_m, &n); CHKERRQ(ierr);

  PetscScalar* pgd = new PetscScalar[n];
  rokko::heisenberg_hamiltonian::fill_diagonal(m->L, m->lattice, pgd);

  int k=0;
  for(int i = Istart; i < Iend; ++i) {
    pd[k] = pgd[i];
    ++k;
  }

  ierr = VecRestoreArray(diag ,&pd); CHKERRQ(ierr);
 
  PetscFunctionReturn(0);
}

