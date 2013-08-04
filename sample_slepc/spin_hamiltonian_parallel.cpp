#include <slepceps.h>
#include <petscblaslapack.h>
#include <rokko/localized_matrix.hpp>
#include <rokko/utility/spin_hamiltonian.hpp>

struct model {
  int L;
  std::vector<std::pair<int, int> > lattice;
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
  PetscMPIInt    size;
  PetscInt       N, nev;
  PetscErrorCode ierr;

  SlepcInitialize(&argc, &argv, (char*)0, 0);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

  int L = 8;
  ierr = PetscOptionsGetInt(NULL,"-L", &L, NULL); CHKERRQ(ierr);
  N = 1 << L;

  // Create Hermitean matrix                                                           
  ierr = MatCreate(PETSC_COMM_WORLD, ctx->T);CHKERRQ(ierr);
  ierr = MatSetSizes(ctx->T,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(ctx->T);CHKERRQ(ierr);
  ierr = MatSetUp(ctx->T);CHKERRQ(ierr);
  
  ierr = MatGetOwnershipRange(ctx->T,&Istart,&Iend);CHKERRQ(ierr);
  value[0]=1.0; value[1]=-2.0; value[2]=1.0;
  for (i=(FirstBlock? Istart+1: Istart); i<(LastBlock? Iend-1: Iend); i++) {
    col[0]=i-1; col[1]=i; col[2]=i+1;


  model m;
  m.L = L;
  for (int i=0; i<L-1; ++i) {
    m.lattice.push_back(std::make_pair(i, i+1));
  }

  int N = 1 << L;
  for (int k1=0; k1<N; ++k1) {
    for (int k2=0; k2<N; ++k2) {
      for (int l=0; l<lattice.size(); ++l) {
        int i = lattice[l].first;
        int j = lattice[l].second;
        //cout << "k=" << k << " i=" << i << " j=" << j << endl;                                                                                                                
        int m1 = 1 << i;
        int m2 = 1 << j;
        int m3 = m1 + m2;
        if ((k2 & m3) == m1) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
          ierr = MatSetValues(ctx->T, 1, &i, 3, col, value, INSERT_VALUES); CHKERRQ(ierr); 
          //          mat(k1,k2) -= 0.25;

        }
        else if ((k2 & m3) == m2) {
          mat(k1,k2) -= 0.25;
        }
        else {
          mat(k1,k2) += 0.25;
          //      cout << "else" << endl;                                                                                                                                       
        }
      }
    }
  }

  ierr = MatCreateShell(PETSC_COMM_WORLD, N, N, N, N, &m, &A); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)())MatMult_myMat);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT_TRANSPOSE,(void(*)())MatMult_myMat);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_GET_DIAGONAL,(void(*)())MatGetDiagonal_myMat);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);

  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);

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
  ierr = MatMult(A, x, y); CHKERRQ(ierr);
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
  rokko::spin_hamiltonian::fill_diagonal(m->L, m->lattice, pd);
  ierr = VecRestoreArray(diag ,&pd); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

