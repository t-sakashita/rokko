#include <iostream>
#include <fstream>
#include <slepceps.h>
#include <petscblaslapack.h>
#include <rokko/rokko.hpp>
#include <rokko/utility/xyz_hamiltonian.hpp>

struct model {
  int L;
  std::vector<std::pair<int, int>> lattice;
  std::vector<std::tuple<double, double, double>> coupling;
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
  if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only");

  /* 
  char str[100];
  ierr = PetscOptionsGetString(NULL,"-filename", str, 100, NULL); CHKERRQ(ierr);
  //std::string filename(str);
  std::cout << "str=" << argv[1] << std::endl;

  std::ifstream ifs(argv[1]); //str);
  //  if (!ifs) {

  //}
  */

  int num_bonds;
  model m;
  std::cin >> m.L >> num_bonds;
  for (int i=0; i<num_bonds; ++i) {
    int j, k;
    std::cin >> j >> k;
    m.lattice.push_back(std::make_pair(j, k));
  }

  for (int i=0; i<num_bonds; ++i) {
    double jx, jy, jz;
    std::cin >> jx >> jy >> jz;
    m.coupling.push_back(std::make_tuple(jx, jy, jz));
  }

  std::cout << "L=" << m.L << " num_bonds=" << num_bonds << std::endl;
  for (int i=0; i<num_bonds; ++i) {
    std::cout << m.lattice[i].first << " " << m.lattice[i].second << " " << m.coupling[i].get<0>() << " " << m.coupling[i].get<1>() << " " << m.coupling[i].get<2>() << std::endl;
  }
  N = 1 << m.L;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

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

  ierr = EPSValuesView(eps,NULL);CHKERRQ(ierr);
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
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
  rokko::xyz_hamiltonian::multiply(m->L, m->lattice, m->coupling, px, py);
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
  rokko::xyz_hamiltonian::fill_diagonal(m->L, m->lattice, m->coupling, pd);
  ierr = VecRestoreArray(diag ,&pd); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

