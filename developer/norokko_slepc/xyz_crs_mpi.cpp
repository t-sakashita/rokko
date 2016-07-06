#include <slepceps.h>
#include <petscblaslapack.h>
#include <rokko/localized_matrix.hpp>
#include <rokko/utility/heisenberg_hamiltonian.hpp>
#include <vector>
#include <boost/tuple/tuple.hpp>

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
  ierr = PetscOptionsGetInt(NULL,NULL,"-L", &L, NULL); CHKERRQ(ierr);
  N = 1 << L;
  std::vector<std::pair<int, int> > lattice;
  std::vector<boost::tuple<double, double, double> > coupling;
  for (int i=0; i<L; ++i) {
    lattice.push_back(std::make_pair(i, (i+1)%L));
    //coupling.push_back(boost::make_tuple(1, 1, 1));
    coupling.push_back(boost::make_tuple(0.8, 0.3, 0.1));
  }

  // Create Hermitean matrix                                                           
  ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
  ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);
  ierr = MatSetUp(A); CHKERRQ(ierr);

  PetscInt Istart, Iend;
  ierr = MatGetOwnershipRange(A, &Istart, &Iend); CHKERRQ(ierr);

  std::vector<PetscInt> cols;
  std::vector<double> values;

  for (int l=0; l<lattice.size(); ++l) {
    for (int k=Istart; k<Iend; ++k) {
      cols.clear();
      values.clear();
      int i = lattice[l].first;
      int j = lattice[l].second;
      double jx = coupling[l].get<0>();
      double jy = coupling[l].get<1>();
      double jz = coupling[l].get<2>();
      double diag_plus = jz / 4.0;
      double diag_minus = - jz/ 4.0;
      double offdiag_plus = (jx + jy) / 4.0;
      double offdiag_minus = (jx - jy) / 4.0;

      int m1 = 1 << i;
      int m2 = 1 << j;
      int m3 = m1 + m2;
      if (((k & m3) == m1) || ((k & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        cols.push_back(k^m3);
        values.push_back(offdiag_plus);
        cols.push_back(k);
        values.push_back(diag_minus);
      } else {
        cols.push_back(k^m3);
        values.push_back(offdiag_minus);
        cols.push_back(k);
        values.push_back(diag_plus);
      }
      ierr = MatSetValues(A, 1, &k, cols.size(), &cols[0], &values[0], ADD_VALUES); CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  //ierr = MatGetLocalSize(A, &n, NULL); CHKERRQ(ierr);
  MatView(A, PETSC_VIEWER_STDOUT_WORLD);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD, &eps); CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps, A, NULL); CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, EPS_HEP); CHKERRQ(ierr);
  ierr = EPSSetDimensions(eps, 10, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSSolve(eps); CHKERRQ(ierr);

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
  ierr = EPSValuesView(eps,NULL); CHKERRQ(ierr);
  ierr = EPSDestroy(&eps); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return 0;
}

