#include <slepceps.h>
#include <petscblaslapack.h>
#include <vector>
#include <rokko/utility/lattice.hpp>
#include <boost/lexical_cast.hpp>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  Mat            A;               /* operator matrix */
  EPS            eps;             /* eigenproblem solver context */
  EPSType        type;
  PetscMPIInt    rank;
  PetscInt       nev;
  PetscErrorCode ierr;
  double init_tick, initend_tick, gen_tick, diag_tick, end_tick;

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  init_tick = MPI_Wtime();
  SlepcInitialize(&argc, &argv, (char*)0, 0);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
  initend_tick = MPI_Wtime();

  int len_ladder = 5;
  if (argc >= 2) len_ladder = boost::lexical_cast<int>(argv[1]);

  int L = 2 * len_ladder;
  std::vector<std::pair<int, int> > lattice;
  rokko::ladder_lattice_1dim(len_ladder, lattice);
  int dim = 1 << L;

  // Create Hermitean matrix
  gen_tick = MPI_Wtime();    
  ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
  ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, dim, dim); CHKERRQ(ierr);
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
      int m1 = 1 << i;
      int m2 = 1 << j;
      int m3 = m1 + m2;
      if (((k & m3) == m1) || ((k & m3) == m2)) {  // when (bit i == 1, bit j == 0) or (bit i == 0, bit j == 1)
        cols.push_back(k^m3);
        values.push_back(0.5);
        cols.push_back(k);
        values.push_back(-0.25);
      } else {
        cols.push_back(k);
        values.push_back(0.25);
      }
      ierr = MatSetValues(A, 1, &k, cols.size(), &cols[0], &values[0], ADD_VALUES); CHKERRQ(ierr);
    }
  }

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  //ierr = MatGetLocalSize(A, &n, NULL); CHKERRQ(ierr);
  //MatView(A, PETSC_VIEWER_STDOUT_WORLD);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create eigensolver context
  */
  diag_tick = MPI_Wtime();
  ierr = EPSCreate(PETSC_COMM_WORLD, &eps); CHKERRQ(ierr);
  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps, A, NULL); CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps, EPS_HEP); CHKERRQ(ierr);
  //ierr = EPSSetDimensions(eps, 10, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(ierr);
  //char routine_name[] = "lobpcg";
  //ierr = EPSSetType(eps,routine_name); CHKERRQ(ierr);
  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = EPSSolve(eps); CHKERRQ(ierr);
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
  ///ierr = EPSPrintSolution(eps,NULL); CHKERRQ(ierr);
  ierr = EPSDestroy(&eps); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = SlepcFinalize();
  MPI_Finalize();
  return 0;
}

