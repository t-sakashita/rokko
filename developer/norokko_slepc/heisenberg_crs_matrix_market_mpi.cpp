#include <slepceps.h>
#include <petscblaslapack.h>
#include <rokko/localized_matrix.hpp>
#include <rokko/utility/heisenberg_hamiltonian.hpp>
#include <vector>

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
  for (int i=0; i<L; ++i) {
    lattice.push_back(std::make_pair(i, (i+1)%L));
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
  MatView(A, PETSC_VIEWER_STDOUT_WORLD);

  PetscInt num_cols, rstart, rend;
  const PetscInt * colus;
  const PetscScalar * vals;

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  std::cout << "%%MatrixMarket matrix coordinate real general" << std::endl;
  MatGetOwnershipRange(A,&rstart,&rend);
  MatInfo info;
  MatGetInfo(A,MAT_LOCAL,&info);
  int num_nnz = info.nz_used;
  std::cout << "num_nnz=" << num_nnz << std::endl;
  for (int global_row=0; global_row<N; ++global_row) {
    if ((global_row >= rstart) && (global_row < rend)) {
      MatGetRow(A, global_row, &num_cols, &colus, &vals);
      for (int i=0; i<num_cols; ++i) {
	std::cout << global_row + 1 << " " << colus[i] + 1 << " " << vals[i] << std::endl;
      }
      MatRestoreRow(A, global_row, &num_cols, &colus, &vals);
    }
    MPI_Barrier(PETSC_COMM_WORLD);
  }
  
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return 0;
}

