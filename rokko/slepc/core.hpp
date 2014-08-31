/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SLEPC_SOLVER_H
#define ROKKO_SLEPC_SOLVER_H

#include <rokko/slepc/distributed_crs_matrix.hpp>
#include <rokko/distributed_mfree.hpp>
#include "distributed_multivector.hpp"
#include <rokko/utility/timer.hpp>

#include <petscvec.h>
#include <slepceps.h>
#include <petscblaslapack.h>

namespace rokko {
namespace slepc {

#undef __FUNCT__
#define __FUNCT__ "MatMult_myMat"
PetscErrorCode MatMult_myMat(Mat A, Vec x, Vec y) {
  PetscFunctionBeginUser;
  PetscErrorCode ierr;
  
  rokko::distributed_mfree_slepc *op_ctx;
  ierr = MatShellGetContext(A, &op_ctx); CHKERRQ(ierr);
  
  PetscScalar const * px;
  PetscScalar * py;

  ierr = VecGetArrayRead(x, &px);  CHKERRQ(ierr);
  ierr = VecGetArray(y, &py);  CHKERRQ(ierr);
  
  op_ctx->multiply(px, py);
  
  ierr = VecRestoreArrayRead(x,&px);  CHKERRQ(ierr);
  ierr = VecRestoreArray(y,&py);  CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatGetDiagonal_myMat"
PetscErrorCode MatGetDiagonal_myMat(Mat A, Vec diag) {
  PetscFunctionBeginUser;
  PetscErrorCode ierr;

  rokko::distributed_mfree_slepc *op_ctx;
  ierr = MatShellGetContext(A, &op_ctx); CHKERRQ(ierr);

  PetscScalar *pd;
  
  ierr = VecGetArray(diag, &pd); CHKERRQ(ierr);
  op_ctx->diagonal(pd);
  ierr = VecRestoreArray(diag ,&pd); CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}


class solver {
public: 
  solver() {}
  ~solver() {}
  void initialize(int& argc, char**& argv) { SlepcInitialize(NULL, NULL, (char*)0, 0); }
  void finalize() {}
  void diagonalize(rokko::distributed_crs_matrix& mat,
                   distributed_multivector_slepc const& ivec,
                   int num_evals, int block_size, int max_iters, double tol) {
    Mat            A = *(reinterpret_cast<slepc::distributed_crs_matrix*>(mat.get_matrix())->get_matrix());
    EPSType        type;
    PetscMPIInt    size;
    PetscInt       nev;

    EPS            eps;             /* eigenproblem solver context */
    ierr = EPSCreate(PETSC_COMM_WORLD, &eps); //CHKERRQ(ierr);

    /* Set operators. In this case, it is a standard eigenvalue problem */
    ierr = EPSSetOperators(eps, A, NULL); //CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps, EPS_HEP); //CHKERRQ(ierr);
    //ierr = EPSSetDimensions(eps, num_evals, block_size, PETSC_DECIDE); //CHKERRQ(ierr);
    ierr = EPSSetDimensions(eps, num_evals, 2 * num_evals, PETSC_DECIDE); //CHKERRQ(ierr);
    ierr = EPSSetTolerances(eps, (PetscScalar) tol, (PetscInt) max_iters);
    /* Set solver parameters at runtime */
    ierr = EPSSetFromOptions(eps); //CHKERRQ(ierr);
    
    /* Solve the eigensystem */       
    ierr = EPSSolve(eps); //CHKERRQ(ierr);
    
    /* Optional: Get some information from the solver and display it */
    ierr = EPSGetType(eps, &type); //CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type); //CHKERRQ(ierr);
    ierr = EPSGetDimensions(eps, &nev, NULL, NULL); //CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev); //CHKERRQ(ierr);

    PetscInt nconv;
    EPSGetConverged(eps, &nconv);

    if (nconv == 0) {
      std::cout << "doesn't converge" << std::endl;
    }

    PetscScalar eval_r, eval_i;
    Vec evec_r, evec_i;
    MatGetVecs(A, NULL, &evec_r);
    MatGetVecs(A, NULL, &evec_i);
 
    for (PetscInt i = 0; i<nconv; ++i) {
      EPSGetEigenpair(eps, i, &eval_r, &eval_i, evec_r, evec_i);
      evals_.push_back((double)eval_r);
      //VecView(evec_r, PETSC_VIEWER_STDOUT_WORLD);
    }

    // Display solution and clean up
    //ierr = EPSPrintSolution(eps, NULL); //CHKERRQ(ierr);

    ierr = EPSDestroy(&eps); //CHKERRQ(ierr);
    ierr = VecDestroy(&evec_r); //CHKERRQ(ierr);
    ierr = VecDestroy(&evec_i); //CHKERRQ(ierr);
    ierr = MatDestroy(&A); //CHKERRQ(ierr);
  }

  void diagonalize(rokko::distributed_crs_matrix& mat,
                   distributed_multivector_anasazi const& ivec,
                   int num_evals, int block_size, int max_iters, double tol) {
  }

  void diagonalize(rokko::distributed_crs_matrix& mat,
                   distributed_multivector_anasazi const& ivec,
                   int num_evals, int block_size, int max_iters, double tol, timer& time) {
  }

  void diagonalize(rokko::distributed_mfree_slepc* const mat,
                   distributed_multivector_slepc const& ivec,
                   int num_evals, int block_size, int max_iters, double tol) {
    Mat            A;
    EPSType        type;
    PetscMPIInt    size;
    PetscInt       nev;

    // define matrix-free type operator
    std::cout << "mat->get_mapping_1d().num_local_rows()=" << mat->get_mapping_1d().num_local_rows() << std::endl;
    std::cout << " mat->get_mapping_1d().get_dim()=" <<  mat->get_mapping_1d().get_dim() << std::endl;
    ierr = MatCreateShell(PETSC_COMM_WORLD, mat->get_mapping_1d().num_local_rows(), mat->get_mapping_1d().num_local_rows(), mat->get_mapping_1d().get_dim(), mat->get_mapping_1d().get_dim(), mat, &A); //CHKERRQ(ierr);
    ierr = MatSetFromOptions(A); //CHKERRQ(ierr);

    ierr = MatShellSetOperation(A, MATOP_MULT, (void(*)())MatMult_myMat); //CHKERRQ(ierr);
    ierr = MatShellSetOperation(A, MATOP_MULT_TRANSPOSE, (void(*)())MatMult_myMat); //CHKERRQ(ierr);
    ierr = MatShellSetOperation(A, MATOP_GET_DIAGONAL, (void(*)())MatGetDiagonal_myMat); //CHKERRQ(ierr);

    EPS            eps;             /* eigenproblem solver context */
    ierr = EPSCreate(PETSC_COMM_WORLD, &eps); //CHKERRQ(ierr);
    /* Set operators. In this case, it is a standard eigenvalue problem */
    ierr = EPSSetOperators(eps, A, NULL); //CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps, EPS_HEP); //CHKERRQ(ierr);
    //ierr = EPSSetDimensions(eps, num_evals, block_size, PETSC_DECIDE); //CHKERRQ(ierr);
    ierr = EPSSetDimensions(eps, num_evals, 2 * num_evals, PETSC_DECIDE); //CHKERRQ(ierr);
    ierr = EPSSetTolerances(eps, (PetscScalar) tol, (PetscInt) max_iters);
    /* Set solver parameters at runtime */
    ierr = EPSSetFromOptions(eps); //CHKERRQ(ierr);
    
    /* Solve the eigensystem */       
    ierr = EPSSolve(eps); //CHKERRQ(ierr);
    
    /* Optional: Get some information from the solver and display it */
    ierr = EPSGetType(eps, &type); //CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type); //CHKERRQ(ierr);
    ierr = EPSGetDimensions(eps, &nev, NULL, NULL); //CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev); //CHKERRQ(ierr);

    PetscInt nconv;
    EPSGetConverged(eps, &nconv);

    if (nconv == 0) {
      std::cout << "doesn't converge" << std::endl;
    }

    PetscScalar eval_r, eval_i;
    Vec evec_r, evec_i;
    MatGetVecs(A, NULL, &evec_r);
    MatGetVecs(A, NULL, &evec_i);
 
    for (PetscInt i = 0; i<nconv; ++i) {
      EPSGetEigenpair(eps, i, &eval_r, &eval_i, evec_r, evec_i);
      evals_.push_back((double)eval_r);
      //VecView(evec_r, PETSC_VIEWER_STDOUT_WORLD);
    }

    // Display solution and clean up
    //ierr = EPSPrintSolution(eps, NULL); //CHKERRQ(ierr);

    ierr = EPSDestroy(&eps); //CHKERRQ(ierr);
    ierr = VecDestroy(&evec_r); //CHKERRQ(ierr);
    ierr = VecDestroy(&evec_i); //CHKERRQ(ierr);
    ierr = MatDestroy(&A); //CHKERRQ(ierr);
  }

  rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(mapping_1d const& map) {
    return new slepc::distributed_crs_matrix();
  }

  std::vector<double> eigenvalues() const { return evals_; }
  distributed_multivector_slepc eigenvectors() const { return evecs_; }
private:
  PetscErrorCode ierr;
  std::vector<double> evals_;
  distributed_multivector_slepc evecs_;
};

} // namespace slepc
} // namespace rokko

#endif // ROKKO_SLEPC_DISTRIBUTED_CRS_MATRIX_H
