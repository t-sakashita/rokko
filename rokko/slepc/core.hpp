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

#ifndef ROKKO_ANASAZI_SOLVER_H
#define ROKKO_ANASAZI_SOLVER_H

#include <rokko/distributed_mfree.hpp>
#include "distributed_crs_matrix.hpp"
#include "distributed_multivector.hpp"

#include <slepceps.h>
#include <petscblaslapack.h>

namespace rokko {

class solver_slepc {
public:
  solver_slepc() {
    SlepcInitialize(NULL, NULL, (char*)0, 0);
  }

  ~solver_slepc() {
    //ierr = SlepcFinalize();
  }

  void diagonalize(distributed_crs_matrix const& mat,
                   distributed_multivector_slepc const& ivec,
                   int num_evals, int block_size, int max_iters, double tol) {
    Mat            A = reinterpret_cast<slepc::distributed_crs_matrix*>(mat.matrix_impl_.get())->matrix_;          
    EPSType        type;
    PetscMPIInt    size;
    PetscInt       nev;

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

  void diagonalize(rokko::distributed_mfree* const mat,
                   distributed_multivector_slepc const& ivec,
                   int num_evals, int block_size, int max_iters, double tol) {
  }

  std::vector<double> eigenvalues() const { return evals_; }
  distributed_multivector_slepc eigenvectors() const { return evecs_; }
private:
  PetscErrorCode ierr;
  EPS            eps;             /* eigenproblem solver context */
  std::vector<double> evals_;
  distributed_multivector_slepc evecs_;
};

} // namespace rokko

#endif // ROKKO_SLEPC_DISTRIBUTED_CRS_MATRIX_H
