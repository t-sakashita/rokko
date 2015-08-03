/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 by Synge Todo <wistaria@comp-phys.org>,
*                       Tatsuya Sakashita <t-sakashita@issp.u-tokyo.ac.jp>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SLEPC_CORE_HPP
#define ROKKO_SLEPC_CORE_HPP

#include <rokko/slepc/distributed_crs_matrix.hpp>
#include <rokko/distributed_mfree.hpp>

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
  solver() { SlepcInitialize(NULL, NULL, (char*)0, 0); }
  ~solver() {}
  void initialize(int& argc, char**& argv) { SlepcInitialize(NULL, NULL, (char*)0, 0); }
  void finalize() {}
  parameters diagonalize(rokko::distributed_crs_matrix& mat, int num_evals, int block_size,
			 int max_iters, double tol) {
    parameters params;
    params.set("num_eigenvalues", num_evals);
    params.set("Block Size", block_size);
    params.set("Maximum Iterations", max_iters);
    params.set("Convergence Tolerance", tol);
    return diagonalize(mat, params);
  }

  parameters diagonalize(rokko::distributed_mfree& mat, int num_evals, int block_size, int max_iters,
			 double tol) {
    parameters params;
    params.set("num_eigenvalues", num_evals);
    params.set("Block Size", block_size);
    params.set("Maximum Iterations", max_iters);
    params.set("Convergence Tolerance", tol);
    return diagonalize(mat, params);
  }

  parameters diagonalize(rokko::distributed_crs_matrix& mat, rokko::parameters const& params) {
    parameters params_out;
    dimension_ = mat.get_dim();
    offset_local_ = mat.start_row();
    num_local_rows_ = mat.num_local_rows();
    int block_size_;
    if (params.defined("Block Size"))  {
      block_size_ = params.get<int>("Block Size");
    }
    else {
      block_size_ = 1;
    }
    double tol;
    if (params.defined("Convergence Tolerance")) tol = params.get<double>("Convergence Tolerance");
    int max_iters;
    if (params.defined("Maximum Iterations")) max_iters = params.get<int>("Maximum Iterations");

    A = reinterpret_cast<slepc::distributed_crs_matrix*>(mat.get_matrix())->get_matrix();
    EPSType        type;
    PetscMPIInt    size;
    PetscInt       num_evals;
    if (params.defined("num_eigenvalues")) {
      num_evals = (PetscInt) params.get<int>("num_eigenvalues");
    } else {
      num_evals = 1;
    }
    ierr = EPSCreate(PETSC_COMM_WORLD, &eps);

    /* Set operators. In this case, it is a standard eigenvalue problem */
    ierr = EPSSetOperators(eps, *A, NULL);
    ierr = EPSSetProblemType(eps, EPS_HEP);
    if (params.defined("routine")) {
      if ((params.type("routine") != typeid(std::string)) && params.type("routine") != typeid(const char*))
	throw "error: routine must be charatcters or string.";
      routine_ = params.get_string("routine");
      ierr = EPSSetType(eps, (EPSType)routine_.c_str());
    }
    ierr = EPSSetDimensions(eps, num_evals, 2 * num_evals, PETSC_DECIDE);
    ierr = EPSSetTolerances(eps, (PetscScalar) tol, (PetscInt) max_iters);
    /* Set solver parameters at runtime */
    ierr = EPSSetFromOptions(eps);

    /* Solve the eigensystem */       
    ierr = EPSSolve(eps);
    
    /* Get some information from the solver and display it */
    ierr = EPSGetType(eps, &type);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);
    ierr = EPSGetDimensions(eps, &num_evals, NULL, NULL);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",num_evals);
    EPSGetConverged(eps, &num_conv_);
    params_out.set("num_conv", num_conv_);
    if (num_conv_ == 0) {
      std::cout << "doesn't converge" << std::endl;
    }
    return params_out;
  }

  parameters diagonalize(rokko::distributed_mfree& mat_in, rokko::parameters const& params) {
    parameters params_out;
    rokko::distributed_mfree* mat = &mat_in;
    // define matrix-free type operator
    dimension_ = mat->get_dim();
    offset_local_ = mat->get_local_offset();
    num_local_rows_ = mat->get_num_local_rows();
    A = new Mat();
    ierr = MatCreateShell(PETSC_COMM_WORLD, mat->get_num_local_rows(), mat->get_num_local_rows(), mat->get_dim(), mat->get_dim(), mat, A);
    ierr = MatSetFromOptions(*A);
    ierr = MatShellSetOperation(*A, MATOP_MULT, (void(*)())MatMult_myMat);
    ierr = MatShellSetOperation(*A, MATOP_MULT_TRANSPOSE, (void(*)())MatMult_myMat);
    ierr = MatShellSetOperation(*A, MATOP_GET_DIAGONAL, (void(*)())MatGetDiagonal_myMat);

    ierr = EPSCreate(PETSC_COMM_WORLD, &eps);
    /* Set operators. In this case, it is a standard eigenvalue problem */
    ierr = EPSSetOperators(eps, *A, NULL);

    int block_size_;
    if (params.defined("Block Size"))  {
      block_size_ = params.get<int>("Block Size");
    }
    else {
      block_size_ = 1;
    }
    double tol;
    if (params.defined("Convergence Tolerance")) tol = params.get<double>("Convergence Tolerance");
    int max_iters;
    if (params.defined("Maximum Iterations")) max_iters = params.get<int>("Maximum Iterations");
    
    EPSType        type;
    //PetscMPIInt    size;
    PetscInt       num_evals;
    if (params.defined("num_eigenvalues")) {
      num_evals = (PetscInt) params.get<int>("num_eigenvalues");
    } else {
      num_evals = 1;
    }

    /* Set operators. In this case, it is a standard eigenvalue problem */
    ierr = EPSSetProblemType(eps, EPS_HEP);
    if (params.defined("routine")) {
      if ((params.type("routine") != typeid(std::string)) && params.type("routine") != typeid(const char*))
	throw "error: routine must be charatcters or string.";
      routine_ = params.get_string("routine");
      ierr = EPSSetType(eps, (EPSType)routine_.c_str());
    }

    ierr = EPSSetDimensions(eps, num_evals, 2 * num_evals, PETSC_DECIDE);
    ierr = EPSSetTolerances(eps, (PetscScalar) tol, (PetscInt) max_iters);
    /* Set solver parameters at runtime */
    ierr = EPSSetFromOptions(eps);

    /* Solve the eigensystem */       
    ierr = EPSSolve(eps);
    
    /* Get some information from the solver and display it */
    ierr = EPSGetType(eps, &type);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);
    ierr = EPSGetDimensions(eps, &num_evals, NULL, NULL);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",num_evals);
    EPSGetConverged(eps, &num_conv_);
    params_out.set("num_conv", num_conv_);
    if (num_conv_ == 0) {
      std::cout << "doesn't converge" << std::endl;
    }
    return params_out;
  }

  double eigenvalue(int i) const {
    PetscScalar eval_r, eval_i;
    EPSGetEigenvalue(eps, i, &eval_r, &eval_i);
    return eval_r;
  }

  void eigenvector(int k, double* vec) const {
    Vec evec_r, evec_i;
    MatGetVecs(*A, NULL, &evec_r);
    MatGetVecs(*A, NULL, &evec_i);
    VecPlaceArray(evec_r, &vec[0]);
    EPSGetEigenvector(eps, k, evec_r, evec_i);
    VecDestroy(&evec_r);
    VecDestroy(&evec_i);
  }

  void eigenvector(int k, std::vector<double>& vec) const {
    if (vec.size() < num_local_rows_) vec.resize(num_local_rows_);
    eigenvector(k, &(vec[0]));
  }

  void eigenvector(int k, distributed_vector& vec) const {
    vec.initialize(dimension_, offset_local_, offset_local_ + num_local_rows_);
    eigenvector(k, vec.get_storage());
  }

  int num_conv() const {
    return num_conv_;
  }

  rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(int row_dim,
    int col_dim) {
    return new slepc::distributed_crs_matrix(row_dim, col_dim);
  }

private:
  int dimension_, offset_local_, num_local_rows_;
  Mat*           A;
  EPS            eps;             /* eigenproblem solver context */
  PetscErrorCode ierr;
  std::string routine_;
  EPSType routine_type;
  int num_conv_;
};

} // namespace slepc
} // namespace rokko

#endif // ROKKO_SLEPC_CORE_HPP
