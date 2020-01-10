/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SLEPC_CORE_HPP
#define ROKKO_SLEPC_CORE_HPP

#include <rokko/slepc.hpp>

#include <rokko/slepc/distributed_crs_matrix.hpp>
#include <rokko/slepc/distributed_mfree.hpp>

namespace rokko {

namespace slepc {

class solver {
public: 
  solver() {
    SlepcInitialize(NULL, NULL, (char*)0, 0);
    MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
  }
  ~solver() {
    if (!A) {
      delete A;
      A = nullptr;
    }
  }

  void initialize(int& argc, char**& argv) {}
  void finalize() { SlepcFinalize(); }

  static EPSType get_routine(rokko::parameters const& params) {
    if (params.defined("routine")) {
      if ((params.type("routine") != typeid(std::string)) && params.type("routine") != typeid(const char*))
        throw std::invalid_argument("slepc::get_routine() : routine must be charatcters or string.");
      std::string routine = params.get_string("routine");
      if (!routine.empty()) {
        return static_cast<EPSType>(routine.c_str());
      } else {
        return EPSKRYLOVSCHUR;
      }
    }
    else {
      return EPSKRYLOVSCHUR;
    }
  }

  static PetscInt get_num_eigvals(rokko::parameters const& params) {
    return params.defined("num_eigvals") ? (PetscInt) params.get<int>("num_eigvals") : 1;
  }

  static PetscInt get_max_block_size(rokko::parameters const& params) {
    return params.defined("max_block_size") ? params.get<int>("max_block_size") : PETSC_DECIDE;
  }

  static PetscReal get_conv_tol(rokko::parameters const& params) {
    return params.defined("conv_tol") ? params.get<double>("conv_tol") : (PetscReal)PETSC_DEFAULT;
  }

  static PetscInt get_max_iters(rokko::parameters const& params) {
    return params.defined("max_iters") ? params.get<int>("max_iters") : PETSC_DECIDE;
  }

  parameters diagonalize(rokko::distributed_crs_matrix& mat, rokko::parameters const& params) {
    PetscErrorCode ierr;
    parameters params_out;
    dimension_ = mat.get_dim();
    offset_local_ = mat.start_row();
    num_local_rows_ = mat.num_local_rows();

    PetscInt num_evals = get_num_eigvals(params);
    PetscInt max_block_size = get_max_block_size(params);
    PetscReal tol = get_conv_tol(params);
    PetscInt max_iters = get_max_iters(params);

    A = static_cast<slepc::distributed_crs_matrix&>(mat.get_matrix()).get_matrix();
    ierr = EPSCreate(PETSC_COMM_WORLD, &eps);

    /* Set operators for a standard eigenvalue problem */
    ierr = EPSSetOperators(eps, *A, NULL);
    ierr = EPSSetProblemType(eps, EPS_HEP);
    ierr = EPSSetType(eps, get_routine(params));
    ierr = EPSSetDimensions(eps, num_evals, max_block_size, PETSC_DECIDE);
    ierr = EPSSetTolerances(eps, tol, max_iters);
    /* Set solver parameters at runtime */
    ierr = EPSSetFromOptions(eps);

    /* Solve the eigensystem */       
    ierr = EPSSolve(eps);

    print_result(params_out);
    if (params.get_bool("verbose") && (myrank == 0)) {
      info_verbose();
    }
    return params_out;
  }

  parameters diagonalize(const rokko::distributed_mfree& mat_in, rokko::parameters const& params) {
    PetscErrorCode ierr;
    parameters params_out;
    rokko::distributed_mfree const*const mat = &mat_in;
    // define matrix-free type operator
    dimension_ = mat->get_dim();
    offset_local_ = mat->get_local_offset();
    num_local_rows_ = mat->get_num_local_rows();
    A = new Mat();
    ierr = MatCreateShell(PETSC_COMM_WORLD, mat->get_num_local_rows(), mat->get_num_local_rows(), mat->get_dim(), mat->get_dim(), const_cast<rokko::distributed_mfree*>(mat), A);
    ierr = MatSetFromOptions(*A);
    ierr = MatShellSetOperation(*A, MATOP_MULT, (void(*)())MatMult_myMat);
    ierr = MatShellSetOperation(*A, MATOP_MULT_TRANSPOSE, (void(*)())MatMult_myMat);
    //ierr = MatShellSetOperation(*A, MATOP_GET_DIAGONAL, (void(*)())MatGetDiagonal_myMat);

    ierr = EPSCreate(PETSC_COMM_WORLD, &eps);
    /* Set operators for a standard eigenvalue problem */
    ierr = EPSSetOperators(eps, *A, NULL);

    PetscInt num_evals = get_num_eigvals(params);
    PetscInt max_block_size = get_max_block_size(params);
    PetscReal tol = get_conv_tol(params);
    PetscInt max_iters = get_max_iters(params);

    ierr = EPSSetProblemType(eps, EPS_HEP);
    ierr = EPSSetType(eps, get_routine(params));
    ierr = EPSSetDimensions(eps, num_evals, max_block_size, PETSC_DECIDE);
    ierr = EPSSetTolerances(eps, tol, max_iters);
    /* Set solver parameters at runtime */
    ierr = EPSSetFromOptions(eps);

    /* Solve the eigensystem */       
    ierr = EPSSolve(eps);

    print_result(params_out);
    if (params.get_bool("verbose") && (myrank == 0)) {
      info_verbose();
    }
    return params_out;
  }

  void print_result(rokko::parameters& params_out) const {
    PetscErrorCode ierr;
    /* Get some information from the solver and display it */
    EPSType type;
    ierr = EPSGetType(eps, &type);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);
    PetscInt num_evals;
    ierr = EPSGetDimensions(eps, &num_evals, NULL, NULL);
    ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",num_evals);
    int num_conv = get_num_conv();
    params_out.set("num_conv", num_conv);
    if (num_conv == 0) {
      std::cerr << "doesn't converge" << std::endl;
    }
  }

  void info_verbose() const {
    PetscErrorCode ierr;
    PetscInt nev, ncv, mpd;
    PetscReal tol;
    PetscInt maxits, its;

    ierr = EPSGetDimensions(eps, &nev, &ncv, &mpd);
    ierr = EPSGetTolerances(eps, &tol, &maxits);
    ierr = EPSGetIterationNumber(eps, &its);
    std::cout << "number of eigenvalues to compute=" << nev << std::endl;
    std::cout << "maximum dimension of the subspace=" << ncv << std::endl;
    std::cout << "maximum dimension allowed for the projected problem=" << mpd << std::endl;
    std::cout << "convergence tolerance=" << tol << std::endl;
    std::cout << "maximum number of iterations=" << maxits << std::endl;
    std::cout << "number of iterations=" << its << std::endl;
  }
  
  double eigenvalue(int i) const {
    PetscScalar eval_r, eval_i;
    EPSGetEigenvalue(eps, i, &eval_r, &eval_i);
    return eval_r;
  }

  void eigenvector(int k, double *const vec) const {
    Vec evec_r, evec_i;
    MatCreateVecs(*A, NULL, &evec_r);
    MatCreateVecs(*A, NULL, &evec_i);
    VecPlaceArray(evec_r, vec);
    EPSGetEigenvector(eps, k, evec_r, evec_i);
    VecDestroy(&evec_r);
    VecDestroy(&evec_i);
  }

  void eigenvector(int k, std::vector<double>& vec) const {
    if (vec.size() < num_local_rows_) vec.resize(num_local_rows_);
    eigenvector(k, vec.data());
  }

  void eigenvector(int k, distributed_vector<double>& vec) const {
    vec.initialize(dimension_, offset_local_, offset_local_ + num_local_rows_);
    eigenvector(k, vec.get_storage());
  }

  int get_num_conv() const {
    PetscErrorCode ierr;
    PetscInt num_conv;
    ierr = EPSGetConverged(eps, &num_conv);
    return static_cast<int>(num_conv);
  }

  rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(int row_dim,
    int col_dim) {
    return new slepc::distributed_crs_matrix(row_dim, col_dim);
  }

  rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(int row_dim,
    int col_dim, int num_entries_per_row) { // fix me: accepting different number of entries for diagonal and off-diagonal blocks in PETSc
    return new slepc::distributed_crs_matrix(row_dim, col_dim, num_entries_per_row);
  }

private:
  int myrank;
  int dimension_, offset_local_, num_local_rows_;
  Mat*           A;
  EPS            eps;             /* eigenproblem solver context */
};

} // namespace slepc
} // namespace rokko

#endif // ROKKO_SLEPC_CORE_HPP
