/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2020 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_SLEPC_SOLVER_HPP
#define ROKKO_SLEPC_SOLVER_HPP

#include <rokko/slepc.hpp>

#include <rokko/slepc/distributed_crs_matrix.hpp>
#include <rokko/slepc/distributed_mfree.hpp>

namespace rokko {

namespace slepc {

class solver {
public:
  static const std::map<std::string,EPSLanczosReorthogType> lanczos_reorthog_map;

  solver() {
    SlepcInitialize(NULL, NULL, (char*)NULL, NULL);
  }
  ~solver() = default;

  void initialize(int& argc, char**& argv) {}
  void finalize() { SlepcFinalize(); }

  static std::string get_routine(rokko::parameters const& params) {
    if (params.defined("routine")) {
      if ((params.type("routine") != typeid(std::string)) && params.type("routine") != typeid(const char*))
        throw std::invalid_argument("slepc::get_routine() : routine must be charatcters or string.");
      std::string routine = params.get_string("routine");
      if (!routine.empty()) {
        return routine;
      } else {
        return EPSKRYLOVSCHUR;
      }
    } else {
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

  static EPSWhich get_wanted_eigenvalues(std::string const& str) {
    if ((str == "largest") || (str == "EPS_LARGEST_MAGNITUDE"))
      return EPS_LARGEST_MAGNITUDE;
    if ((str == "smallest") || (str == "EPS_SMALLEST_MAGNITUDE"))
      return EPS_SMALLEST_MAGNITUDE;
    if ((str == "largest_real") || (str == "EPS_LARGEST_REAL"))
      return EPS_LARGEST_REAL;
    if ((str == "smallest_real") || (str == "EPS_SMALLEST_REAL"))
      return EPS_SMALLEST_REAL;
    else
      throw std::invalid_argument("get_wanted_eigenvalues: invalid parameter");
  }

  void set_lanczos_reorthog(rokko::parameters const& params) {
    if (params.defined("reorthog")) {
      std::string str = params.get<std::string>("reorthog");
      if (lanczos_reorthog_map.find(str) == lanczos_reorthog_map.end())
        throw std::invalid_argument(str + " is not lanczos reorthog in SLEPc");
      PetscErrorCode ierr = EPSLanczosSetReorthog(eps, lanczos_reorthog_map.at(str));
    }
  }

  void set_block_size(rokko::parameters const& params) {
    PetscErrorCode ierr;
    if (params.defined("block_size")) {
      int block_size = params.get<int>("block_size");
      if (get_routine(params) == "lobpcg") {
        //ierr = EPSLOBPCGSetBlockSize(eps, block_size);
      } else if (get_routine(params) == "jd") {
        ierr = EPSJDSetBlockSize(eps, block_size);
      }
    }
  }

  void set_restart(rokko::parameters const& params) {
    PetscErrorCode ierr;
    if (params.defined("restart")) {
      double restart = params.get<double>("restart");
      if (get_routine(params) == "lobpcg") {
        ierr = EPSLOBPCGSetRestart(eps, restart);
      } else if (get_routine(params) == "krylovschur") {
        ierr = EPSKrylovSchurSetRestart(eps, restart);
      } else if (get_routine(params) == "jd") {
        //ierr = EPSJDSetRestart(eps, restart);
      }
    }
  }

  void set_locking(rokko::parameters const& params) {
    PetscErrorCode ierr;
    if (params.defined("locking")) {
      PetscBool locking = params.get<bool>("locking") ? PETSC_TRUE : PETSC_FALSE;
      if (get_routine(params) == "lobpcg") {
        ierr = EPSLOBPCGSetLocking(eps, locking);
      } else if (get_routine(params) == "krylovschur") {
        ierr = EPSKrylovSchurSetLocking(eps, locking);
      }
    }
  }

  void set_wanted_eigenvalues(rokko::parameters const& params) {
    std::string str = params.defined("wanted_eigenvalues") ? params.get_string("wanted_eigenvalues") : std::string{};
    if (!str.empty()) {
      EPSSetWhichEigenpairs(eps, get_wanted_eigenvalues(str));
    }
  }

  parameters diagonalize(const rokko::slepc::distributed_crs_matrix& mat, rokko::parameters const& params) {
    dimension_ = mat.get_dim();
    offset_local_ = mat.start_row();
    num_local_rows_ = mat.get_num_local_rows();
    comm_ = mat.get_map().get_mpi_comm().get_comm();

    return diagonalize_common(mat.get_matrix(), params);
  }

  parameters diagonalize(const rokko::distributed_crs_matrix& mat, rokko::parameters const& params) {
    if (mat.get_solver_name() != "slepc") {
      throw std::invalid_argument("rokko::slepc::solver::diagonalize() : " + mat.get_solver_name() + "'s distributed_crs_matrix is given.");
    }
    return diagonalize(*std::static_pointer_cast<const rokko::slepc::distributed_crs_matrix>(mat.get_ptr()), params);
  }

  parameters diagonalize(const rokko::distributed_mfree& mat, rokko::parameters const& params) {
    Mat A;
    PetscErrorCode ierr;
    ierr = MatCreateShell(mat.get_comm(), mat.get_num_local_rows(), mat.get_num_local_rows(), mat.get_dim(), mat.get_dim(), const_cast<rokko::distributed_mfree*>(&mat), &A);
    ierr = MatSetFromOptions(A);
    ierr = MatShellSetOperation(A, MATOP_MULT, (void(*)())MatMult);
    ierr = MatShellSetOperation(A, MATOP_MULT_TRANSPOSE, (void(*)())MatMult);
    ierr = MatShellSetOperation(A, MATOP_GET_DIAGONAL, (void(*)())MatGetDiagonal);

    dimension_ = mat.get_dim();
    offset_local_ = mat.get_local_offset();
    num_local_rows_ = mat.get_num_local_rows();
    comm_ = mat.get_comm();

    return diagonalize_common(A, params);
  }

  #undef __FUNCT__
  #define __FUNCT__ "slepc::solver::MatMult"
  static PetscErrorCode MatMult(Mat A, Vec x, Vec y) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    rokko::distributed_mfree *op_ctx;
    ierr = MatShellGetContext(A, &op_ctx); CHKERRQ(ierr);

    PetscScalar const* px;
    PetscScalar * py;

    ierr = VecGetArrayRead(x, &px);  CHKERRQ(ierr);
    ierr = VecGetArray(y, &py);  CHKERRQ(ierr);

    op_ctx->multiply(px, py);

    ierr = VecRestoreArrayRead(x,&px);  CHKERRQ(ierr);
    ierr = VecRestoreArray(y,&py);  CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  #undef __FUNCT__
  #define __FUNCT__ "slepc::solver::MatGetDiagonal"
  static PetscErrorCode MatGetDiagonal(Mat A, Vec diag) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;

    rokko::distributed_mfree *op_ctx;
    ierr = MatShellGetContext(A, &op_ctx); CHKERRQ(ierr);

    PetscScalar *pd;

    ierr = VecGetArray(diag, &pd); CHKERRQ(ierr);
    op_ctx->fill_diagonal(pd);
    ierr = VecRestoreArray(diag ,&pd); CHKERRQ(ierr);

    PetscFunctionReturn(0);
  }

  parameters diagonalize_common(Mat const& A, rokko::parameters const& params) {
    parameters params_out;

    PetscInt num_evals = get_num_eigvals(params);
    PetscInt max_block_size = get_max_block_size(params);
    PetscReal tol = get_conv_tol(params);
    PetscInt max_iters = get_max_iters(params);

    PetscErrorCode ierr;
    ierr = EPSCreate(comm_, &eps);

    // Set operators for a standard eigenvalue problem
    ierr = EPSSetOperators(eps, A, NULL);
    ierr = EPSSetProblemType(eps, EPS_HEP);
    ierr = EPSSetType(eps, get_routine(params).c_str());
    ierr = EPSSetDimensions(eps, num_evals, max_block_size, PETSC_DECIDE);
    ierr = EPSSetTolerances(eps, tol, max_iters);
    set_wanted_eigenvalues(params);

    set_block_size(params);
    set_restart(params);
    set_locking(params);

    if (get_routine(params) == "lanczos")
      set_lanczos_reorthog(params);

    // Set solver parameters at runtime
    ierr = EPSSetFromOptions(eps);

    // Solve the eigensystem
    ierr = EPSSolve(eps);

    set_output_parameters(params_out);

    int myrank;
    MPI_Comm_rank(comm_, &myrank);
    if (params.get_bool("verbose") && (myrank == 0))
      info_verbose();

    return params_out;
  }

  void set_output_parameters(rokko::parameters& params_out) const {
    int num_conv = get_num_conv();
    params_out.set("num_conv", num_conv);
  }

  void info_verbose() const {
    PetscErrorCode ierr;

    EPSType type;
    ierr = EPSGetType(eps, &type);

    PetscInt nev, ncv, mpd;
    PetscReal tol;
    PetscInt maxits, its;

    ierr = EPSGetDimensions(eps, &nev, &ncv, &mpd);
    ierr = EPSGetTolerances(eps, &tol, &maxits);
    ierr = EPSGetIterationNumber(eps, &its);

    std::cout << "input parameters:" << std::endl;
    std::cout << "  solution method: " << type << std::endl;
    std::cout << "  number of eigenvalues to compute=" << nev << std::endl;
    std::cout << "  maximum dimension of the subspace=" << ncv << std::endl;
    std::cout << "  maximum dimension allowed for the projected problem=" << mpd << std::endl;
    std::cout << "  convergence tolerance=" << tol << std::endl;
    std::cout << "  maximum number of iterations=" << maxits << std::endl;
    std::cout << "results:" << std::endl;
    std::cout << "  number of converged eigenpairs=" << get_num_conv() << std::endl;
    std::cout << "  number of iterations=" << its << std::endl;
  }
  
  double eigenvalue(int i) const {
    PetscScalar eval_r, eval_i;
    EPSGetEigenvalue(eps, i, &eval_r, &eval_i);
    return eval_r;
  }

  void eigenvector(int k, double *const vec) const {
    Vec evec_r;
    VecCreateMPIWithArray(comm_, PETSC_DECIDE, num_local_rows_, dimension_, vec, &evec_r);
    EPSGetEigenvector(eps, k, evec_r, NULL);
    VecDestroy(&evec_r);
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

  const EPS& get_object() const {
    return eps;
  }

  EPS& get_object() {
    return eps;
  }

private:
  int dimension_, offset_local_, num_local_rows_;
  MPI_Comm comm_;
  EPS eps;  // eigenproblem solver context
};

const std::map<std::string,EPSLanczosReorthogType> solver::lanczos_reorthog_map{
  {"local", EPS_LANCZOS_REORTHOG_LOCAL}, {"full", EPS_LANCZOS_REORTHOG_FULL},
  {"selective", EPS_LANCZOS_REORTHOG_SELECTIVE}, {"periodic", EPS_LANCZOS_REORTHOG_PERIODIC},
  {"partial", EPS_LANCZOS_REORTHOG_PARTIAL}, {"delayed", EPS_LANCZOS_REORTHOG_DELAYED} };

} // namespace slepc

} // namespace rokko

#endif // ROKKO_SLEPC_SOLVER_HPP
