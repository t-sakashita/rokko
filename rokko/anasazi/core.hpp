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

#ifndef ROKKO_ANASAZI_CORE_HPP
#define ROKKO_ANASAZI_CORE_HPP

#include <tuple>

#include <rokko/distributed_vector.hpp>
#include <rokko/distributed_mfree.hpp>
#include <rokko/anasazi/mapping_1d.hpp>
#include <rokko/anasazi/distributed_crs_matrix.hpp>
#include <rokko/anasazi/distributed_mfree.hpp>
#include <rokko/parameters.hpp>
#include <rokko/solver_parameters.hpp>

#include <AnasaziEpetraAdapter.hpp>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>

#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziLOBPCGSolMgr.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <AnasaziBlockDavidsonSolMgr.hpp>
#include <AnasaziRTRSolMgr.hpp>
#include <Teuchos_RCPDecl.hpp>

namespace rokko {

namespace anasazi {

class solver {
public:
  using eigenproblem_t = Anasazi::BasicEigenproblem<double, Epetra_MultiVector, Epetra_Operator>;
  using solvermanager_t = Anasazi::SolverManager<double, Epetra_MultiVector, Epetra_Operator>;

  static const std::vector<std::string> names;
  
  solver() {
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  }
  ~solver() {}
  void initialize(int& argc, char**& argv) {}
  void finalize() {}

  solvermanager_t* create_solver_manager(std::string const& routine, Teuchos::ParameterList& pl) {
    if ((routine == "LOBPCG") || (routine == ""))
      return new Anasazi::LOBPCGSolMgr<double, Epetra_MultiVector, Epetra_Operator>(problem_, pl);
    else if (routine == "BlockKrylovSchur")
      return new Anasazi::BlockKrylovSchurSolMgr<double, Epetra_MultiVector, Epetra_Operator>(problem_, pl);
    else if (routine == "BlockDavidson")
      return new Anasazi::BlockDavidsonSolMgr<double, Epetra_MultiVector, Epetra_Operator>(problem_, pl);
    else if (routine == "RTR")
      return new Anasazi::RTRSolMgr<double, Epetra_MultiVector, Epetra_Operator>(problem_, pl);
    else {
      std::stringstream msg;
      msg << "anasazi::solver::create_solver_manager : " << routine << " is not a solver in Anasazi" << std::endl;
      msg << "list of Anasazi solvers:" << std::endl;
      for (int i=0; i<names.size(); ++i) {
        msg << names[i] << " " << std::endl;
      }
      throw std::invalid_argument(msg.str());
    }
  }

  static std::tuple<int,int> retrieve_number_size(rokko::parameters const& params) {
    int num_eigvals = params.defined("num_eigvals") ?
      params.get<int>("num_eigvals") : 1;
    int max_block_size = params.defined("max_block_size") ?
      params.get<int>("max_block_size") : num_eigvals;  // fix me : This default value must depend on eigenalgorithm such as 2 * num_eigvals
    return {num_eigvals, max_block_size};
  }

  static std::string get_routine(rokko::parameters const& params) {
    std::string routine;

    if (params.defined("routine")) {
      if ((params.type("routine") == typeid(std::string)) && params.type("routine") == typeid(const char*))
        throw std::invalid_argument("anasazi::solver::diagonalize() : routine must be charatcters or string");
      routine = params.get_string("routine");
    } else {  // default
      routine = "LOBPCG";
    }
    return routine;
  }

  static Teuchos::ParameterList set_anasazi_parameters(rokko::parameters const& params) {
    Teuchos::ParameterList pl;

    if (params.defined("block_size"))  // if block size is provided by common key name "block_size"
      pl.set("Block Size", params.get<int>("block_size"));
    if (params.defined("conv_tol"))
      pl.set("Convergence Tolerance", params.get<double>("conv_tol"));
    if (params.defined("max_iters"))
      pl.set("Maximum Iterations", params.get<int>("max_iters"));
    //if (!params.defined("Which")) pl_.set("Which", "LM");

    if (params.get_bool("verbose")) {
      pl.set( "Verbosity", Anasazi::Errors | Anasazi::Warnings | Anasazi::IterationDetails | Anasazi::FinalSummary | Anasazi::Debug | Anasazi::OrthoDetails );
    }

    std::list<std::string> keys = params.keys();
    for(auto const& key : keys) {
      if (!is_rokko_solver_key(key)) {
        if (params.type(key) == typeid(int))
          pl.set(key, params.get<int>(key));
        else if (params.type(key) == typeid(double))
          pl.set(key, params.get<double>(key));
        else if (params.type(key) == typeid(std::string))
          pl.set(key, params.get<std::string>(key));
        else if (params.type(key) == typeid(const char*))
          pl.set(key, params.get<const char*>(key));
      }
    }

    return pl;
  }

  parameters diagonalize(rokko::distributed_crs_matrix& mat, rokko::parameters const& params) {
    int num_eigvals, max_block_size;
    std::tie(num_eigvals, max_block_size) = retrieve_number_size(params);
    Teuchos::ParameterList pl = set_anasazi_parameters(params);

    map_ = new mapping_1d(mat.get_dim());
    multivector_ = Teuchos::rcp(new Epetra_MultiVector(map_->get_epetra_map(), max_block_size));
    multivector_->Random();
    problem_ = Teuchos::rcp(new eigenproblem_t(reinterpret_cast<anasazi::distributed_crs_matrix*>(mat.get_matrix())->get_matrix(), multivector_));
    problem_->setHermitian(true);
    problem_->setNEV(num_eigvals);
    problem_->setProblem();

    routine_ = get_routine(params);
    solvermanager_t* solvermanager = create_solver_manager(routine_, pl);
    
    bool boolret = problem_->setProblem();
    if (!boolret) {
      std::cerr << "setProblem()_error" << std::endl;
    }

    Anasazi::ReturnType returnCode = solvermanager->solve();
    if (returnCode == Anasazi::Unconverged) {
      std::cout << "solvermanager.solve() does not converge." << std::endl;
    }

    num_conv_ = problem_->getSolution().numVecs;
    parameters params_out;
    params_out.set("num_conv", num_conv_);
    return params_out;
  }

  parameters diagonalize(rokko::distributed_mfree& mat_in, rokko::parameters const& params) {
    rokko::distributed_mfree const*const mat = &mat_in;

    int num_eigvals, max_block_size;
    std::tie(num_eigvals, max_block_size) = retrieve_number_size(params);
    Teuchos::ParameterList pl = set_anasazi_parameters(params);

    map_ = new mapping_1d(mat->get_dim());
    Teuchos::RCP<anasazi_mfree_operator> anasazi_op_ = Teuchos::rcp(new anasazi_mfree_operator(mat, *map_));
    multivector_ = Teuchos::rcp(new Epetra_MultiVector(map_->get_epetra_map(), max_block_size));
    multivector_->Random();
    problem_ = Teuchos::rcp(new eigenproblem_t(anasazi_op_, multivector_));
    problem_->setHermitian(true);
    problem_->setNEV(num_eigvals);
    problem_->setProblem();

    routine_ = get_routine(params);
    solvermanager_t* solvermanager = create_solver_manager(routine_, pl);
    
    bool boolret = problem_->setProblem();
    if (!boolret) {
      std::cerr << "setProblem()_error" << std::endl;
    }

    Anasazi::ReturnType returnCode = solvermanager->solve();
    if (returnCode == Anasazi::Unconverged) {
      std::cout << "solvermanager.solve() does not converge." << std::endl;
    }
    num_conv_ = problem_->getSolution().numVecs;
    parameters params_out;
    params_out.set("num_conv", num_conv_);
    return params_out;
  }

  rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(int row_dim,
    int col_dim) {
    return new anasazi::distributed_crs_matrix(row_dim, col_dim);
  }
  rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(int row_dim,
    int col_dim, int num_entries_per_row) {
    return new anasazi::distributed_crs_matrix(row_dim, col_dim, num_entries_per_row);
  }
  double eigenvalue(int i) const { return problem_->getSolution().Evals[i].realpart; }

  void eigenvector(int k, std::vector<double>& vec) const {
    if (vec.size() < map_->get_num_local_rows()) vec.resize(map_->get_num_local_rows());
    eigenvector(k, vec.data());
  }
  void eigenvector(int k, double *const vec) const {
    double const*const vec_pt = (*problem_->getSolution().Evecs)[k];
    std::copy(vec_pt, vec_pt + map_->get_num_local_rows(), vec);
  }
  void eigenvector(int k, distributed_vector<double>& vec) const {
    vec.initialize(map_->get_dim(), map_->get_epetra_map().MinMyGID(),
                   map_->get_epetra_map().MaxMyGID() + 1);
    eigenvector(k, vec.get_storage());
  }

  int num_conv() const { return num_conv_; }

private:
  //std::list<std::string> anasazi_keys = { "Which", "Maximum Iterations", "Convergence Tolerance" };
  int myrank;
  mapping_1d* map_;
  Teuchos::RCP<Epetra_MultiVector> multivector_;
  Teuchos::RCP<eigenproblem_t> problem_;
  //std::vector<Anasazi::Value<double>> evals_;
  std::string routine_;
  int num_conv_;
};

const std::vector<std::string> solver::names{ "LOBPCG", "BlockKrylovSchur", "BlockDavidson", "RTR" };

} // namespace anasazi

} // namespace rokko

#endif // ROKKO_ANASAZI_CORE_HPP
