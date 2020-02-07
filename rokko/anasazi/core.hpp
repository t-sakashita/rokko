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
  using value_type = double;
  using eigenproblem_t = Anasazi::BasicEigenproblem<value_type, Epetra_MultiVector, Epetra_Operator>;
  using solvermanager_t = Anasazi::SolverManager<value_type, Epetra_MultiVector, Epetra_Operator>;

  static const std::vector<std::string> names;
  
  solver() {
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  }
  ~solver() = default;
  void initialize(int& argc, char**& argv) {}
  void finalize() {}

  static std::unique_ptr<solvermanager_t> create_solver_manager(std::string const& routine, Teuchos::RCP<eigenproblem_t> problem, Teuchos::ParameterList& pl) {
    if ((routine == "LOBPCG") || (routine == ""))
      return std::make_unique<Anasazi::LOBPCGSolMgr<value_type, Epetra_MultiVector, Epetra_Operator>>(problem, pl);
    else if (routine == "BlockKrylovSchur")
      return std::make_unique<Anasazi::BlockKrylovSchurSolMgr<value_type, Epetra_MultiVector, Epetra_Operator>>(problem, pl);
    else if (routine == "BlockDavidson")
      return std::make_unique<Anasazi::BlockDavidsonSolMgr<value_type, Epetra_MultiVector, Epetra_Operator>>(problem, pl);
    else if (routine == "RTR")
      return std::make_unique<Anasazi::RTRSolMgr<value_type, Epetra_MultiVector, Epetra_Operator>>(problem, pl);
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
      pl.set("Convergence Tolerance", params.get<value_type>("conv_tol"));
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
        else if (params.type(key) == typeid(value_type))
          pl.set(key, params.get<value_type>(key));
        else if (params.type(key) == typeid(std::string))
          pl.set(key, params.get<std::string>(key));
        else if (params.type(key) == typeid(const char*))
          pl.set(key, params.get<const char*>(key));
      }
    }

    return pl;
  }

  parameters diagonalize(const rokko::distributed_crs_matrix& mat, rokko::parameters const& params) {
    return diagonalize(*static_cast<const rokko::anasazi::distributed_crs_matrix*>(mat.get_ptr()->get_impl()), params);
  }

  parameters diagonalize(const rokko::anasazi::distributed_crs_matrix& mat, rokko::parameters const& params) {
    int num_eigvals, max_block_size;
    std::tie(num_eigvals, max_block_size) = retrieve_number_size(params);
    Teuchos::ParameterList pl = set_anasazi_parameters(params);

    map_ = mat.get_map();
    Teuchos::RCP<Epetra_MultiVector> multivector = Teuchos::rcp(new Epetra_MultiVector(map_->get_epetra_map(), max_block_size));
    multivector->Random();
    problem_ = Teuchos::rcp(new eigenproblem_t(mat.get_matrix(), multivector));
    problem_->setHermitian(true);
    problem_->setNEV(num_eigvals);
    problem_->setProblem();

    std::unique_ptr<solvermanager_t> solvermanager = create_solver_manager(get_routine(params), problem_, pl);
    
    bool boolret = problem_->setProblem();
    if (!boolret) {
      throw std::invalid_argument("anasazi::diagonalize : Return value from setProblem() is false");
    }

    Anasazi::ReturnType returnCode = solvermanager->solve();
    if (returnCode == Anasazi::Unconverged) {
      std::cout << "solvermanager.solve() does not converge." << std::endl;
    }

    parameters params_out;
    params_out.set("num_conv", get_num_conv());
    return params_out;
  }

  parameters diagonalize(const rokko::distributed_mfree& mat_in, rokko::parameters const& params) {
    rokko::distributed_mfree const*const mat = &mat_in;

    int num_eigvals, max_block_size;
    std::tie(num_eigvals, max_block_size) = retrieve_number_size(params);
    Teuchos::ParameterList pl = set_anasazi_parameters(params);

    map_ = new mapping_1d(mat->get_dim());
    Teuchos::RCP<anasazi_mfree_operator> anasazi_op = Teuchos::rcp(new anasazi_mfree_operator(mat, *map_));
    Teuchos::RCP<Epetra_MultiVector> multivector = Teuchos::rcp(new Epetra_MultiVector(map_->get_epetra_map(), max_block_size));
    multivector->Random();
    problem_ = Teuchos::rcp(new eigenproblem_t(anasazi_op, multivector));
    problem_->setHermitian(true);
    problem_->setNEV(num_eigvals);
    problem_->setProblem();

    std::unique_ptr<solvermanager_t> solvermanager = create_solver_manager(get_routine(params), problem_, pl);
    
    bool boolret = problem_->setProblem();
    if (!boolret) {
      throw std::invalid_argument("anasazi::diagonalize : Return value from setProblem() is false");
    }

    Anasazi::ReturnType returnCode = solvermanager->solve();
    if (returnCode == Anasazi::Unconverged) {
      std::cout << "solvermanager.solve() does not converge." << std::endl;
    }
    parameters params_out;
    params_out.set("num_conv", get_num_conv());
    return params_out;
  }

  value_type eigenvalue(int i) const { return problem_->getSolution().Evals[i].realpart; }

  void eigenvector(int k, std::vector<value_type>& vec) const {
    if (vec.size() < map_->get_num_local_rows()) vec.resize(map_->get_num_local_rows());
    eigenvector(k, vec.data());
  }
  void eigenvector(int k, value_type *const vec) const {
    value_type const*const vec_pt = (*problem_->getSolution().Evecs)[k];
    std::copy(vec_pt, vec_pt + map_->get_num_local_rows(), vec);
  }
  void eigenvector(int k, distributed_vector<value_type>& vec) const {
    vec.initialize(map_->get_dim(), map_->get_epetra_map().MinMyGID(),
                   map_->get_epetra_map().MaxMyGID() + 1);
    eigenvector(k, vec.get_storage());
  }

  int get_num_conv() const {
    return problem_->getSolution().numVecs;
  }

private:
  //std::list<std::string> anasazi_keys = { "Which", "Maximum Iterations", "Convergence Tolerance" };
  int myrank;
  const rokko::anasazi::mapping_1d* map_;
  Teuchos::RCP<eigenproblem_t> problem_;
};

const std::vector<std::string> solver::names{ "LOBPCG", "BlockKrylovSchur", "BlockDavidson", "RTR" };

} // namespace anasazi

} // namespace rokko

#endif // ROKKO_ANASAZI_CORE_HPP
