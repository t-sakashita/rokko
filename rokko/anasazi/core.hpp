/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2014 Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_ANASAZI_CORE_HPP
#define ROKKO_ANASAZI_CORE_HPP

#include <rokko/distributed_vector.hpp>
#include <rokko/distributed_mfree.hpp>
#include <rokko/anasazi/mapping_1d.hpp>
#include <rokko/anasazi/distributed_crs_matrix.hpp>
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
#include <Teuchos_RCPDecl.hpp>

namespace rokko {

namespace anasazi {

class anasazi_mfree_operator : public Epetra_Operator {
public:
  anasazi_mfree_operator(rokko::distributed_mfree* const op, mapping_1d const* map) :
    op_(op), ep_map(map->get_epetra_map()), ep_comm(map->get_epetra_comm()) {}
  ~anasazi_mfree_operator() {}
  virtual int SetUseTranspose(bool UseTranspose) { return 0; };
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    const int numvectors = X.NumVectors();
    Y.PutScalar(0);
    for (int i=0; i<numvectors; ++i) {
      const double* x = X[i];
      double* y = Y[i];
      op_->multiply(x, y);
    }
    return 0;
  }
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const { return 0; }
  virtual double NormInf() const { return 0; }
  virtual const char * Label() const { return "Anasazi matrix_free"; }
  virtual bool UseTranspose() const { return false; }
  virtual bool HasNormInf() const  { return false; }
  virtual const Epetra_Comm & Comm() const { return ep_comm; }
  virtual const Epetra_Map & OperatorDomainMap() const { return ep_map; }
  virtual const Epetra_Map & OperatorRangeMap() const { return ep_map; }
private:
  distributed_mfree* op_;
  MPI_Comm comm_;
  Epetra_MpiComm ep_comm;
  Epetra_Map ep_map;
};

static const char* const anasazi_solvers[] = { "LOBPCG", "BlockKrylovSchur", "BlockDavidson" };

class solver {
public:
  typedef Anasazi::BasicEigenproblem<double, Epetra_MultiVector, Epetra_Operator> eigenproblem_t;
  typedef Anasazi::LOBPCGSolMgr<double, Epetra_MultiVector, Epetra_Operator> solvermanager_lobpcg_t;
  typedef Anasazi::SolverManager<double, Epetra_MultiVector, Epetra_Operator> solvermanager_t;

  solver() {
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  }
  ~solver() {}
  void initialize(int& argc, char**& argv) {}
  void finalize() {}

  solvermanager_t* create_solver_manager(std::string const& routine) {
    if (routine == "LOBPCG")
      return new Anasazi::LOBPCGSolMgr<double, Epetra_MultiVector, Epetra_Operator>(problem_, pl_);
    if (routine == "BlockKrylovSchur")
      return new Anasazi::BlockKrylovSchurSolMgr<double, Epetra_MultiVector, Epetra_Operator>(problem_, pl_);
    if (routine == "BlockDavidson")
      return new Anasazi::BlockDavidsonSolMgr<double, Epetra_MultiVector, Epetra_Operator>(problem_, pl_);
    else {
      std::stringstream msg;
      msg << "anasazi::solver::create_solver_manager : " << routine << " is not a solver in Anasazi" << std::endl;
      msg << "list of Anasazi solvers:" << std::endl;
      for (int i=0; i<ARRAY_SIZE(anasazi_solvers); ++i) {
	std::cerr << anasazi_solvers[i] << " " << std::endl;
      }
      BOOST_THROW_EXCEPTION(std::invalid_argument(msg.str()));
    }
  }

  void set_anasazi_parameters(rokko::parameters const& params) {
    std::list<std::string> keys = params.keys();
    BOOST_FOREACH(std::string const& key, keys) {
      if (!is_rokko_solver_key(key)) {
	if (params.type(key) == typeid(int)) {
	  pl_.set(key, params.get<int>(key));
	}
	if (params.type(key) == typeid(double)) {
	  pl_.set(key, params.get<double>(key));
	}
	if (params.type(key) == typeid(std::string)) {
	  pl_.set(key, params.get<std::string>(key));
	}
	if (params.type(key) == typeid(const char*)) {
	  pl_.set(key, params.get<const char*>(key));
	}
      }
    }
  }

  parameters diagonalize(rokko::distributed_crs_matrix& mat, rokko::parameters const& params) {
    parameters params_out;
    map_ = new mapping_1d(mat.get_dim());

    set_anasazi_parameters(params);
    if (params.defined("block_size"))  // if block size is provided by common key name "block_size"
      pl_.set("Block Size", params.get<int>("block_size"));
    int num_eigvals;
    if (params.defined("num_eigvals")) num_eigvals = params.get<int>("num_eigvals");
    else num_eigvals = 1;
    int max_block_size;
    if (params.defined("max_block_size")) max_block_size = params.get<int>("max_block_size");
    else max_block_size = num_eigvals;  // fix me : This default value must depend on eigenalgorithm such as 2 * num_eigvals
    if (params.defined("conv_tol"))
      pl_.set("Convergence Tolerance", params.get<double>("conv_tol"));    
    if (params.defined("max_iters"))
      pl_.set("Maximum Iterations", params.get<int>("max_iters"));
    //if (!params.defined("Which")) pl_.set("Which", "LM");

    if (params.get_bool("verbose")) {
      pl_.set( "Verbosity", Anasazi::Errors | Anasazi::Warnings | Anasazi::IterationDetails | Anasazi::FinalSummary | Anasazi::Debug | Anasazi::OrthoDetails );
    }
    multivector_ = Teuchos::rcp(new Epetra_MultiVector(map_->get_epetra_map(), max_block_size));
    multivector_->Random();
    problem_ = Teuchos::rcp(new eigenproblem_t(reinterpret_cast<anasazi::distributed_crs_matrix*>(mat.get_matrix())->get_matrix(), multivector_));
    problem_->setHermitian(true);
    problem_->setNEV(num_eigvals);
    problem_->setProblem();

    if (params.defined("routine")) {
      if ((params.type("routine") != typeid(std::string)) && params.type("routine") != typeid(const char*))
	BOOST_THROW_EXCEPTION(std::invalid_argument("anasazi::solver::diagonalize() : routine must be charatcters or string"));
      routine_ = params.get_string("routine");
    } else {
      routine_ = "LOBPCG";
    }
    solvermanager_t* solvermanager = create_solver_manager(routine_);
    
    bool boolret = problem_->setProblem();
    if (boolret != true) {
      std::cout << "setProblem()_error" << std::endl;
    }

    Anasazi::ReturnType returnCode = solvermanager->solve();
    if (returnCode == Anasazi::Unconverged) {
      std::cout << "solvermanager.solve() does not converge." << std::endl;
    }

    num_conv_ = problem_->getSolution().numVecs;
    params_out.set("num_conv", num_conv_);
    return params_out;
  }

  parameters diagonalize(rokko::distributed_mfree& mat_in, rokko::parameters const& params) {
    parameters params_out;
    rokko::distributed_mfree* mat = &mat_in;
    map_ = new mapping_1d(mat->get_dim());

    set_anasazi_parameters(params);
    if (params.defined("block_size"))  // if block size is provided by common key name "block_size"
      pl_.set("Block Size", params.get<int>("block_size"));
    int num_eigvals;
    if (params.defined("num_eigvals")) num_eigvals = params.get<int>("num_eigvals");
    else num_eigvals = 1;
    int max_block_size;
    if (params.defined("max_block_size")) max_block_size = params.get<int>("max_block_size");
    else max_block_size = num_eigvals;  // fix me : This default value must depend on eigenalgorithm such as 2 * num_eigvals
    if (params.defined("conv_tol"))
      pl_.set("Convergence Tolerance", params.get<double>("conv_tol"));    
    if (params.defined("max_iters"))
      pl_.set("Maximum Iterations", params.get<int>("max_iters"));
    //if (!params.defined("Which")) pl_.set("Which", "LM");

    if (params.get_bool("verbose")) {
      pl_.set( "Verbosity", Anasazi::Errors | Anasazi::Warnings | Anasazi::IterationDetails | Anasazi::FinalSummary | Anasazi::Debug | Anasazi::OrthoDetails );
    }
    
    Teuchos::RCP<anasazi_mfree_operator> anasazi_op_ = Teuchos::rcp(new anasazi_mfree_operator(mat, map_));
    multivector_ = Teuchos::rcp(new Epetra_MultiVector(map_->get_epetra_map(), max_block_size));
    multivector_->Random();
    problem_ = Teuchos::rcp(new eigenproblem_t(anasazi_op_, multivector_));
    problem_->setHermitian(true);
    problem_->setNEV(num_eigvals);
    problem_->setProblem();

    if (params.defined("routine")) {
      if ((params.type("routine") == typeid(std::string)) && params.type("routine") == typeid(const char*))
	BOOST_THROW_EXCEPTION(std::invalid_argument("anasazi::solver::diagonalize() : routine must be charatcters or string"));
      routine_ = params.get_string("routine");
    } else {
      routine_ = "LOBPCG";
    }
    solvermanager_t* solvermanager = create_solver_manager(routine_);
    
    bool boolret = problem_->setProblem();
    if (boolret != true) {
      std::cout << "setProblem()_error" << std::endl;
    }

    Anasazi::ReturnType returnCode = solvermanager->solve();
    if (returnCode == Anasazi::Unconverged) {
      std::cout << "solvermanager.solve() does not converge." << std::endl;
    }
    num_conv_ = problem_->getSolution().numVecs;
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
    eigenvector(k, &(vec[0]));
  }
  void eigenvector(int k, double* vec) const {
    double* vec_pt = (*problem_->getSolution().Evecs)[k];
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
  //std::vector<Anasazi::Value<double> > evals_;
  Teuchos::ParameterList pl_;
  int max_block_size;
  std::string routine_;
  int num_conv_;
};

} // namespace anasazi

} // namespace rokko

#endif // ROKKO_ANASAZI_CORE_HPP
