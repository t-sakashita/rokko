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

#include <rokko/anasazi/distributed_crs_matrix.hpp>
#include <rokko/distributed_mfree.hpp>
//#include "distributed_crs_matrix.hpp"
#include "distributed_multivector.hpp"
#include <rokko/utility/timer.hpp>

#include <AnasaziEpetraAdapter.hpp>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>

#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <Teuchos_RCPDecl.hpp>

namespace rokko {

namespace anasazi {

class anasazi_mfree_operator : public Epetra_Operator {
public:
  anasazi_mfree_operator(rokko::distributed_mfree* op) : op_(op), ep_comm(Epetra_MpiComm(MPI_COMM_WORLD)), ep_map(op->get_mapping_1d().get_epetra_map()) {}

  ~anasazi_mfree_operator() {};

  virtual int SetUseTranspose(bool UseTranspose) { return 0; };

  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    const int numvectors = X.NumVectors();

    Y.PutScalar(0);
    for (int i=0; i<numvectors; ++i) {
      const double* x = X[i];
      double* y = Y[i];
      op_->multiply(x, y);
    }
    //std::cout << "X=" << X << std::endl;
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
  mutable std::vector<double> buffer_;
  int L_;
  std::vector<std::pair<int, int> > lattice_;
  Epetra_MpiComm ep_comm;
  Epetra_Map ep_map;
};

class solver {
public:
  typedef Anasazi::BasicEigenproblem<double, Epetra_MultiVector, Epetra_Operator> eigenproblem_t;
  typedef Anasazi::BlockKrylovSchurSolMgr<double, Epetra_MultiVector, Epetra_Operator> solvermanager_t;
  solver() {}
  ~solver() {}
  void initialize(int& argc, char**& argv) {}
  void finalize() {}
  void diagonalize(distributed_crs_matrix const& mat,
                   distributed_multivector_anasazi const& ivec,
                   int num_evals, int block_size, int max_iters, double tol, timer& time) {
    problem_ = Teuchos::rcp(new eigenproblem_t(reinterpret_cast<anasazi::distributed_crs_matrix*>(mat->get_matrix().get_matrix(), ivec.get_pointer()));
    problem_->setHermitian(true);
    problem_->setNEV(num_evals);
    problem_->setProblem();

    Teuchos::ParameterList pl;
    pl.set("Which", "SR");
    pl.set("Block Size", block_size);
    pl.set("Maximum Iterations", max_iters);
    pl.set("Convergence Tolerance", tol);
    solvermanager_t solvermanager(problem_, pl);
    solvermanager.solve();

    Anasazi::Eigensolution<double, Epetra_MultiVector> sol = problem_->getSolution();
    evals_ = sol.Evals;
    evecs_.get_pointer() = sol.Evecs;
  }

  void diagonalize(rokko::distributed_mfree* const mat,
                   distributed_multivector_anasazi const& ivec,
                   int num_evals, int block_size, int max_iters, double tol) {
    Teuchos::RCP<anasazi_mfree_operator> anasazi_op_ = Teuchos::rcp( new anasazi_mfree_operator(mat) );
    problem_ = Teuchos::rcp(new eigenproblem_t(anasazi_op_, ivec.get_pointer()));
    problem_->setHermitian(true);
    problem_->setNEV(num_evals);
    problem_->setProblem();

    Teuchos::ParameterList pl;
    pl.set("Which", "SR");
    pl.set("Block Size", block_size);
    pl.set("Maximum Iterations", max_iters);
    pl.set("Convergence Tolerance", tol);
    solvermanager_t solvermanager(problem_, pl);
    bool boolret = problem_->setProblem();
    if (boolret != true) {
      std::cout << "setProblem()_error" << std::endl;
    }
    Anasazi::ReturnType returnCode = solvermanager.solve();
    if (returnCode == Anasazi::Unconverged) {
      std::cout << "solvermanager.solve()_error" << std::endl;
    }

    Anasazi::Eigensolution<double, Epetra_MultiVector> sol = problem_->getSolution();
    evals_ = sol.Evals;
    evecs_.get_pointer() = sol.Evecs;
  }

  rokko::detail::distributed_crs_matrix_base* create_distributed_crs_matrix(mapping_1d const& map) {
    return new anasazi::distributed_crs_matrix(map);
  }

  std::vector<Anasazi::Value<double> > eigenvalues() const { return evals_; }
  distributed_multivector_anasazi eigenvectors() const { return evecs_; }
private:
  Teuchos::RCP<eigenproblem_t> problem_;
  std::vector<Anasazi::Value<double> > evals_;
  distributed_multivector_anasazi evecs_;
};

} // namespace anasazi

} // namespace rokko

#endif // ROKKO_ANASAZI_DISTRIBUTED_CRS_MATRIX_H
