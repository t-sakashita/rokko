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

#include "distributed_crs_matrix.hpp"
#include "distributed_multivector.hpp"

#include <AnasaziBasicEigenproblem.hpp>
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include <Teuchos_RCPDecl.hpp>

namespace rokko {

class solver_anasazi {
public:
  typedef Anasazi::BasicEigenproblem<double, Epetra_MultiVector, Epetra_Operator> eigenproblem_t;
  typedef Anasazi::BlockKrylovSchurSolMgr<double, Epetra_MultiVector, Epetra_Operator> solvermanager_t;
  solver_anasazi() {}
  ~solver_anasazi() {}
  void diagonalize(distributed_crs_matrix_anasazi const& mat,
                   distributed_multivector_anasazi const& ivec,
                   int num_evals, int block_size, int max_iters, double tol) {
    problem_ = Teuchos::rcp(new eigenproblem_t(mat.get_pointer(), ivec.get_pointer()));
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
  std::vector<Anasazi::Value<double> > eigenvalues() const { return evals_; }
  distributed_multivector_anasazi eigenvectors() const { return evecs_; }
private:
  Teuchos::RCP<eigenproblem_t> problem_;
  std::vector<Anasazi::Value<double> > evals_;
  distributed_multivector_anasazi evecs_;
};

} // namespace rokko

#endif // ROKKO_ANASAZI_DISTRIBUTED_CRS_MATRIX_H
