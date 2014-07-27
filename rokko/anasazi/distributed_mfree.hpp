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

#ifndef ROKKO_ANASAZI_DISTRIBUTED_MFREE_H
#define ROKKO_ANASAZI_DISTRIBUTED_MFREE_H

#include <rokko/mapping_1d.hpp>
#include <rokko/distributed_mfree.hpp>

#include <Epetra_CrsMatrix.h>
#include <AnasaziEpetraAdapter.hpp>
#include <Teuchos_RCPDecl.hpp>

namespace rokko {
namespace anasazi {

class distributed_mfree_operator : public Epetra_Operator {
public:
  distributed_mfree_operator(const rokko::distributed_operator& op, mapping_1d const& map) : op_(op), ep_comm(MPI_COMM_WORLD), ep_map(map.get_epetra_map()) {}

  ~distributed_mfree_operator() {};

  virtual int SetUseTranspose(bool UseTranspose) { return 0; };

  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    const int numvectors = X.NumVectors();

    Y.PutScalar(0);
    for (int i=0; i<numvectors; ++i) {
      const double* x = X[i];
      double* y = Y[i];
      op_.multiply(x, y);
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
  const rokko::distributed_operator& op_;
  MPI_Comm comm_;
  mutable std::vector<double> buffer_;
  int L_;
  std::vector<std::pair<int, int> > lattice_;
  Epetra_MpiComm ep_comm;
  Epetra_Map ep_map;
};

class distributed_mfree : public rokko::detail::df_matrix_base {
public:
  typedef rokko::detail::df_matrix_base super_type;
  distributed_mfree() : super_type() {}
  ~distributed_mfree() {}
  void initialize(mapping_1d const& map) {
    map_ = map;
  }
  void define_operator(rokko::distributed_operator& op) {
    op_ = op;
    anasazi_op_ = Teuchos::rcp( new rokko::anasazi::distributed_mfree_operator(op_, map_) );
  }
private:
  mapping_1d map_;
  std::vector<int> rows_;
public:
  rokko::distributed_operator op_;
  Teuchos::RCP<rokko::anasazi::distributed_mfree_operator> anasazi_op_;
};

} // namespace anasazi
} // namespace rokko

#endif // ROKKO_ANASAZI_DISTRIBUTED_MFREE_H
