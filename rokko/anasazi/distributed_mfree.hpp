/*****************************************************************************
*
* Rokko: Integrated Interface for libraries of eigenvalue decomposition
*
* Copyright (C) 2012-2019 by Rokko Developers https://github.com/t-sakashita/rokko
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef ROKKO_ANASAZI_DISTRIBUTED_MFREE_HPP
#define ROKKO_ANASAZI_DISTRIBUTED_MFREE_HPP

#include <rokko/distributed_mfree.hpp>
#include <rokko/anasazi/mapping_1d.hpp>

#include <AnasaziEpetraAdapter.hpp>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>

namespace rokko {

namespace anasazi {

class anasazi_mfree_operator : public Epetra_Operator {
public:
  anasazi_mfree_operator(rokko::distributed_mfree const& op, mapping_1d const& map) :
    op_(op), ep_map(map.get_epetra_map()), ep_comm(map.get_epetra_comm()) {}
  ~anasazi_mfree_operator() = default;
  virtual int SetUseTranspose(bool UseTranspose) override { return 0; };
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override {
    const int numvectors = X.NumVectors();
    Y.PutScalar(0);
    for (int i=0; i<numvectors; ++i) {
      const double* x = X[i];
      double* y = Y[i];
      op_.multiply(x, y);
    }
    return 0;
  }
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override { return 0; }
  virtual double NormInf() const override { return 0; }
  virtual const char * Label() const override { return "Anasazi matrix_free"; }
  virtual bool UseTranspose() const override { return false; }
  virtual bool HasNormInf() const override { return false; }
  virtual const Epetra_Comm & Comm() const override { return ep_comm; }
  virtual const Epetra_Map & OperatorDomainMap() const override { return ep_map; }
  virtual const Epetra_Map & OperatorRangeMap() const override { return ep_map; }

private:
  const distributed_mfree& op_;
  Epetra_MpiComm ep_comm;
  Epetra_Map ep_map;
};

} // end namespace anasazi

} // end namespace rokko

#endif // ROKKO_ANASAZI_DISTRIBUTED_MFREE_HPP
