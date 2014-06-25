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

#ifndef ROKKO_ANASAZI_DISTRIBUTED_CRS_MATRIX_H
#define ROKKO_ANASAZI_DISTRIBUTED_CRS_MATRIX_H

#include <rokko/mapping_1d.hpp>

#include <Epetra_CrsMatrix.h>
#include <AnasaziEpetraAdapter.hpp>
#include <Teuchos_RCPDecl.hpp>

namespace rokko {

class distributed_crs_matrix_anasazi {
public:
  distributed_crs_matrix_anasazi(mapping_1d const& map) : map_(map) {
    matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_.get_epetra_map(), map_.dimension()));
  }
  ~distributed_crs_matrix_anasazi() {}
  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) {
    matrix_->InsertGlobalValues(row, cols.size(), &values[0], &cols[0]);
  }
  void finish_insert() {
    matrix_->FillComplete();
    matrix_->SetTracebackMode(1);
  }
  Teuchos::RCP<Epetra_CrsMatrix> get_pointer() const { return matrix_; }
private:
  mapping_1d const& map_;
  std::vector<int> rows_;
  Teuchos::RCP<Epetra_CrsMatrix> matrix_;  
};

} // namespace rokko

#endif // ROKKO_ANASAZI_DISTRIBUTED_CRS_MATRIX_H
