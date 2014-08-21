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
//#include <rokko/distributed_crs_matrix.hpp>

#include <Epetra_CrsMatrix.h>
#include <AnasaziEpetraAdapter.hpp>
#include <Teuchos_RCPDecl.hpp>

namespace rokko {
namespace anasazi {


class distributed_crs_matrix {
public:
  distributed_crs_matrix(rokko::anasazi::solver& solver_in) {}
  ~distributed_crs_matrix() {}
  void initialize(mapping_1d const& map, ) {
    map_ = map;
    matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_.get_epetra_map(), map_.dimension()));
  }
  void insert(int row, std::vector<int> const& cols, std::vector<double> const& values) {
    matrix_->InsertGlobalValues(row, cols.size(), &values[0], &cols[0]);
  }
  void complete() {
    matrix_->FillComplete();
    matrix_->SetTracebackMode(1);
  }
private:
  mapping_1d map_;
  std::vector<int> rows_;
public:
  Teuchos::RCP<Epetra_CrsMatrix> matrix_;  
};

} // namespace anasazi
} // namespace rokko

#endif // ROKKO_ANASAZI_DISTRIBUTED_CRS_MATRIX_H
