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

#ifndef ROKKO_ANASAZI_DISTRIBUTED_MULTIVECTOR_H
#define ROKKO_ANASAZI_DISTRIBUTED_MULTIVECTOR_H

#include <rokko/mapping_1d.hpp>

#include <Epetra_MultiVector.h>
#include <AnasaziEpetraAdapter.hpp>
#include <Teuchos_RCPDecl.hpp>

namespace rokko {

class distributed_multivector_anasazi {
public:
  distributed_multivector_anasazi() {}
  distributed_multivector_anasazi(mapping_1d const& map, int blocksize) : map_(map) {
    multivector_ = Teuchos::rcp(new Epetra_MultiVector(map_.get_epetra_map(), blocksize));
  }
  ~distributed_multivector_anasazi() {}
  void init_random() {
    multivector_->Random();
  }
  Teuchos::RCP<Epetra_MultiVector> get_pointer() const { return multivector_; }
private:
  mapping_1d map_;
  Teuchos::RCP<Epetra_MultiVector> multivector_;  
};

} // namespace rokko

#endif // ROKKO_ANASAZI_DISTRIBUTED_MULTIVECTOR_H
