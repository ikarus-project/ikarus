//
// Created by Alex on 10.05.2021.
//

#pragma once

#include <dune/common/fvector.hh>

#include <Eigen/Core>

namespace Ikarus {

  template <typename ScalarType, int size>
  Dune::FieldVector<ScalarType, size> toFieldVector(const Eigen::Vector<ScalarType, size>& vec) {
    Dune::FieldVector<ScalarType, size> fieldvec;
    for (int i = 0; i < size; ++i)
      fieldvec[i] = vec[i];
    return fieldvec;
  }

}  // namespace Ikarus
