// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <assert.h>

#include <dune/common/diagonalmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/scaledidmatrix.hh>

#include <Eigen/Core>

namespace Ikarus {

  /** \brief Creates a Dune::FieldVector from a given Eigen::Vector */
  template <typename ScalarType, int size>
  Dune::FieldVector<ScalarType, size> toDune(const Eigen::Vector<ScalarType, size>& vec) {
    Dune::FieldVector<ScalarType, size> fieldVector;
    for (int i = 0; i < size; ++i)
      fieldVector[i] = vec[i];
    return fieldVector;
  }

  /** \brief Creates a Dune::FieldVector from a given Eigen::Matrix. The matrix has fixed dynamic size. The matrix needs
   * to have a single column. */
  template <typename ScalarType, int rows>
  Dune::FieldVector<ScalarType, rows> toDune(const Eigen::Matrix<ScalarType, rows, 0>& vec) {
    assert(vec.cols() == 1 && "The passed matrix needs to have a single column.");
    Dune::FieldVector<ScalarType, rows> fieldVector{0.0};

    for (int i = 0; i < vec.rows(); ++i)
      fieldVector[i] = vec(i, 0);
    return fieldVector;
  }

  /** \brief Creates a Dune::FieldMatrix from a given Eigen::Matrix. The matrix has fixed dynamic size  **/
  template <typename ScalarType, int rows, int cols>
  Dune::FieldMatrix<ScalarType, rows, cols> toDune(const Eigen::Matrix<ScalarType, rows, cols>& mat) {
    Dune::FieldMatrix<ScalarType, rows, cols> fieldMatrix{0.0};

    for (int i = 0; i < mat.rows(); ++i)
      for (int j = 0; j < mat.cols(); ++j)
        fieldMatrix[i][j] = mat(i, j);
    return fieldMatrix;
  }

  /** \brief Views a dune fieldvector as an Eigen::Vector as Map, no copies take place! */
  template <typename ScalarType, int size>
  Eigen::Map<const Eigen::Vector<ScalarType, size>> toEigenMap(const Dune::FieldVector<ScalarType, size>& vec) {
    return {vec.data(), size};
  }

  /** \brief Views a const dune fieldvector as a const Eigen::Vector, no copies take place! */
  template <typename ScalarType, int size>
  Eigen::Map<Eigen::Vector<ScalarType, size>> toEigenMap(Dune::FieldVector<ScalarType, size>& vec) {
    return {vec.data(), size};
  }

}  // namespace Ikarus
