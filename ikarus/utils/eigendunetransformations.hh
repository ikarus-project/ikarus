// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file eigendunetransformations.hh
 * \brief Helper for transform between Dune linear algebra types and Eigen
 */

#pragma once

#include <assert.h>

#include <dune/common/diagonalmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/scaledidmatrix.hh>

#include <Eigen/Core>

namespace Ikarus {

  /**
   * @brief Create Eigen::Vector to Dune::FieldVector.
   * @tparam ScalarType The scalar type of the vectors.
   * @tparam size The size of the vectors.
   * @param vec The Eigen::Vector to be converted.
   * @return Dune::FieldVector<ScalarType, size> representing the converted vector.
   */
  template <typename ScalarType, int size>
  Dune::FieldVector<ScalarType, size> toDune(const Eigen::Vector<ScalarType, size>& vec) {
    Dune::FieldVector<ScalarType, size> fieldVector;
    for (int i = 0; i < size; ++i)
      fieldVector[i] = vec[i];
    return fieldVector;
  }

  /**
   * @brief Convert Eigen::Matrix to Dune::FieldVector.
   * \details The matrix has fixed row size. The matrix needs
   * to have a single column.
   * @tparam ScalarType The scalar type of the vectors.
   * @tparam rows The number of rows in the matrix.
   * @param vec The Eigen::Matrix to be converted.
   * @return Dune::FieldVector<ScalarType, rows> representing the converted matrix.
   */
  template <typename ScalarType, int rows>
  Dune::FieldVector<ScalarType, rows> toDune(const Eigen::Matrix<ScalarType, rows, 0>& vec) {
    assert(vec.cols() == 1 && "The passed matrix needs to have a single column.");
    Dune::FieldVector<ScalarType, rows> fieldVector{0.0};

    for (int i = 0; i < vec.rows(); ++i)
      fieldVector[i] = vec(i, 0);
    return fieldVector;
  }

  /**
   * @brief Convert Eigen::Matrix to Dune::FieldMatrix.
   * \details The matrix has fixed rows and column size
   * @tparam ScalarType The scalar type of the matrix.
   * @tparam rows The number of rows in the matrix.
   * @tparam cols The number of columns in the matrix.
   * @param mat The Eigen::Matrix to be converted.
   * @return Dune::FieldMatrix<ScalarType, rows, cols> representing the converted matrix.
   */
  template <typename ScalarType, int rows, int cols>
  Dune::FieldMatrix<ScalarType, rows, cols> toDune(const Eigen::Matrix<ScalarType, rows, cols>& mat) {
    Dune::FieldMatrix<ScalarType, rows, cols> fieldMatrix{0.0};

    for (int i = 0; i < mat.rows(); ++i)
      for (int j = 0; j < mat.cols(); ++j)
        fieldMatrix[i][j] = mat(i, j);
    return fieldMatrix;
  }

  /**
   * @brief View a Dune::FieldVector as an Eigen::Vector using Map, no copies take place.
   * @tparam ScalarType The scalar type of the vector.
   * @tparam size The size of the vector.
   * @param vec The Dune::FieldVector to be viewed as Eigen::Vector.
   * @return Eigen::Map<const Eigen::Vector<ScalarType, size>> representing the viewed vector.
   */
  template <typename ScalarType, int size>
  Eigen::Map<const Eigen::Vector<ScalarType, size>> toEigenMap(const Dune::FieldVector<ScalarType, size>& vec) {
    return {vec.data(), size};
  }

  /**
   * @brief View a constant Dune::FieldVector as a constant Eigen::Vector, no copies take place.
   * @tparam ScalarType The scalar type of the vector.
   * @tparam size The size of the vector.
   * @param vec The Dune::FieldVector to be viewed as Eigen::Vector.
   * @return Eigen::Map<Eigen::Vector<ScalarType, size>> representing the viewed vector.
   */
  template <typename ScalarType, int size>
  Eigen::Map<Eigen::Vector<ScalarType, size>> toEigenMap(Dune::FieldVector<ScalarType, size>& vec) {
    return {vec.data(), size};
  }

}  // namespace Ikarus
