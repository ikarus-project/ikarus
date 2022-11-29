// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include <dune/common/diagonalmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <Eigen/Core>

namespace Ikarus {

  /** \brief Creates a Dune::FieldVector from a given Eigen::Vector */
  template <typename ScalarType, int size>
  Dune::FieldVector<ScalarType, size> toDune(const Eigen::Vector<ScalarType, size>& vec) {
    Dune::FieldVector<ScalarType, size> fieldvec;
    for (int i = 0; i < size; ++i)
      fieldvec[i] = vec[i];
    return fieldvec;
  }

  /** \brief Creates a Dune::FieldVector from a given Eigen::Matrix. The matrix has fixed dynamic size. The matrix needs
   * to have a single column. */
  template <typename ScalarType, int maxRows, int maxCols>
  Dune::FieldVector<ScalarType, maxRows> toDune(
      const Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, 0, maxRows, maxCols>& vec) {
    assert(vec.cols() == 1 && "The passed matrix needs to have a single column.");
    Dune::FieldVector<ScalarType, maxRows> fieldvec{0.0};

    for (int i = 0; i < vec.rows(); ++i)
      fieldvec[i] = vec(i, 0);
    return fieldvec;
  }

  /** \brief Creates a Dune::FieldMatrix from a given Eigen::Matrix. The matrix has fixed dynamic size  **/
  template <typename ScalarType, int maxRows, int maxCols>
  Dune::FieldMatrix<ScalarType, maxRows, maxCols> toDune(
      const Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic, 0, maxRows, maxCols>& mat) {
    Dune::FieldMatrix<ScalarType, maxRows, maxCols> fieldmat{0.0};

    for (int i = 0; i < mat.rows(); ++i)
      for (int j = 0; j < mat.cols(); ++j)
        fieldmat[i][j] = mat(i, j);
    return fieldmat;
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

  /** \brief Creates a Eigen::Vector from a given Dune::FieldVector  */
  template <typename ScalarType, int size>
  Eigen::Vector<ScalarType, size> toEigen(const Dune::FieldVector<ScalarType, size>& vec) {
    Eigen::Vector<ScalarType, size> eigenVector;
    for (int i = 0; i < size; ++i)
      eigenVector(i) = vec[i];
    return eigenVector;
  }

  /** \brief Creates a Eigen::Matrix from a given Dune::FieldMatrix  */
  template <typename ScalarType, int size1, int size2>
  Eigen::Matrix<ScalarType, size1, size2> toEigen(const Dune::FieldMatrix<ScalarType, size1, size2>& mat) {
    Eigen::Matrix<ScalarType, size1, size2> eigenmatrix;
    for (int i = 0; i < size1; ++i)
      for (int j = 0; j < size2; ++j)
        eigenmatrix(i, j) = mat[i][j];
    return eigenmatrix;
  }

  /** \brief Creates a Eigen::Matrix from a given Dune::DiagonalMatrix. This should return Eigen::DiagonalMatrix but
   * Eigen::DiagonalMatrix does not contain e.g. a transpose method. And therefore we would need to specialize user
   * code. Maybe someone wants to do a PR at Eigen? */
  template <typename ScalarType, int size1>
  Eigen::Matrix<ScalarType, size1, size1> toEigen(const Dune::DiagonalMatrix<ScalarType, size1>& mat) {
    Eigen::Matrix<ScalarType, size1, size1> eigenmatrix;
    eigenmatrix.setZero();
    for (int i = 0; i < size1; ++i)
      eigenmatrix(i, i) = mat[i][i];
    return eigenmatrix;
  }

}  // namespace Ikarus
