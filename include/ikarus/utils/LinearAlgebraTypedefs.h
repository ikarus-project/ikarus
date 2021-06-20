//
// Created by Alex on 10.05.2021.
//

#pragma once

#include <dune/common/fvector.hh>

#include <Eigen/Core>

namespace Ikarus {

  ///** \brief Dynamic matrix */
  template <typename ScalarType> using DynMatrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;

  ///** \brief Fixed size matrix */
  template <typename ScalarType, int rowSize, int colSize> using FixedMatrix
      = Eigen::Matrix<ScalarType, rowSize, colSize>;

  ///** \brief Dynamic size vector */
  template <typename ScalarType> using DynVector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

  ///** \brief Fixed size vector */
  template <typename ScalarType, int size> using FixedVector = Eigen::Matrix<ScalarType, size, 1>;

  ///** \brief Dynamic size eigen array */
  using DynArrayXi = Eigen::ArrayXi;

  ///** \brief Fixed size eigen array */
  template <typename ScalarType, int rowSize, int colSize> using FixedArray
      = Eigen::Array<ScalarType, rowSize, colSize>;

  ///** \brief Fixed size vector with doubles */
  template <int size> using FixedVectord = FixedVector<double, size>;

  ///** \brief Fixed size vector with int components */
  template <int size> using FixedVectori = FixedVector<int, size>;

  ///** \brief Fixed size matrix with doubles */
  template <int rowSize, int colSize> using FixedMatrixd = FixedMatrix<double, rowSize, colSize>;

  ///** \brief Fixed size matrix with int components */
  template <int rowSize, int colSize> using FixedMatrixi = FixedMatrix<int, rowSize, colSize>;

  ///** \brief Partial fixed size double  matrix with \f$rowSize \times n\f$ components */
  template <int rowSize> using FixedRowHybridMatrixd = FixedMatrixd<rowSize, Eigen::Dynamic>;

  ///** \brief Partial fixed size double  matrix with \colSize \times n\f$ components */
  template <int colSize> using FixedColHybridMatrixd = FixedMatrixd<Eigen::Dynamic, colSize>;

  ///** \brief Partial fixed size double  matrix with \f$2 \times n\f$ components */
  using HybridMatrix2Xd = FixedRowHybridMatrixd<2>;
  ///** \brief Partial fixed size double  matrix with \f$3 \times n\f$ components */
  using HybridMatrix3Xd = FixedRowHybridMatrixd<3>;
  ///** \brief Partial fixed size double  matrix with \f$4 \times n\f$ components */
  using HybridMatrix4Xd = FixedRowHybridMatrixd<4>;
  ///** \brief Partial fixed size double  matrix with \f$n \times 2\f$ components */
  using HybridMatrixX2d = FixedColHybridMatrixd<2>;
  ///** \brief Partial fixed size double  matrix with \f$n \times 3\f$ components */
  using HybridMatrixX3d = FixedColHybridMatrixd<3>;
  ///** \brief Partial fixed size double  matrix with \f$n \times 4\f$ components */
  using HybridMatrixX4d = FixedColHybridMatrixd<4>;

  ///** \brief Partial fixed size int matrix with \f$2 \times n\f$ components */
  using HybridMatrix2Xi = FixedMatrixi<2, Eigen::Dynamic>;
  ///** \brief Partial fixed size int matrix with \f$3 \times n\f$ components */
  using HybridMatrix3Xi = FixedMatrixi<3, Eigen::Dynamic>;
  ///** \brief Partial fixed size int matrix with \f$4 \times n\f$ components */
  using HybridMatrix4Xi = FixedMatrixi<4, Eigen::Dynamic>;
  ///** \brief Partial fixed size int matrix with \f$n \times 2\f$ components */
  using HybridMatrixX2i = FixedMatrixi<Eigen::Dynamic, 2>;
  ///** \brief Partial fixed size int matrix with \f$n \times 3\f$ components */
  using HybridMatrixX3i = FixedMatrixi<Eigen::Dynamic, 3>;
  ///** \brief Partial fixed size int matrix with \f$n \times 4\f$ components */
  using HybridMatrixX4i = FixedMatrixi<Eigen::Dynamic, 4>;

  ///** \brief Dynamic size matrix with doubles */
  using DynMatrixd = DynMatrix<double>;
  ///** \brief Dynamic size vector with doubles */
  using DynVectord = DynVector<double>;

  ///** \brief Fixed size \f$2 \times 1\f$ vector  with doubles */
  using FixedVector2d = FixedVectord<2>;
  ///** \brief Fixed size \f$3 \times 1\f$ vector  with doubles */
  using FixedVector3d = FixedVectord<3>;
  ///** \brief Fixed size \f$4 \times 1\f$ vector  with doubles */
  using FixedVector4d = FixedVectord<4>;

  ///** \brief Fixed size \f$2 \times 1\f$ vector  with ints */
  using FixedVector2i = FixedVectori<2>;
  ///** \brief Fixed size \f$3 \times 1\f$ vector  with ints */
  using FixedVector3i = FixedVectori<3>;
  ///** \brief Fixed size \f$4 \times 1\f$ vector  with ints */
  using FixedVector4i = FixedVectori<4>;

  ///** \brief Fixed size \f$2 \times 2\f$ matrix  with doubles */
  using FixedMatrix2d = FixedMatrixd<2, 2>;
  ///** \brief Fixed size \f$3 \times 3\f$ matrix  with doubles */
  using FixedMatrix3d = FixedMatrixd<3, 3>;
  ///** \brief Fixed size \f$4 \times 4\f$ matrix  with doubles */
  using FixedMatrix4d = FixedMatrixd<4, 4>;

  ///** \brief Fixed size \f$2 \times 2\f$ matrix  with ints */
  using FixedMatrix2i = FixedMatrixi<2, 2>;
  ///** \brief Fixed size \f$3 \times 3\f$ matrix  with ints */
  using FixedMatrix3i = FixedMatrixi<3, 3>;
  ///** \brief Fixed size \f$4 \times 4\f$ matrix  with ints */
  using FixedMatrix4i = FixedMatrixi<4, 4>;

  constexpr int Dynamic = Eigen::Dynamic;  // This is just -1 to indicate dynmic sized arrays

  template <typename ScalarType, int size>
  Dune::FieldVector<ScalarType, size> toFieldVector(const FixedVector<ScalarType, size>& vec) {
    Dune::FieldVector<ScalarType, size> fieldvec;
    for (int i = 0; i < size; ++i) fieldvec[i] = vec[i];
    return fieldvec;
  }

}  // namespace Ikarus
