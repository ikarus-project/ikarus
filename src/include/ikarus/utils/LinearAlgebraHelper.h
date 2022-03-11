

#pragma once
#include <dune/istl/bvector.hh>

#include <Eigen/Dense>

namespace Ikarus::LinearAlgebra {
  template <typename Derived>
  auto orthonormalizeMatrixColumns(const Eigen::MatrixBase<Derived>& A) {
    // Gram Schmidt Ortho
    auto Q = A.eval();

    Q.col(0).normalize();

    for (int colIndex = 1; colIndex < Q.cols(); colIndex++) {
      Q.col(colIndex) -= Q.leftCols(colIndex) * (Q.leftCols(colIndex).transpose() * A.col(colIndex));
      Q.col(colIndex).normalize();
    }

    return Q;
  }

  template <typename ValueType>
  auto viewAsFlatEigenVector(Dune::BlockVector<ValueType>& blockedVector) {
    Eigen::Map<Eigen::VectorX<typename ValueType::field_type>> vec(&blockedVector.begin()->begin().operator*(),
                                                                   blockedVector.size() * blockedVector[0].size());

    return vec;
  }
  template <typename ValueType>
  auto viewAsFlatEigenVector(const Dune::BlockVector<ValueType>& blockedVector) {
    Eigen::Map<const Eigen::VectorX<typename ValueType::field_type>> vec(
        &blockedVector.begin()->begin().operator*(), blockedVector.size() * blockedVector[0].size());

    return vec;
  }
}  // namespace Ikarus::LinearAlgebra

template <typename Derived>
requires(!std::floating_point<Derived>) auto norm(const Eigen::MatrixBase<Derived>& v) { return v.norm(); }

/** \brief Helper Free Function to have the same interface as for Eigen Vector Types */
auto norm(const std::floating_point auto& v) { return std::abs(v); }

template <typename Derived>
Derived sym(const Eigen::MatrixBase<Derived>& A) {
  return 0.5 * (A.transpose() + A.transpose());
}