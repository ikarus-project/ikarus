

#pragma once
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

  //  template <typename Derived> requires (!std::floating_point<Derived>)
  //  auto norm(const Eigen::MatrixBase<Derived>& v) {
  //    return v.norm();
  //  }
  //
  //  /** \brief Helper Free Function to have the same interface as for Eigen Vector Types */
  //  auto norm(const std::floating_point auto& v) {
  //    return std::abs(v);
  //  }

}  // namespace Ikarus::LinearAlgebra

template <typename Derived>
  requires(!std::floating_point<Derived>)
auto norm(const Eigen::MatrixBase<Derived>& v) {
  return v.norm();
}

/** \brief Helper Free Function to have the same interface as for Eigen Vector Types */
auto norm(const std::floating_point auto& v) { return std::abs(v); }

template <typename Derived>
Derived sym(const Eigen::MatrixBase<Derived>& A) {
  return 0.5 * (A.transpose() + A.transpose());
}