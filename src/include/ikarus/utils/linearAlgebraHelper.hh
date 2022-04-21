

#pragma once
#include <dune/istl/bvector.hh>

#include <Eigen/Dense>
#include <ikarus/manifolds/manifoldInterface.hh>


namespace Ikarus::LinearAlgebra {

  /** \brief Orthonormalizes all Matrix columns */
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

  /** \brief View Dune::BlockVector as a Eigen::Vector */
  template <typename ValueType>
  auto viewAsFlatEigenVector(Dune::BlockVector<ValueType>& blockedVector) {
    Eigen::Map<Eigen::VectorX<typename ValueType::field_type>> vec(&blockedVector.begin()->begin().operator*(),
                                                                   blockedVector.size() * blockedVector[0].size());

    return vec;
  }

  /** \brief View Dune::BlockVector as a Eigen::Matrix with dynamic rows and fixed columns depending on the size of  the
   * ValueType*/
  template <typename ValueType>
  auto viewAsEigenMatrixAsDynFixed(Dune::BlockVector<ValueType>& blockedVector) {
    Eigen::Map<Eigen::Matrix<typename ValueType::field_type, Eigen::Dynamic, ValueType::valueSize, Eigen::RowMajor>>
        vec(&blockedVector.begin()->begin().operator*(), blockedVector.size(), blockedVector[0].size());

    return vec;
  }

  /** \brief Const view Dune::BlockVector as a Eigen::Matrix with dynamic rows and fixed columns depending on the size
   * of  the ValueType */
  template <typename ValueType>
  auto viewAsEigenMatrixAsDynFixed(const Dune::BlockVector<ValueType>& blockedVector) {
    Eigen::Map<
        const Eigen::Matrix<typename ValueType::field_type, Eigen::Dynamic, ValueType::valueSize, Eigen::RowMajor>>
        vec(&blockedVector.begin()->begin().operator*(), blockedVector.size(), blockedVector[0].size());

    return vec;
  }

  /** \brief View Dune::BlockVector as a Eigen::Matrix with fixed rows depending on the size of  the ValueType and
   * dynamics columns */
  template <typename ValueType>
  auto viewAsEigenMatrixFixedDyn(Dune::BlockVector<ValueType>& blockedVector) {
    Eigen::Map<Eigen::Matrix<typename ValueType::field_type, ValueType::valueSize, Eigen::Dynamic>> vec(
        &blockedVector.begin()->begin().operator*(), blockedVector[0].size(), blockedVector.size());

    return vec;
  }

  /** \brief Const view Dune::BlockVector as a Eigen::Matrix with fixed rows depending on the size of  the ValueType and
   * dynamics columns */
  template <typename ValueType>
  auto viewAsEigenMatrixFixedDyn(const Dune::BlockVector<ValueType>& blockedVector) {
    Eigen::Map<const Eigen::Matrix<typename ValueType::field_type, ValueType::valueSize, Eigen::Dynamic>> vec(
        &blockedVector.begin()->begin().operator*(), blockedVector[0].size(), blockedVector.size());

    return vec;
  }

  /** \brief View Dune::BlockVector as a Eigen::Vector */
  template <typename ValueType>
  auto viewAsFlatEigenVector(const Dune::BlockVector<ValueType>& blockedVector) {
    Eigen::Map<const Eigen::VectorX<typename ValueType::field_type>> vec(
        &blockedVector.begin()->begin().operator*(), blockedVector.size() * blockedVector[0].size());

    return vec;
  }
}  // namespace Ikarus::LinearAlgebra

/** \brief Adding free norm function to Eigen types */
template <typename Derived>
requires(!std::floating_point<Derived>) auto norm(const Eigen::MatrixBase<Derived>& v) { return v.norm(); }

/** \brief Helper Free Function to have the same interface as for Eigen Vector Types */
auto norm(const std::floating_point auto& v) { return std::abs(v); }

/** \brief Returns the symmetric part of a matrix*/
template <typename Derived>
Derived sym(const Eigen::MatrixBase<Derived>& A) {
  return 0.5 * (A + A.transpose());
}

/** \brief Returns the skew part of a matrix*/
template <typename Derived>
Derived skew(const Eigen::MatrixBase<Derived>& A) {
  return 0.5 * (A - A.transpose());
}

/** \brief Evaluates Eigen expressions */
template <typename Derived>
auto eval(const Eigen::EigenBase<Derived>& A) {

  if constexpr(static_cast<bool>(Eigen::internal::is_diagonal<Derived>::ret)) //workaround needed since Eigen::DiagonalWrapper does not has a eval function
  {
    using Scalar = typename Derived::Scalar;
    using namespace Eigen;
    constexpr int diag_size = EIGEN_SIZE_MIN_PREFER_DYNAMIC(Derived::RowsAtCompileTime, Derived::ColsAtCompileTime);
    constexpr int max_diag_size = EIGEN_SIZE_MIN_PREFER_FIXED(Derived::MaxRowsAtCompileTime, Derived::MaxColsAtCompileTime);

    return Eigen::DiagonalMatrix<Scalar,diag_size,max_diag_size>(A.derived().diagonal());}
  else
    return A.derived().eval();
}

/** \brief Does nothing if type is not an Eigen type but our manifolds type instead*/
auto eval(const Ikarus::Concepts::Manifold auto& A) {
    return A;
}

