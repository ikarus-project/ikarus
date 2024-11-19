// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file dynamics.hh
 * \brief Helper for
 */

#pragma once
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/dynamics/generaleigensolver.hh>
#include <ikarus/utils/makeenum.hh>

namespace Ikarus::Dynamics {

namespace Impl {
  inline auto lumpingSparse = []<typename ScalarType>(Eigen::SparseMatrix<ScalarType>& mat) -> void {
    for (auto i : Dune::range(mat.rows())) {
      auto sum           = mat.row(i).sum();
      mat.coeffRef(i, i) = sum;
    }
    // Deletes all entries expect for main diagonal (entries are really getting deleted)
    mat.prune([](int i, int j, auto) { return i == j; });
  };

  inline auto lumpingDense = []<typename ScalarType>(Eigen::MatrixX<ScalarType>& mat) -> void {
    for (auto i : Dune::range(mat.rows())) {
      auto sum  = mat.row(i).sum();
      mat(i, i) = sum;
    }
    mat = mat.diagonal().asDiagonal();
  };

  template <typename ScalarType = double>
  inline auto rowSumLumpingSparse() {
    return
        [](const auto&, const auto&, auto, auto, Eigen::SparseMatrix<ScalarType>& mat) { return lumpingSparse(mat); };
  }

  template <typename ScalarType = double>
  inline auto rowSumLumpingDense() {
    return [](const auto&, const auto&, auto, auto, Eigen::MatrixX<ScalarType>& mat) { return lumpingDense(mat); };
  }
} // namespace Impl

template <Concepts::FlatAssembler AS>
auto makeLumpedFlatAssembler(const std::shared_ptr<AS>& assembler) {
  constexpr auto isSparse = Concepts::SparseEigenMatrix<typename AS::MatrixType>;
  using ScalarType        = typename AS::MatrixType::Scalar;
  auto lumpedAssembler    = Ikarus::makeAssemblerManipulator(*assembler);
  if constexpr (isSparse)
    lumpedAssembler->bind(Ikarus::Dynamics::Impl::rowSumLumpingSparse<ScalarType>());
  else
    lumpedAssembler->bind(Ikarus::Dynamics::Impl::rowSumLumpingDense<ScalarType>());
  return lumpedAssembler;
}
} // namespace Ikarus::Dynamics
