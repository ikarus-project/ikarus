// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file masslumping.hh
 * \brief Helper for dune-functions
 */

#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/makeenum.hh>

namespace Ikarus::Dynamics {

namespace Impl {
  inline auto lumpingSparse = []<typename ScalarType>(Eigen::SparseMatrix<ScalarType>& mat) -> void {
    for (auto i : Dune::range(mat.rows())) {
      auto sum           = mat.row(i).sum();
      mat.coeffRef(i, i) = sum;
    }
    // Deletes all entries expect for main diagonal
    mat.prune([](int i, int j, auto) { return i == j; });
  };

  inline auto lumpingDense = []<typename ScalarType>(Eigen::MatrixX<double>& mat) -> void {
    for (auto i : Dune::range(mat.rows())) {
      auto sum  = mat.row(i).sum();
      mat(i, i) = sum;
    }
    mat = mat.diagonal().asDiagonal();
  };
} // namespace Impl

template <typename ScalarType = double>
inline auto rowSumLumpingSparse() {
  return [](const auto&, const auto&, auto, auto, Eigen::SparseMatrix<ScalarType>& mat) {
    return Impl::lumpingSparse.operator()(mat);
  };
}

template <typename ScalarType = double>
inline auto rowSumLumpingDense() {
  return [](const auto&, const auto&, auto, auto, Eigen::MatrixX<ScalarType>& mat) {
    return Impl::lumpingDense.operator()(mat);
  };
}

template <Concepts::FlatAssembler AS>
auto sparseLumpedAssembler(const std::shared_ptr<AS>& assembler) {
  auto assMrowLumped = Ikarus::makeAssemblerManipulator(*assembler);
  assMrowLumped->bind(Ikarus::Dynamics::rowSumLumpingSparse());
  return assMrowLumped;
}

} // namespace Ikarus::Dynamics
