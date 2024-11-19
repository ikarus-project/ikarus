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
#include <ikarus/io/vtkwriter.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/dynamics/generaleigensolver.hh>
#include <ikarus/utils/makeenum.hh>

namespace Ikarus::Dynamics {

namespace Impl {
  inline auto lumpingSparseImpl = []<typename ScalarType>(Eigen::SparseMatrix<ScalarType>& mat) -> void {
    for (auto i : Dune::range(mat.rows())) {
      auto sum           = mat.row(i).sum();
      mat.coeffRef(i, i) = sum;
    }
    // Deletes all entries expect for main diagonal (entries are really getting deleted)
    mat.prune([](int i, int j, auto) { return i == j; });
  };

  inline auto lumpingDenseImpl = []<typename ScalarType>(Eigen::MatrixX<ScalarType>& mat) -> void {
    for (auto i : Dune::range(mat.rows())) {
      auto sum  = mat.row(i).sum();
      mat(i, i) = sum;
    }
    mat = mat.diagonal().asDiagonal();
  };

} // namespace Impl

template <typename ScalarType = double>
inline auto lumpingSparse() {
  return
      [](const auto&, const auto&, auto, auto, Eigen::SparseMatrix<ScalarType>& mat) { return Impl::lumpingSparseImpl(mat); };
}

template <typename ScalarType = double>
inline auto lumpingDense() {
  return [](const auto&, const auto&, auto, auto, Eigen::MatrixX<ScalarType>& mat) { return Impl::lumpingDenseImpl(mat); };
}

template <Concepts::FlatAssembler AS>
auto makeLumpedFlatAssembler(const std::shared_ptr<AS>& assembler) {
  constexpr auto isSparse = Concepts::SparseEigenMatrix<typename AS::MatrixType>;
  using ScalarType        = typename AS::MatrixType::Scalar;
  auto lumpedAssembler    = makeAssemblerManipulator(*assembler);
  if constexpr (isSparse)
    lumpedAssembler->bind(lumpingSparse<ScalarType>());
  else
    lumpedAssembler->bind(lumpingDense<ScalarType>());
  return lumpedAssembler;
}

template <Concepts::EigenValueSolver Eigensolver, Concepts::FlatAssembler Assembler>
void writeEigenformsToVTK(const Eigensolver& solver, std::shared_ptr<Assembler> assembler, const std::string& filename,
                          std::optional<Eigen::Index> nev_ = std::nullopt) {
  auto nev          = nev_.value_or(solver.nev());
  auto eigenvectors = solver.eigenvectors();
  auto basis        = assembler->basis();

  auto writer = Ikarus::Vtk::Writer(assembler);
  for (auto i : Dune::range(nev)) {
    auto evG = assembler->createFullVector(eigenvectors.col(i));
    writer.addInterpolation(std::move(evG), basis, "EF " + std::to_string(i));
  }
  writer.write("eigenformen");
}

} // namespace Ikarus::Dynamics
