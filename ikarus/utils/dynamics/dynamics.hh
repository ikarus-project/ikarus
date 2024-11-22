// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file dynamics.hh
 * \brief Helper for
 */

#pragma once
#include <dune/vtk/pvdwriter.hh>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/io/vtkwriter.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/makeenum.hh>

namespace Ikarus::Dynamics {

namespace Impl { // namespace Impl

  template <typename ScalarType = double>
  inline auto lumpingSparse() {
    return [](const auto&, const auto&, auto, auto, Eigen::SparseMatrix<ScalarType>& mat) {
      for (auto i : Dune::range(mat.rows())) {
        auto sum           = mat.row(i).sum();
        mat.coeffRef(i, i) = sum;
      }
      // Deletes all entries expect for main diagonal (entriese are really getting deleted)
      mat.prune([](int i, int j, auto) { return i == j; });
    };
  }

  template <typename ScalarType = double>
  inline auto lumpingDense() {
    return [](const auto&, const auto&, auto, auto, Eigen::MatrixX<ScalarType>& mat) {
      for (auto i : Dune::range(mat.rows())) {
        auto sum  = mat.row(i).sum();
        mat(i, i) = sum;
      }
      mat = mat.diagonal().asDiagonal();
    };
  };
} // namespace Impl

template <Concepts::FlatAssembler AS>
auto makeLumpedFlatAssembler(const std::shared_ptr<AS>& assembler) {
  constexpr auto isSparse = Concepts::SparseEigenMatrix<typename AS::MatrixType>;
  using ScalarType        = typename AS::MatrixType::Scalar;
  auto lumpedAssembler    = makeAssemblerManipulator(*assembler);
  if constexpr (isSparse)
    lumpedAssembler->bind(Impl::lumpingSparse<ScalarType>());
  else
    lumpedAssembler->bind(Impl::lumpingDense<ScalarType>());
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
  writer.write(filename);
}

template <Concepts::EigenValueSolver Eigensolver, Concepts::FlatAssembler Assembler>
void writeEigenformsToPVD(const Eigensolver& solver, std::shared_ptr<Assembler> assembler, const std::string& filename,
                          std::optional<Eigen::Index> nev_ = std::nullopt) {
  auto nev          = nev_.value_or(solver.nev());
  auto eigenvectors = solver.eigenvectors();
  auto basis        = assembler->basis();

  auto writer    = std::make_shared<decltype(Ikarus::Vtk::Writer(assembler))>(Ikarus::Vtk::Writer(assembler));
  auto pvdWriter = Dune::Vtk::PvdWriter(writer);

  Eigen::VectorXd evG(assembler->size());
  writer->addInterpolation(evG, basis, "EF");

  for (auto i : Dune::range(nev)) {
    evG = assembler->createFullVector(eigenvectors.col(i));
    pvdWriter.writeTimestep(i, filename, "data", false);
  }
  pvdWriter.write(filename);
}

} // namespace Ikarus::Dynamics
