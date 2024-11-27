// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file modalanalysishelper.hh
 * \brief Helper for modal analysis related things
 */

#pragma once
#include <dune/vtk/pvdwriter.hh>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/io/vtkwriter.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/makeenum.hh>
#include <ikarus/utils/modalanalysis/lumpingschemes.hh>

namespace Ikarus::Dynamics {

template <Concepts::FlatAssembler AS, typename LumpingScheme = LumpingSchemes::RowSumLumping>
auto makeLumpedFlatAssembler(const std::shared_ptr<AS>& assembler) {
  constexpr auto isSparse = Concepts::SparseEigenMatrix<typename AS::MatrixType>;
  using ScalarType        = typename AS::MatrixType::Scalar;
  auto lumpedAssembler    = makeAssemblerManipulator(*assembler);
  lumpedAssembler->bind(LumpingScheme{});
  return lumpedAssembler;
}

template <Concepts::EigenValueSolver Eigensolver, Concepts::FlatAssembler Assembler>
void writeEigenmodesToVTK(const Eigensolver& solver, std::shared_ptr<Assembler> assembler, const std::string& filename,
                          std::optional<Eigen::Index> nev_ = std::nullopt) {
  auto nev          = nev_.value_or(solver.nev());
  auto eigenvectors = solver.eigenvectors();
  auto basis        = assembler->basis();

  auto writer = Ikarus::Vtk::Writer(assembler);
  for (auto i : Dune::range(nev)) {
    auto evG = assembler->createFullVector(eigenvectors.col(i));
    writer.addInterpolation(std::move(evG), basis, "Eigenmode " + std::to_string(i));
  }
  writer.write(filename);
}

template <typename Derived, Concepts::FlatAssembler Assembler>
void writeEigenmodesAsTimeSeries(const Eigen::MatrixBase<Derived>& eigenvectors, std::shared_ptr<Assembler> assembler,
                                 const std::string& filename) {
  auto nev   = eigenvectors.cols();
  auto basis = assembler->basis();

  auto writer    = std::make_shared<decltype(Ikarus::Vtk::Writer(assembler))>(Ikarus::Vtk::Writer(assembler));
  auto pvdWriter = Dune::Vtk::PvdWriter(writer);

  Eigen::VectorXd evG(assembler->size());
  writer->addInterpolation(evG, basis, "Eigenmode");

  for (auto i : Dune::range(nev)) {
    evG = assembler->createFullVector(eigenvectors.col(i));
    pvdWriter.writeTimestep(i + 1, filename, "data", false);
  }
  pvdWriter.write(filename);
}

template <Concepts::EigenValueSolver Eigensolver, Concepts::FlatAssembler Assembler>
void writeEigenmodesAsTimeSeries(const Eigensolver& solver, std::shared_ptr<Assembler> assembler,
                                 const std::string& filename, std::optional<Eigen::Index> nev_ = std::nullopt) {
  auto nev          = nev_.value_or(solver.nev());
  auto eigenvectors = solver.eigenvectors();
  writeEigenmodesAsTimeSeries(eigenvectors, assembler, filename);
}

} // namespace Ikarus::Dynamics
