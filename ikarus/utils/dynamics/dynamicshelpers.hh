// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file dynamicshelpers.hh
 * \brief Helper for Dynamic related things
 */

#pragma once
#include <dune/vtk/pvdwriter.hh>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/io/vtkwriter.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/dynamics/lumpingschemes.hh>
#include <ikarus/utils/makeenum.hh>

namespace Ikarus::Dynamics {

template <Concepts::FlatAssembler AS>
auto makeLumpedFlatAssembler(const std::shared_ptr<AS>& assembler) {
  constexpr auto isSparse = Concepts::SparseEigenMatrix<typename AS::MatrixType>;
  using ScalarType        = typename AS::MatrixType::Scalar;
  auto lumpedAssembler    = makeAssemblerManipulator(*assembler);
  lumpedAssembler->bind(LumpingSchemes::RowSumLumping{});
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
    writer.addInterpolation(std::move(evG), basis, "EF " + std::to_string(i));
  }
  writer.write(filename);
}

template <Concepts::EigenValueSolver Eigensolver, Concepts::FlatAssembler Assembler>
void writeEigenmodesToPVD(const Eigensolver& solver, std::shared_ptr<Assembler> assembler, const std::string& filename,
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
