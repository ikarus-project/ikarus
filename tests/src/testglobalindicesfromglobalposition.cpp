// SPDX-FileCopyrightText: 2021-2026 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/uggrid.hh>

#include <ikarus/utils/functionhelper.hh>

using Dune::TestSuite;

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/functionhelper.hh>
#include <ikarus/utils/init.hh>

auto testGlobalIndicesFromGlobalPosition() {
  TestSuite t("Test global indices from global position");
  using namespace Ikarus;
  constexpr int gridDim = 2;
  Dune::GridFactory<Dune::UGGrid<gridDim>> gridFactory;

  Dune::FieldVector<double, gridDim> bottomLeft  = {0.0, 0.0};
  Dune::FieldVector<double, gridDim> bottomRight = {1.0, 0.0};
  Dune::FieldVector<double, gridDim> topLeft     = {3.3e-12, 0.99};
  Dune::FieldVector<double, gridDim> topRight    = {1.0, 1.0};

  gridFactory.insertVertex(bottomLeft);  // 0
  gridFactory.insertVertex(bottomRight); // 1
  gridFactory.insertVertex(topLeft);     // 2
  gridFactory.insertVertex(topRight);    // 3

  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {0, 1, 2, 3});
  auto grid = gridFactory.createGrid();

  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis0 = Ikarus::makeBasis(gridView, lagrange<1>());
  auto basis1 = Ikarus::makeBasis(gridView, power<5>(lagrange<1>(), FlatInterleaved()));
  auto basis2 = Ikarus::makeBasis(
      gridView, composite(power<2>(lagrange<1>(), FlatInterleaved{}), lagrange<0>(), BlockedLexicographic{}));

  const auto& flatBasis0 = basis0.flat();
  const auto& flatBasis1 = basis1.flat();
  const auto& flatBasis2 = basis2.flat();

  const auto bottomLeftIndices0 = Ikarus::utils::globalIndexFromGlobalPosition(flatBasis0, bottomLeft);
  const auto bottomLeftIndices1 = Ikarus::utils::globalIndexFromGlobalPosition(flatBasis1, bottomLeft);
  const auto bottomLeftIndices2 = Ikarus::utils::globalIndexFromGlobalPosition(flatBasis2, bottomLeft);

  const auto topLeftIndices0 = Ikarus::utils::globalIndexFromGlobalPosition(flatBasis0, topLeft);
  const auto topLeftIndices1 = Ikarus::utils::globalIndexFromGlobalPosition(flatBasis1, topLeft);
  const auto topLeftIndices2 = Ikarus::utils::globalIndexFromGlobalPosition(flatBasis2, topLeft);

  const auto bottomRightIndices0 = Ikarus::utils::globalIndexFromGlobalPosition(flatBasis0, bottomRight);
  const auto bottomRightIndices1 = Ikarus::utils::globalIndexFromGlobalPosition(flatBasis1, bottomRight);
  const auto bottomRightIndices2 = Ikarus::utils::globalIndexFromGlobalPosition(flatBasis2, bottomRight);

  const auto topRightIndices0 = Ikarus::utils::globalIndexFromGlobalPosition(flatBasis0, topRight);
  const auto topRightIndices1 = Ikarus::utils::globalIndexFromGlobalPosition(flatBasis1, topRight);
  const auto topRightIndices2 = Ikarus::utils::globalIndexFromGlobalPosition(flatBasis2, topRight);

  checkScalars(t, bottomLeftIndices0[0], std::size_t{0}, " Incorrect bottom left scalar index");
  checkScalars(t, bottomRightIndices0[0], std::size_t{1}, " Incorrect bottom right scalar index");
  checkScalars(t, topLeftIndices0[0], std::size_t{2}, " Incorrect top left scalar index");
  checkScalars(t, topRightIndices0[0], std::size_t{3}, " Incorrect top right scalar index");

  std::array<std::size_t, 5> expBL1 = {0, 1, 2, 3, 4};
  std::array<std::size_t, 5> expBR1 = {5, 6, 7, 8, 9};
  std::array<std::size_t, 5> expTL1 = {10, 11, 12, 13, 14};
  std::array<std::size_t, 5> expTR1 = {15, 16, 17, 18, 19};

  for (const auto i : Dune::range(5)) {
    checkScalars(t, bottomLeftIndices1[i][0], expBL1[i], " Incorrect bottom left power index");
    checkScalars(t, bottomRightIndices1[i][0], expBR1[i], " Incorrect bottom right power index");
    checkScalars(t, topLeftIndices1[i][0], expTL1[i], " Incorrect top left power index");
    checkScalars(t, topRightIndices1[i][0], expTR1[i], " Incorrect top right power index");
  }

  std::array<std::size_t, 2> expBL2 = {0, 1};
  std::array<std::size_t, 2> expBR2 = {2, 3};
  std::array<std::size_t, 2> expTL2 = {4, 5};
  std::array<std::size_t, 2> expTR2 = {6, 7};

  for (const auto i : Dune::range(2)) {
    checkScalars(t, bottomLeftIndices2[i][0], expBL2[i], " Incorrect bottom left composite index");
    checkScalars(t, bottomRightIndices2[i][0], expBR2[i], " Incorrect bottom right composite index");
    checkScalars(t, topLeftIndices2[i][0], expTL2[i], " Incorrect top left composite index");
    checkScalars(t, topRightIndices2[i][0], expTR2[i], " Incorrect top right composite index");
  }

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;
  t.subTest(testGlobalIndicesFromGlobalPosition());
  return t.exit();
}
