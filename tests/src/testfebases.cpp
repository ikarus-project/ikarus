// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <ikarus/finiteelements/febases.hh>

using Dune::TestSuite;

#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>

template <typename GridView, typename PreBasis>
auto getBasis(const GridView& gridView, const PreBasis& preBasis) {
  const auto basis      = Ikarus::makeBasis(gridView, preBasis);
  const auto& localView = basis.flat().localView();
  auto element          = elements(gridView).begin();
  return Ikarus::FEBases<decltype(basis)>(basis, *element);
}

template <typename TestSuiteType, typename... Bases>
void checkBases(TestSuiteType& t, const std::vector<int>& expectedNumberOfChildren,
                const std::vector<int>& expectedSize, const std::string& basisName, const Bases&... bases) {
  assert(sizeof...(bases) == expectedNumberOfChildren.size() && "Input size mismatch for number of children.");
  assert(sizeof...(bases) == expectedSize.size() && "Input size mismatch for local view size.");
  int i = 0;
  (
      [&] {
        checkScalars(t, bases.numberOfChildren(), expectedNumberOfChildren[i],
                     " Number of children is incorrect for i = " + std::to_string(i) + basisName);
        checkScalars(t, static_cast<int>(bases.size()), expectedSize[i],
                     " Size is incorrect for i = " + std::to_string(i) + basisName);
        std::vector<typename std::remove_cvref_t<decltype(bases.localView())>::MultiIndex> dofs;
        bases.globalFlatIndices(dofs);
        checkScalars(t, dofs.size(), bases.size(),
                     " Size of dofs and args are not equal for i = " + std::to_string(i) + basisName);
        ++i;
      }(),
      ...);
}

auto FEBasesTest() {
  TestSuite t("FEBasesTest");
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox       = {4, 2};
  std::array<int, 2> elementsPerDirection = {1, 1};
  const auto grid                         = std::make_shared<Grid>(bbox, elementsPerDirection);

  const auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;

  const auto firstOrderLagrangePreBasis  = lagrange<1>();
  const auto secondOrderLagrangePreBasis = lagrange<2>();

  const auto firstOrderLagrangePowerPreBasis  = power<2>(firstOrderLagrangePreBasis);
  const auto secondOrderLagrangePowerPreBasis = power<5>(secondOrderLagrangePreBasis);

  const auto scalarScalarCompositePreBasis = composite(firstOrderLagrangePreBasis, secondOrderLagrangePreBasis);
  const auto scalarPowerCompositePreBasis  = composite(firstOrderLagrangePreBasis, secondOrderLagrangePowerPreBasis);

  const auto powerScalarCompositePreBasis = composite(firstOrderLagrangePowerPreBasis, secondOrderLagrangePreBasis);
  const auto powerPowerCompositePreBasis = composite(firstOrderLagrangePowerPreBasis, secondOrderLagrangePowerPreBasis);

  const auto scalarCompositePreBasis = composite(firstOrderLagrangePreBasis);
  const auto powerCompositePreBasis  = composite(firstOrderLagrangePowerPreBasis);

  /// Types of scalar FE bases
  const auto scalar1 = getBasis(gridView, firstOrderLagrangePreBasis);
  const auto scalar2 = getBasis(gridView, secondOrderLagrangePreBasis);

  /// Types of power FE bases
  const auto power1 = getBasis(gridView, firstOrderLagrangePowerPreBasis);
  const auto power2 = getBasis(gridView, secondOrderLagrangePowerPreBasis);

  /// Types of composite FE bases
  const auto composite1 = getBasis(gridView, scalarScalarCompositePreBasis);
  const auto composite2 = getBasis(gridView, scalarPowerCompositePreBasis);
  const auto composite3 = getBasis(gridView, powerScalarCompositePreBasis);
  const auto composite4 = getBasis(gridView, powerPowerCompositePreBasis);
  const auto composite5 = getBasis(gridView, scalarCompositePreBasis);
  const auto composite6 = getBasis(gridView, powerCompositePreBasis);

  checkBases(t, {0, 0}, {4, 9}, " Scalar Basis", scalar1, scalar2);
  checkBases(t, {2, 5}, {8, 45}, " Power Basis", power1, power2);
  checkBases(t, {2, 2, 2, 2, 1, 1}, {13, 49, 17, 53, 4, 8}, " Composite Basis", composite1, composite2, composite3,
             composite4, composite5, composite6);

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(FEBasesTest());
  return t.exit();
}
