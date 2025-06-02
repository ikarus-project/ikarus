// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcantileverbeam.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/displacementpressure.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mixin.hh>

using namespace Ikarus;
using Dune::TestSuite;

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t;
  auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = 100.0, .nu = 0.3});
  Materials::StVenantKirchhoff matSVK(matParameter);
  Materials::NeoHooke matNH(matParameter);
  auto matBK = Materials::makeBlatzKo(40.0);

  // auto reducedMats = Dune::makeTupleVector(planeStrain(matSVK), planeStrain(matNH), planeStrain(matBK));

  using Grid     = Dune::YaspGrid<2>;
  const double L = 50;

  Dune::FieldVector<double, 2> bbox       = {L, L};
  std::array<int, 2> elementsPerDirection = {1, 1};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);
  auto gridView                           = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(
      gridView, composite(power<2>(lagrange<1>(), FlatInterleaved{}), lagrange<0>(), BlockedLexicographic{}));

  auto sk      = skills(displacementPressure(planeStrain(matNH)));
  using FEType = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;

  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  auto nDOF = basis.flat().size();

  auto& fe = fes.front();
  Eigen::VectorXd d;
  d.setZero(nDOF);
  double lambda = 0.0;
  d << 2, 4, 3.25, -1.2, 0.003, 6, 3, 2.864, 1e-4;
  auto req = typename FEType::Requirement(d, lambda);

  Eigen::MatrixXd k;
  k.setZero(9, 9);

  calculateMatrix(fe, req, Ikarus::MatrixAffordance::stiffness, k);

  std::cout << k << std::endl;

  // Ikarus::DirichletValues dirichletValues(basis.flat());

  // auto sparseFlatAssembler = makeSparseFlatAssembler(fes, dirichletValues);
  // sparseFlatAssembler->bind(req, Ikarus::AffordanceCollections::elastoStatics, DBCOption::Full);

  // auto& K = sparseFlatAssembler->matrix();

  // std::cout << K.toDense() << std::endl;

  return t.exit();
}
