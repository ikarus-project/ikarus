// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#include <config.h>

#include "tests/src/testhelpers.hh"

#include <vector>

#include <dune/common/float_cmp.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/modalanalysis/modalanalysis.hh>

using Dune::TestSuite;

static auto modalAnalysisTest() {
  TestSuite t("ModalAnalysis test");
  using Grid = Dune::YaspGrid<2>;

  const double Lx                         = 4.0;
  const double Ly                         = 4.0;
  Dune::FieldVector<double, 2> bbox       = {Lx, Ly};
  std::array<int, 2> elementsPerDirection = {5, 5};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);

  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>()));

  auto linMat  = Ikarus::LinearElasticity(Ikarus::toLamesFirstParameterAndShearModulus({.emodul = 100, .nu = 0.2}));
  auto skills_ = Ikarus::skills(Ikarus::linearElastic(Ikarus::planeStress(linMat)));

  using FEType = decltype(Ikarus::makeFE(basis, skills_));
  std::vector<FEType> fes;
  for (auto&& element : elements(gridView)) {
    fes.emplace_back(Ikarus::makeFE(basis, skills_));
    fes.back().bind(element);
  }

  auto dirichletValues = Ikarus::DirichletValues(basis.flat());
  dirichletValues.fixBoundaryDOFs([](auto& dirichFlags, auto&& indexGlobal) { dirichFlags[indexGlobal] = true; });

  auto mA = Ikarus::Dynamics::ModalAnalysis(std::move(fes), dirichletValues);
  mA.compute();

  auto frequencies = mA.angularFrequencies();

  mA.bindLumpingScheme<Ikarus::Dynamics::LumpingSchemes::RowSumLumping>();
  mA.compute();

  auto frequenciesLumped = mA.angularFrequencies();
  t.check(frequencies.sum() > frequenciesLumped.sum()) << testLocation();

  mA.unBindLumpingScheme();
  mA.compute();
  auto frequencies2 = mA.angularFrequencies();
  t.check(isApproxSame(frequencies, frequencies2, 1e-14)) << testLocation();

  mA.writeEigenModes("eigenforms", 20);

  auto req = FEType::Requirement();
  typename FEType::Requirement::SolutionVectorType d;
  d.setZero(basis.flat().size());
  req.insertGlobalSolution(d);

  auto massAssembler = Ikarus::makeSparseFlatAssembler(fes, dirichletValues);
  massAssembler->bind(req, Ikarus::AffordanceCollections::dynamics, Ikarus::DBCOption::Reduced);
  auto lumpedAssembler = Ikarus::Dynamics::makeLumpedFlatAssembler(massAssembler);

  auto M       = massAssembler->matrix();
  auto MLumped = lumpedAssembler->matrix();

  t.check(MLumped.nonZeros() == massAssembler->reducedSize()) << testLocation();
  t.check(MLumped.coeff(0, 0) == M.row(0).sum()) << testLocation();

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(modalAnalysisTest());
  return t.exit();
}
