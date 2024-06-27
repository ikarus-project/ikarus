// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/vtk/datacollectors/discontinuousdatacollector.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/vtkreader.hh>

#include <Eigen/Core>

#include <ikarus/assembler/simpleassemblers.hh>
#include "ikarus/finiteelements/feresulttypes.hh"
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mechanics/loads/volume.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/io/resultfunction.hh>
#include <ikarus/io/vtkwriter.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/init.hh>
using Dune::TestSuite;


auto runTestCase() {
  TestSuite t("Test ResultFunction");
  std::string fileName = "ResultFunctionTest";

  using Grid     = Dune::UGGrid<2>;
  using GridView = Grid::LeafGridView;

  constexpr double Lx                                    = 4.0;
  constexpr double Ly                                    = 4.0;
  const Dune::FieldVector<double, 2> bbox                = {Lx, Ly};
  const std::array<unsigned int, 2> elementsPerDirection = {2, 2};

  auto grid     = Dune::StructuredGridFactory<Grid>::createCubeGrid({0, 0}, bbox, elementsPerDirection);
  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>()));

  Ikarus::DirichletValues dirichletValues(basis.flat());
  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[1]) < 1e-8)
      dirichletFlags[localView.index(localIndex)] = true;
  });

  auto vL      = []([[maybe_unused]] auto& globalCoord, auto& lamb) { return Eigen::Vector2d{0, -1}; };
  auto skills_ = Ikarus::skills(Ikarus::linearElastic({.emodul = 100, .nu = 0.2}), Ikarus::volumeLoad<2>(vL));

  using LinearElastic = decltype(Ikarus::makeFE(basis, skills_));
  std::vector<LinearElastic> fes;

  for (auto&& element : elements(gridView)) {
    fes.emplace_back(Ikarus::makeFE(basis, skills_));
    fes.back().bind(element);
  }

  /// Create a sparse assembler
  auto sparseAssembler = makeSparseFlatAssembler(fes, dirichletValues);

  Eigen::VectorXd D_Glob = Eigen::VectorXd::Zero(basis.flat().size());

  auto req        = LinearElastic::Requirement();
  auto lambdaLoad = 1.0;
  req.insertGlobalSolution(D_Glob).insertParameter(lambdaLoad);

  sparseAssembler->bind(req);
  sparseAssembler->bind(Ikarus::DBCOption::Full);

  auto nonLinOp = Ikarus::NonLinearOperatorFactory::op(
      sparseAssembler,
      Ikarus::AffordanceCollection(Ikarus::VectorAffordance::forces, Ikarus::MatrixAffordance::stiffness));

  const auto& K    = nonLinOp.derivative();
  const auto& Fext = nonLinOp.value();

  auto linSolver = Ikarus::LinearSolver(Ikarus::SolverTypeTag::sd_CholmodSupernodalLLT);
  linSolver.compute(K);
  linSolver.solve(D_Glob, -Fext);

  return t;
}

int main(const int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(runTestCase());

  return t.exit();
}