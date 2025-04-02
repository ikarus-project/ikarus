// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcommon.hh"
#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <Eigen/Core>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/nonlinearelastic.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverfactory.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/functionhelper.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/listener/controllogger.hh>
#include <ikarus/utils/listener/controlvtkwriter.hh>
#include <ikarus/utils/listener/nonlinearsolverlogger.hh>

using namespace Ikarus;
using Dune::TestSuite;

template <typename ES, typename MAT>
auto cantileverBeamResults() {
  static_assert(MAT::isReduced,
                "cantileverBeamResults are available only for a reduced material (planeStress or planeStrain).");
  using namespace Materials;
  using MATU = typename MAT::Underlying;
  if constexpr (std::same_as<ES, EAS::GreenLagrangeStrain>) {
    if constexpr (std::is_same_v<MATU, StVenantKirchhoff>) {
      constexpr int expectedIterations = 80;
      constexpr double expectedMaxDisp = 4.459851990257645;
      return std::make_pair(expectedIterations, expectedMaxDisp);
    } else if constexpr (std::is_same_v<MATU, NeoHooke>) {
      constexpr int expectedIterations = 80;
      constexpr double expectedMaxDisp = 4.479930218997457;
      return std::make_pair(expectedIterations, expectedMaxDisp);
    } else if constexpr (std::is_same_v<MATU, Hyperelastic<Deviatoric<BlatzKo>, Volumetric<VF0>>>) {
      constexpr int expectedIterations = 80;
      constexpr double expectedMaxDisp = 4.6877262928164365;
      return std::make_pair(expectedIterations, expectedMaxDisp);
    } else
      static_assert(Dune::AlwaysFalse<MATU>::value,
                    "Expected results are not available for the given material and enhancedStrain type.");
  } else
    static_assert(Dune::AlwaysFalse<MATU>::value,
                  "Expected results are not available for the provided enhanced strain type.");
}

template <typename ES, typename MAT>
auto cantileverBeamTest(const MAT& reducedMat) {
  static_assert(MAT::isReduced,
                "cantileverBeamTest is only valid for a reduced material (planeStress or planeStrain).");
  TestSuite t("Cantilever Beam for Nonlinear EAS element with enhanced strain type as " + ES::name() +
              " and material type as " + Dune::className<MAT>());
  constexpr int gridDim = 2;
  using Grid            = Dune::YaspGrid<gridDim>;
  const double L        = 10;
  const double h        = 2;

  Dune::FieldVector<double, gridDim> bbox       = {L, h};
  std::array<int, gridDim> elementsPerDirection = {10, 1};
  auto grid                                     = std::make_shared<Grid>(bbox, elementsPerDirection);
  auto gridView                                 = grid->leafGridView();
  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));

  auto sk      = skills(nonLinearElastic(reducedMat), eas<ES>(4));
  using FEType = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;

  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  Ikarus::DirichletValues dirichletValues(basis.flat());

  // fix left edge (x=0)
  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[0]) < 1e-8)
      dirichletFlags[localView.index(localIndex)] = true;
  });

  auto sparseFlatAssembler = makeSparseFlatAssembler(fes, dirichletValues);
  auto sparseAssemblerAM   = makeAssemblerManipulator(*sparseFlatAssembler);

  Eigen::VectorXd d;
  d.setZero(basis.flat().size());
  double lambda = 0.0;

  auto req = typename FEType::Requirement(d, lambda);

  sparseAssemblerAM->bind(req, Ikarus::AffordanceCollections::elastoStatics, DBCOption::Full);

  // Apply constant point load at the top right corner
  Dune::FieldVector<double, gridDim> topRightPos{L, h};
  const auto globalIndices = utils::globalIndexFromGlobalPosition(basis.flat(), topRightPos);
  auto pointLoad           = [&](const auto&, const auto& par, auto, auto, Eigen::VectorXd& vec) -> void {
    auto loadFactor = par.parameter();
    vec[globalIndices[1]] -= -loadFactor * 1.0;
  };
  sparseAssemblerAM->bind(pointLoad);

  auto linSolver = LinearSolver(SolverTypeTag::sd_UmfPackLU);

  AffordanceCollection elastoStaticsNoScalar(VectorAffordance::forces, MatrixAffordance::stiffness);

  auto nonOp =
      DifferentiableFunctionFactory::op(sparseAssemblerAM, elastoStaticsNoScalar, sparseAssemblerAM->dBCOption());

  constexpr double tol = 1e-10;

  auto nrConfig =
      Ikarus::NewtonRaphsonConfig<decltype(linSolver)>{.parameters = {.tol = tol}, .linearSolver = linSolver};
  NonlinearSolverFactory nrFactory(nrConfig);
  auto nr = nrFactory.create(sparseAssemblerAM);

  // Only when creating the control routine via the Factory, the elements get registered for correction update
  // automatically.
  auto lc = ControlRoutineFactory::create(LoadControlConfig{20, 0.0, 1.0}, nr, sparseFlatAssembler);

  auto nonLinearSolverObserver = NonLinearSolverLogger();
  auto controlLogger           = ControlLogger();
  auto vtkWriter = ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>>(basis.flat(), d);
  vtkWriter.setFileNamePrefix("CantileverNonlinearEAS");
  vtkWriter.setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 2);

  nonLinearSolverObserver.subscribeTo(lc.nonLinearSolver());
  vtkWriter.subscribeTo(lc);
  controlLogger.subscribeTo(lc);

  const auto controlInfo = lc.run(req);
  d                      = req.globalSolution();
  lambda                 = req.parameter();

  double expectedLambda                            = 1.0;
  const auto [expectedIterations, expectedMaxDisp] = cantileverBeamResults<ES, MAT>();

  t.check(controlInfo.success);
  checkScalars(t, controlInfo.totalIterations, expectedIterations, " Number of iterations");
  const auto maxDisp = std::ranges::max(d.cwiseAbs());

  checkScalars(t, maxDisp, expectedMaxDisp, " Max. displacement", tol);
  checkScalars(t, lambda, expectedLambda, " Load factor", tol);

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t;
  auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = 100.0, .nu = 0.3});
  Materials::StVenantKirchhoff matSVK(matParameter);
  Materials::NeoHooke matNH(matParameter);
  auto matBK         = Materials::makeBlatzKo(40.0);
  auto reducedMatSVK = planeStrain(matSVK);
  auto reducedMatNH  = planeStrain(matNH);
  auto reducedMatBK  = planeStrain(matBK);

  const auto testNonlinearEASFunc = [&]<typename ES>() {
    t.subTest(cantileverBeamTest<ES>(reducedMatSVK));
    t.subTest(cantileverBeamTest<ES>(reducedMatNH));
    t.subTest(cantileverBeamTest<ES>(reducedMatBK));
  };

  testNonlinearEASFunc.operator()<EAS::GreenLagrangeStrain>();

  return t.exit();
}
