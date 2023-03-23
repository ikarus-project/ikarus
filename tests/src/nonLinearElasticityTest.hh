// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

//
#include <config.h>

#include "common.hh"
#include "testHelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include "spdlog/spdlog.h"

#include <Eigen/Core>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlRoutines/loadControl.hh>
#include <ikarus/finiteElements/mechanics/nonLinearElasticityFE.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/solver/nonLinearSolver/trustRegion.hh>
#include <ikarus/utils/algorithms.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>

using Dune::TestSuite;

template <typename Grid, typename Material>
auto NonLinearElasticityLoadControlNRandTR(const Material& mat) {
  TestSuite t("NonLinearElasticityLoadControlNRandTR" + Dune::className(Grid{}));
  auto grid     = createGrid<Grid>();
  auto gridView = grid->leafGridView();

  using namespace Ikarus;

  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));

  auto g          = basis.localView();
  auto volumeLoad = []([[maybe_unused]] auto& globalCoord, auto& lamb) {
    Eigen::Vector2d fext;
    fext.setZero();
    fext[1] = 2 * lamb;
    fext[0] = lamb;
    return fext;
  };

  auto reducedMat = plainStress(mat, 1e-8);

  std::vector<Ikarus::NonLinearElasticityFE<decltype(basis), decltype(reducedMat)>> fes;

  for (auto& element : elements(gridView))
    fes.emplace_back(basis, element, reducedMat, nullptr, nullptr, volumeLoad);
  auto basisP = std::make_shared<const decltype(basis)>(basis);
  Ikarus::DirichletValues dirichletValues(basisP);
  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[1]) < 1e-8) dirichletFlags[localView.index(localIndex)] = true;
  });

  auto sparseAssembler = SparseFlatAssembler(fes, dirichletValues);

  Eigen::VectorXd d;
  d.setZero(basis.size());
  double lambda = 0.0;

  auto req = FErequirements().addAffordance(Ikarus::AffordanceCollections::elastoStatics);

  auto residualFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return sparseAssembler.getVector(req);
  };

  auto KFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return sparseAssembler.getMatrix(req);
  };

  auto energyFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return sparseAssembler.getScalar(req);
  };

  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(energyFunction, residualFunction, KFunction),
                                            parameter(d, lambda));
  //  t.check(checkGradient(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = true})) << "checkGradient Failed";
  //  t.check(checkHessian(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = true})) << "checkHessian Failed";

  const double gradTol = 1e-8;

  auto tr = Ikarus::makeTrustRegion(nonLinOp);
  tr->setup({.verbosity = 1,
             .maxiter   = 1000,
             .grad_tol  = gradTol,
             .corr_tol  = 1e-16,  // everything should converge to the gradient tolerance
             .useRand   = false,
             .rho_reg   = 1e8,
             .Delta0    = 1});

  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<decltype(basis)>>(basis, d, 2);
  vtkWriter->setFileNamePrefix("Test2Dsolid");
  vtkWriter->setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 2);

  auto lc = Ikarus::LoadControl(tr, 1, {0, 50});
  lc.subscribeAll(vtkWriter);
  const auto controlInfo = lc.run();
  nonLinOp.template update<0>();
  const auto maxDisp = std::ranges::max(d);
  const double energyExpected
      = (std::is_same_v<Grid, Grids::Yasp>)
            ? -2.9593431593780032962
            : ((std::is_same_v<Grid, Grids::Alu>) ? -2.9530594665063669702
                                                  : /* std::is_same_v<Grid, Grids::Iga> */ -1.4533281398929942529);
  const double maxDispExpected
      = (std::is_same_v<Grid, Grids::Yasp>)
            ? 0.11291304159624337977
            : ((std::is_same_v<Grid, Grids::Alu>) ? 0.1123397197762363714
                                                  : /* std::is_same_v<Grid, Grids::Iga> */ 0.061647849558021668159);
  std::cout << std::setprecision(20) << nonLinOp.value() << std::endl;
  std::cout << "Maxdisp: " << maxDisp << std::endl;
  if constexpr (std::is_same_v<Material, Ikarus::StVenantKirchhoff<>>) {
    t.check(Dune::FloatCmp::eq(energyExpected, nonLinOp.value()), "energyExpected == nonLinOp.value()")
        << "energyExpected: " << energyExpected << "\nnonLinOp.value(): " << nonLinOp.value();

    t.check(std::abs(maxDispExpected - maxDisp) < 1e-12, "maxDispExpected-maxDisp")
        << "maxDispExpected: \n"
        << maxDispExpected << "\nmaxDisp: \n"
        << maxDisp;
  } else {  // using a Neohooke material yields a lower energy and larger displacements
    t.check(Dune::FloatCmp::gt(energyExpected, nonLinOp.value()), "energyExpected > nonLinOp.value()")
        << "energyExpected: " << energyExpected << "\nnonLinOp.value(): " << nonLinOp.value();

    t.check(maxDispExpected < maxDisp, "maxDispExpected<maxDisp") << "maxDispExpected: \n"
                                                                  << maxDispExpected << "\nmaxDisp: \n"
                                                                  << maxDisp;
  }

  nonLinOp.template update<1>();
  t.check(controlInfo.success, "Successful result");
  t.check(gradTol >= nonLinOp.derivative().norm(), "Gradient Tolerance should be larger than actual tolerance");
  return t;
}
