// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testCommon.hh"
#include "testHelpers.hh"

#include <dune/common/parametertreeparser.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/iga/nurbsbasis.hh>

#include "spdlog/spdlog.h"

#include <Eigen/Core>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlRoutines/loadControl.hh>
#include <ikarus/controlRoutines/pathFollowingTechnique.hh>
#include <ikarus/finiteElements/mechanics/fesettings.hh>
#include <ikarus/finiteElements/mechanics/kirchhoffloveshell.hh>
#include <ikarus/io/resultFunction.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/solver/nonLinearSolver/trustRegion.hh>
#include <ikarus/utils/algorithms.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <ikarus/utils/observer/nonLinearSolverLogger.hh>

using Dune::TestSuite;

auto checkFEByAutoDiff() {
  TestSuite t("Check calculateScalarImpl() and calculateVectorImpl() by Automatic Differentiation of Kirchhoff-Love shell");

  constexpr auto dimworld        = 3;
  const std::array<int, 2> order = {2, 2};

  const std::array<std::vector<double>, 2> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints
      = {
      {{.p = {0, 0.0, 0}, .w = 1}, {.p = {5, 0.0, 0}, .w = 1}, {.p = {10, 0.0, 0}, .w = 1}},
      {{.p = {0, 0.5, 0}, .w = 1}, {.p = {5, 0.5, 0}, .w = 4}, {.p = {10, 0.5, 0}, .w = 1}},
      {{.p = {0, 2, 0}, .w = 1}, {.p = {5, 2, 0}, .w = 1}, {.p = {10, 2, 0}, .w = 1}}};

  std::array<int, 2> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<2, dimworld>;

  Dune::IGA::NURBSPatchData<2, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
//  for (int i = 0; i < 2; ++i)
//    patchData = degreeElevate(patchData, i, 1);
  auto grid = std::make_shared<Grid>(patchData);

  grid->globalRefine(2);
  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis       = Ikarus::makeBasis(gridView, power<3>(nurbs()));
  auto element     = gridView.template begin<0>();
  auto nDOF        = basis.flat().size();
  const double tol = 1e-10;

  auto volumeLoad = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    VectorType fExt;
    fExt.setZero();
    fExt[1] = 2 * lamb;
    return fExt;
  };

  auto neumannBoundaryLoad = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    VectorType fExt;
    fExt.setZero();
    fExt[0] = lamb / 40;
    return fExt;
  };

  /// We artificially apply a Neumann load on the complete boundary
  Dune::BitSetVector<1> neumannVertices(gridView.size(2), true);

  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);

  Ikarus::FESettings feSettings;
  feSettings.addOrAssign("youngs_modulus", 1000.0);
  feSettings.addOrAssign("poissons_ratio", 0.3);
  feSettings.addOrAssign("thickness", 0.1);
  feSettings.addOrAssign("simulationFlag", 0);
  using KLSHELL = Ikarus::KirchhoffLoveShell<decltype(basis)>;
  KLSHELL fe(basis, *element, feSettings, volumeLoad, &neumannBoundary, neumannBoundaryLoad);
  using AutoDiffBasedFE = Ikarus::AutoDiffFE<KLSHELL>;
  AutoDiffBasedFE feAutoDiff(fe);

  Eigen::VectorXd d;
  d.setRandom(nDOF);
  double lambda = 7.3;

  auto req = Ikarus::FErequirements().addAffordance(Ikarus::AffordanceCollections::elastoStatics);
  req.insertGlobalSolution(Ikarus::FESolutions::displacement, d)
      .insertParameter(Ikarus::FEParameter::loadfactor, lambda);

  Eigen::MatrixXd K, KAutoDiff;
  K.setZero(nDOF, nDOF);
  KAutoDiff.setZero(nDOF, nDOF);

  Eigen::VectorXd R, RAutoDiff;
  R.setZero(nDOF);
  RAutoDiff.setZero(nDOF);

//  fe.calculateMatrix(req, K);
//  feAutoDiff.calculateMatrix(req, KAutoDiff);

  fe.calculateVector(req, R);
  feAutoDiff.calculateVector(req, RAutoDiff);

//  t.check(K.isApprox(KAutoDiff, tol),
//          "Mismatch between the stiffness matrices obtained from explicit implementation and the one based on "
//          "automatic differentiation:" << K <<"\n"<< KAutoDiff<<"\n"<< K-KAutoDiff);

  t.check(R.isApprox(RAutoDiff, tol))<<
      "Mismatch between the residual vectors obtained from explicit implementation and the one based on "
      "automatic differentiation:" << R <<"\n"<< RAutoDiff<<"\n"<< R-RAutoDiff;

  t.check(Dune::FloatCmp::eq(fe.calculateScalar(req), feAutoDiff.calculateScalar(req), tol),
          "Mismatch between the energies obtained from explicit implementation and the one based on "
          "automatic differentiation");

  return t;
}

auto NonLinearElasticityLoadControlNRandTRforKLShell() {
  TestSuite t("NonLinearElasticityLoadControlNRandTRforKLShell ");
  constexpr auto dimworld        = 3;
  const std::array<int, 2> order = {1, 1};

  const std::array<std::vector<double>, 2> knotSpans = {{{0, 0, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0, 0}, .w = 1}, {.p = {10, 0, 0}, .w = 1}}, {{.p = {0, 2, 0}, .w = 1}, {.p = {10, 2, 0}, .w = 1}}};

  std::array<int, 2> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<2, dimworld>;

  Dune::IGA::NURBSPatchData<2, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  for (int i = 0; i < 2; ++i)
    patchData = degreeElevate(patchData, i, 1);

  Dune::ParameterTree parameterSet;
  Dune::ParameterTreeParser::readINITree("/workspaces/lex/Dokumente/ikarusVSC/ikarus/tests/src/shell.parset",
                                         parameterSet);

  const auto E              = parameterSet.get<double>("E");
  const auto nu             = parameterSet.get<double>("nu");
  const auto thickness      = parameterSet.get<double>("thickness");
  const auto loadFactor     = parameterSet.get<double>("loadFactor");
  const auto simulationFlag = parameterSet.get<int>("simulationFlag");
  const auto refine         = parameterSet.get<int>("refine");
  auto grid = std::make_shared<Grid>(patchData);

  grid->globalRefineInDirection(0,refine);
  auto gridView = grid->leafGridView();

  using GridView = decltype(gridView);
  using namespace Ikarus;
  using namespace Dune::Functions::BasisFactory;

  // const double E         = 100000;
  // const double nu        = 0.0;
  // const double thickness = 0.001;
  Ikarus::FESettings feSettings;
  feSettings.addOrAssign("youngs_modulus", E);
  feSettings.addOrAssign("poissons_ratio", nu);
  feSettings.addOrAssign("thickness", thickness);
  feSettings.addOrAssign("simulationFlag", simulationFlag);
  auto basis      = Ikarus::makeBasis(gridView, power<3>(nurbs(), FlatInterleaved()));
  auto volumeLoad = [thickness, loadFactor]([[maybe_unused]] auto& globalCoord, auto& lamb) {
    Eigen::Vector3d fext;
    fext.setZero();
    //    fext[1]= 2 * Dune::power(thickness, 3) * lamb / 10;
    fext[2] = 2 * Dune::power(thickness, 3) * lamb * loadFactor;
    return fext;
  };

  using ElementType = Ikarus::AutoDiffFE<Ikarus::KirchhoffLoveShell<decltype(basis)>>;
  std::vector<ElementType> fes;

  for (auto& element : elements(gridView))
    fes.emplace_back(basis, element, feSettings, volumeLoad);

  auto basisP = std::make_shared<const decltype(basis)>(basis);
  Ikarus::DirichletValues dirichletValues(basisP->flat());

  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[0]) < 1e-8) dirichletFlags[localView.index(localIndex)] = true;
  });

  dirichletValues.fixDOFs([&](auto& basis, auto&& dirichletFlags) {
    Dune::Functions::forEachBoundaryDOF(Dune::Functions::subspaceBasis(basis, 2),
                                        [&](auto&& localIndex, auto&& localView, auto&& intersection) {
                                          if (std::abs(intersection.geometry().center()[0]) > 10 - 1e-8)
                                            dirichletFlags[localView.index(localIndex)] = true;
                                        });
  });

  dirichletValues.fixDOFs([&](auto& basis, auto&& dirichletFlags) {
    Dune::Functions::forEachBoundaryDOF(Dune::Functions::subspaceBasis(basis, 1),
                                        [&](auto&& localIndex, auto&& localView, auto&& intersection) {
                                          if (std::abs(intersection.geometry().center()[0]) > 10 - 1e-8)
                                            dirichletFlags[localView.index(localIndex)] = true;
                                        });
  });

  auto sparseAssembler = SparseFlatAssembler(fes, dirichletValues);

  Eigen::VectorXd d;
  d.setZero(basis.flat().size());
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

  auto nonLinOp = Ikarus::NonLinearOperator(functions(residualFunction, KFunction), parameter(d, lambda));

  const double gradTol = 1e-16;

  // auto tr = Ikarus::makeTrustRegion(nonLinOp);
  // tr->setup({.verbosity = 1,
  //            .maxiter   = 1000,
  //            .grad_tol  = gradTol,
  //            .corr_tol  = 1e-8,  // everything should converge to the gradient tolerance
  //            .useRand   = false,
  //            .rho_reg   = 1e8,
  //            .Delta0    = 1});
  auto linSolver               = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::sd_UmfPackLU);
  auto tr                      = Ikarus::makeNewtonRaphson(nonLinOp, std::move(linSolver));
  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
  tr->subscribeAll(nonLinearSolverObserver);

  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>>>(
      basis.flat(), d, 2);
  vtkWriter->setFileNamePrefix("Test2DSolid");
  vtkWriter->setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 3);

  auto lc = Ikarus::LoadControl(tr, 1, {0, 1});
  lc.subscribeAll(vtkWriter);
  const auto controlInfo = lc.run();

  std::cout << std::setprecision(16) << std::ranges::max(d) << std::endl;
  t.check(Dune::FloatCmp::eq(0.2957393081676369, std::ranges::max(d)))
      << std::setprecision(16) << "The maximum displacement is " << std::ranges::max(d);
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  //  const double E             = materialParameters.get<double>("E");
  //  const double nu            = materialParameters.get<double>("nu");
  checkFEByAutoDiff();
//  NonLinearElasticityLoadControlNRandTRforKLShell();
}
