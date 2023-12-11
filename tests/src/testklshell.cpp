// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcommon.hh"
#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#if HAVE_DUNE_IGA
#  include <dune/iga/nurbsbasis.hh>
#endif
#include "spdlog/spdlog.h"

#include <Eigen/Core>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/finiteelements/mechanics/kirchhoffloveshell.hh>
#include <ikarus/io/resultfunction.hh>
#include <ikarus/linearalgebra/dirichletvalues.hh>
#include <ikarus/linearalgebra/nonlinearoperator.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/solver/nonlinearsolver/trustregion.hh>
#include <ikarus/utils/algorithms.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controlvtkwriter.hh>

using Dune::TestSuite;

static auto NonLinearElasticityLoadControlNRandTRforKLShell() {
  TestSuite t("NonLinearElasticityLoadControlNRandTRforKLShell");
  constexpr auto dimworld        = 3;
  const std::array<int, 2> order = {1, 1};

  const std::array<std::vector<double>, 2> knotSpans = {{{0, 0, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0, 0}, .w = 1}, {.p = {10, 0, 0}, .w = 1}}, {{.p = {0, 2, 0}, .w = 1}, {.p = {10, 2, 0}, .w = 1}}};

  std::array<int, 2> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<2, dimworld>;

  Dune::IGA::NURBSPatchData<2, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  for (int i = 0; i < 2; ++i)
    patchData = degreeElevate(patchData, i, 1);

  auto grid = std::make_shared<Grid>(patchData);
  grid->globalRefine(2);
  auto gridView = grid->leafGridView();

  using namespace Ikarus;
  using namespace Dune::Functions::BasisFactory;
  const double E         = 1000;
  const double nu        = 0.0;
  const double thickness = 0.1;
  auto basis             = Ikarus::makeBasis(gridView, power<3>(nurbs(), FlatInterleaved()));
  auto volumeLoad        = [thickness]([[maybe_unused]] auto& globalCoord, auto& lamb) {
    Eigen::Vector3d fext;
    fext.setZero();
    fext[2] = 2 * Dune::power(thickness, 3) * lamb / 10;
    return fext;
  };

  using ElementType = Ikarus::AutoDiffFE<Ikarus::KirchhoffLoveShell<decltype(basis)>>;
  std::vector<ElementType> fes;

  for (auto& element : elements(gridView))
    fes.emplace_back(basis, element, E, nu, thickness, volumeLoad);

  auto basisP = std::make_shared<const decltype(basis)>(basis);
  Ikarus::DirichletValues dirichletValues(basisP->flat());

  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[0]) < 1e-8) dirichletFlags[localView.index(localIndex)] = true;
  });

  dirichletValues.fixDOFs([&](auto& basisL, auto&& dirichletFlags) {
    Dune::Functions::forEachBoundaryDOF(Dune::Functions::subspaceBasis(basisL, 2),
                                        [&](auto&& localIndex, auto&& localView, auto&& intersection) {
                                          if (std::abs(intersection.geometry().center()[0]) > 10 - 1e-8)
                                            dirichletFlags[localView.index(localIndex)] = true;
                                        });
  });

  dirichletValues.fixDOFs([&](auto& basisL, auto&& dirichletFlags) {
    Dune::Functions::forEachBoundaryDOF(Dune::Functions::subspaceBasis(basisL, 1),
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

  auto nonLinOp
      = Ikarus::NonLinearOperator(functions(energyFunction, residualFunction, KFunction), parameter(d, lambda));

  const double gradTol = 1e-8;

  auto tr = Ikarus::makeTrustRegion(nonLinOp);
  tr->setup({.verbosity = 1,
             .maxiter   = 1000,
             .grad_tol  = gradTol,
             .corr_tol  = 1e-16,  // everything should converge to the gradient tolerance
             .useRand   = false,
             .rho_reg   = 1e8,
             .Delta0    = 1});

  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>>>(
      basis.flat(), d, 2);
  vtkWriter->setFileNamePrefix("TestKLShell");
  vtkWriter->setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 3);

  auto lc = Ikarus::LoadControl(tr, 1, {0, 1});
  lc.subscribeAll(vtkWriter);
  const auto controlInfo = lc.run();

  t.check(controlInfo.success);
  std::cout << std::setprecision(16) << std::ranges::max(d) << std::endl;
  t.check(Dune::FloatCmp::eq(0.2957393081676369, std::ranges::max(d)))
      << std::setprecision(16) << "The maximum displacement is " << std::ranges::max(d);
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  NonLinearElasticityLoadControlNRandTRforKLShell();
}
