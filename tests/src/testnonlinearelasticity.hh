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

#include "spdlog/spdlog.h"

#include <Eigen/Core>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mechanics/nonlinearelastic.hh>
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

struct OwnResultFunction {
  template <typename ElementTypeT, typename FERequirements, int size, typename ScalarTypeT>
  double operator()(const ElementTypeT& /* fe */, const Ikarus::ResultRequirements<FERequirements>& /* req */,
                    const Dune::FieldVector<ScalarTypeT, size>& /* pos */, [[maybe_unused]] int /* comp */) const {
    return 7;
  }
  static std::string name() { return "Seven"; }
  static int ncomps() { return 1; }
};

template <typename Grid, typename Material>
auto NonLinearElasticityLoadControlNRandTR(const Material& mat) {
  TestSuite t("NonLinearElasticityLoadControlNRandTR " + Dune::className(Grid{}) + " " + Dune::className(mat));
  auto grid     = createGrid<Grid>();
  auto gridView = grid->leafGridView();

  using GridView = decltype(gridView);
  using namespace Ikarus;
  using namespace Dune::Functions::BasisFactory;

  auto basis      = Ikarus::makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));
  auto volumeLoad = []([[maybe_unused]] auto& globalCoord, auto& lamb) {
    Eigen::Vector2d fext;
    fext.setZero();
    fext[1] = 2 * lamb;
    fext[0] = lamb;
    return fext;
  };

  auto reducedMat = planeStress(mat, 1e-8);

  using ElementType = Ikarus::NonLinearElastic<decltype(basis), decltype(reducedMat)>;
  std::vector<ElementType> fes;

  for (auto& element : elements(gridView))
    fes.emplace_back(basis, element, reducedMat, volumeLoad);

  auto basisP = std::make_shared<const decltype(basis)>(basis);
  Ikarus::DirichletValues dirichletValues(basisP->flat());

  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[1]) < 1e-8) dirichletFlags[localView.index(localIndex)] = true;
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
  vtkWriter->setFileNamePrefix("Test2DSolid");
  vtkWriter->setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 2);

  auto lc = Ikarus::LoadControl(tr, 1, {0, 50});
  lc.subscribeAll(vtkWriter);
  const auto controlInfo = lc.run();
  nonLinOp.template update<0>();
  const auto maxDisp = std::ranges::max(d);
  double energyExpected;
  if (std::is_same_v<Grid, Grids::Yasp>)
    energyExpected = -2.9605187645668578078;
  else if (std::is_same_v<Grid, Grids::Alu>)
    energyExpected = -2.9530594665063669702;
  else /* std::is_same_v<Grid, Grids::Iga> */
    energyExpected = -1.4533281398929942529;

  double maxDispExpected;
  if (std::is_same_v<Grid, Grids::Yasp>)
    maxDispExpected = 0.11293260007792008115;
  else if (std::is_same_v<Grid, Grids::Alu>)
    maxDispExpected = 0.1123397197762363714;
  else /* std::is_same_v<Grid, Grids::Iga> */
    maxDispExpected = 0.061647849558021668159;

  std::cout << std::setprecision(20) << nonLinOp.value() << std::endl;
  std::cout << "Maxdisp: " << maxDisp << std::endl;
  if constexpr (std::is_same_v<Material, Ikarus::StVenantKirchhoff>) {
    t.check(Dune::FloatCmp::eq(energyExpected, nonLinOp.value()), "energyExpected == nonLinOp.value()")
        << "energyExpected: " << energyExpected << "\nnonLinOp.value(): " << nonLinOp.value();

    t.check(std::abs(maxDispExpected - maxDisp) < 1e-12, "maxDispExpected-maxDisp")
        << "\nmaxDispExpected: \n"
        << maxDispExpected << "\nmaxDisp: \n"
        << maxDisp;
  } else {  // using a Neohooke material yields a lower energy and larger displacements
    t.check(Dune::FloatCmp::gt(energyExpected, nonLinOp.value()), "energyExpected > nonLinOp.value()")
        << "energyExpected: " << energyExpected << "\nnonLinOp.value(): " << nonLinOp.value();

    t.check(maxDispExpected < maxDisp, "maxDispExpected<maxDisp") << "maxDispExpected: \n"
                                                                  << maxDispExpected << "\nmaxDisp: \n"
                                                                  << maxDisp;
  }

  Dune::Vtk::VtkWriter<GridView> vtkWriter2(gridView);
  auto resReq = Ikarus::ResultRequirements()
                    .insertGlobalSolution(Ikarus::FESolutions::displacement, d)
                    .insertParameter(Ikarus::FEParameter::loadfactor, lambda)
                    .addResultRequest(ResultType::PK2Stress);
  auto resultFunction = std::make_shared<ResultFunction<ElementType>>(&fes, resReq);

  t.check(resultFunction->name() == "PK2Stress")
      << "Test resultName: " << resultFunction->name() << "should be PK2Stress";
  t.check(resultFunction->ncomps() == 3) << "Test result comps: " << resultFunction->ncomps() << "should be 3";

  vtkWriter2.addPointData(Dune::Vtk::Function<GridView>(resultFunction));

  auto resultFunction2 = ResultFunction(&fes, resReq, ResultEvaluators::PrincipalStress{});

  auto resultFunction2S = std::make_shared<decltype(resultFunction2)>(resultFunction2);
  t.check(resultFunction2S->name() == "PrincipalStress")
      << "Test resultName: " << resultFunction2S->name() << "should be PrincipalStress";
  t.check(resultFunction2S->ncomps() == 2) << "Test result comps: " << resultFunction2S->ncomps() << "should be 2";
  vtkWriter2.addPointData(Dune::Vtk::Function<GridView>(resultFunction2S));
  auto resultFunction3  = ResultFunction(&fes, resReq, ResultEvaluators::VonMises{});
  auto resultFunction3S = std::make_shared<decltype(resultFunction3)>(resultFunction3);
  t.check(resultFunction3S->name() == "VonMises")
      << "Test resultName: " << resultFunction2S->name() << "should be VonMises";
  t.check(resultFunction3S->ncomps() == 1) << "Test result comps: " << resultFunction2S->ncomps() << "should be 1";
  vtkWriter2.addPointData(Dune::Vtk::Function<GridView>(resultFunction3S));

  auto resultFunction4  = ResultFunction(&fes, resReq, OwnResultFunction{});
  auto resultFunction4S = std::make_shared<decltype(resultFunction4)>(resultFunction4);

  vtkWriter2.addPointData(Dune::Vtk::Function<GridView>(resultFunction4S));

  vtkWriter2.write("EndResult" + Dune::className<Grid>());
  auto resReq5 = Ikarus::ResultRequirements()
                     .insertGlobalSolution(Ikarus::FESolutions::displacement, d)
                     .insertParameter(Ikarus::FEParameter::loadfactor, lambda)
                     .addResultRequest(ResultType::linearStress);
  auto resultFunction5 = std::make_shared<ResultFunction<ElementType>>(&fes, resReq5);
  try {
    vtkWriter2.addPointData(Dune::Vtk::Function<GridView>(resultFunction5));
    t.check(false) << "resultFunction5 should have failed for requesting linearStress here";
  } catch (const Dune::NotImplemented&) {
  }

  nonLinOp.template update<1>();
  t.check(controlInfo.success, "Successful result");
  t.check(gradTol >= nonLinOp.derivative().norm(), "Gradient Tolerance should be larger than actual tolerance");
  return t;
}

/**
 * For zero displacements, the stiffness matrix from geometrically linear and non-linear case should be the same,
 * there by ensuring correctness of the Green-Lagrange strain
 */
template <int gridDim, typename Material>
auto GreenLagrangeStrainTest(const Material& mat) {
  static_assert(gridDim == 2 or gridDim == 3, "Only supported for the 2D and 3D case");
  TestSuite t("GreenLagrangeStrainTest with gridDim = " + std::to_string(gridDim));

  const double tol = 1e-8;
  auto grid        = createUGGridFromCorners<gridDim>(CornerDistortionFlag::randomlyDistorted);
  auto gridView    = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis   = Ikarus::makeBasis(gridView, power<gridDim>(lagrange<1>()));
  auto element = gridView.template begin<0>();
  auto nDOF    = basis.flat().size();

  using NonLinearElasticity = Ikarus::NonLinearElastic<decltype(basis), decltype(mat)>;
  using LinearElasticity    = Ikarus::LinearElastic<decltype(basis)>;

  NonLinearElasticity fe(basis, *element, mat);
  LinearElasticity feLE(basis, *element, 1000.0, 0.0);

  Eigen::VectorXd d;
  d.setZero(nDOF);
  double lambda = 0.0;

  auto req = Ikarus::FErequirements().addAffordance(Ikarus::AffordanceCollections::elastoStatics);
  req.insertGlobalSolution(Ikarus::FESolutions::displacement, d)
      .insertParameter(Ikarus::FEParameter::loadfactor, lambda);

  Eigen::MatrixXd K, KLE;
  K.setZero(nDOF, nDOF);
  KLE.setZero(nDOF, nDOF);

  fe.calculateMatrix(req, K);
  feLE.calculateMatrix(req, KLE);

  t.check(K.isApprox(KLE, tol),
          "Mismatch between linear and non-linear stiffness matrix for zero displacements with gridDim = "
              + std::to_string(gridDim));
  return t;
}

template <typename Material, int gridDim = 2>
auto SingleElementTest(const Material& mat) {
  static_assert(gridDim == 2, "Single element test is applicable only for the 2D case");
  TestSuite t("Single element test for non-linear Q1 element");

  auto grid     = createUGGridFromCorners<gridDim>(CornerDistortionFlag::fixedDistorted);
  auto gridView = grid->leafGridView();

  using namespace Ikarus;

  using namespace Dune::Functions::BasisFactory;
  auto basis       = Ikarus::makeBasis(gridView, power<gridDim>(lagrange<1>()));
  auto element     = gridView.template begin<0>();
  auto nDOF        = basis.flat().size();
  const double tol = 1e-10;

  using NonLinearElastic = Ikarus::NonLinearElastic<decltype(basis), decltype(mat)>;
  NonLinearElastic fe(basis, *element, mat);

  Eigen::VectorXd d;
  d.setZero(nDOF);
  double lambda = 0.0;

  d << 2, 4, 3.25, -1.2, 0.003, 6, 3, 2.864;

  auto req = Ikarus::FErequirements().addAffordance(Ikarus::AffordanceCollections::elastoStatics);
  req.insertGlobalSolution(Ikarus::FESolutions::displacement, d)
      .insertParameter(Ikarus::FEParameter::loadfactor, lambda);

  Eigen::MatrixXd K;
  K.setZero(nDOF, nDOF);
  fe.calculateMatrix(req, K);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> essaK;
  essaK.compute(K);
  auto eigenValuesComputed = essaK.eigenvalues();

  /// The eigen values are applicable only if 2x2 Gauss integration points are used
  Eigen::VectorXd eigenValuesExpected;
  eigenValuesExpected.setZero(basis.flat().size());
  eigenValuesExpected << 1e-16, 1e-16, 1845.6296388251504753, 14192.4707553121224317, 19964.32719133414782,
      29973.7943273325380486, 46641.183728849332812, 95447.6156712376251918;
  for (size_t i = 0; i < basis.flat().size(); ++i) {
    if (abs(eigenValuesComputed[i]) > tol) {
      t.check(Dune::FloatCmp::eq(abs(eigenValuesComputed[i]), eigenValuesExpected[i], tol),
              "Mismatch in the " + std::to_string(i + 1)
                  + "-th eigen value in single element test for four node non-linear 2D element");
    }
  }
  return t;
}

template <int gridDim, typename Material>
auto checkFEByAutoDiff(const Material& mat) {
  TestSuite t("Check calculateScalarImpl() and calculateVectorImpl() by Automatic Differentiation for gridDim = "
              + std::to_string(gridDim));

  auto grid     = createUGGridFromCorners<gridDim>(CornerDistortionFlag::randomlyDistorted);
  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis       = Ikarus::makeBasis(gridView, power<gridDim>(lagrange<1>()));
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

  using NonLinearElasticity = Ikarus::NonLinearElastic<decltype(basis), decltype(mat)>;
  NonLinearElasticity fe(basis, *element, mat, volumeLoad, &neumannBoundary, neumannBoundaryLoad);
  using AutoDiffBasedFE = Ikarus::AutoDiffFE<NonLinearElasticity>;
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

  fe.calculateMatrix(req, K);
  feAutoDiff.calculateMatrix(req, KAutoDiff);

  fe.calculateVector(req, R);
  feAutoDiff.calculateVector(req, RAutoDiff);

  t.check(K.isApprox(KAutoDiff, tol),
          "Mismatch between the stiffness matrices obtained from explicit implementation and the one based on "
          "automatic differentiation with gridDim = "
              + std::to_string(gridDim));

  t.check(R.isApprox(RAutoDiff, tol),
          "Mismatch between the residual vectors obtained from explicit implementation and the one based on "
          "automatic differentiation with gridDim = "
              + std::to_string(gridDim));

  t.check(Dune::FloatCmp::eq(fe.calculateScalar(req), feAutoDiff.calculateScalar(req), tol),
          "Mismatch between the energies obtained from explicit implementation and the one based on "
          "automatic differentiation with gridDim = "
              + std::to_string(gridDim));

  return t;
}
