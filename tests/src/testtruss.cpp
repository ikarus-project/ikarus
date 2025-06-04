// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "checkfebyautodiff.hh"
#include "testcommon.hh"
#include "testhelpers.hh"

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/controlroutines/pathfollowing.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/truss.hh>
#include <ikarus/solver/eigenvaluesolver/generalizedeigensolverfactory.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphsonwithscalarsubsidiaryfunction.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverfactory.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/listener/controllogger.hh>
#include <ikarus/utils/listener/controlvtkwriter.hh>
#include <ikarus/utils/listener/genericlistener.hh>
#include <ikarus/utils/listener/nonlinearsolverlogger.hh>

using namespace Ikarus;
using Dune::TestSuite;

Eigen::Matrix<double, 4, 4> rotationMatrix2D(double theta) {
  const auto c = std::cos(theta);
  const auto s = std::sin(theta);
  Eigen::Matrix2d rotMat;
  rotMat.setZero();
  rotMat << c, s, -s, c;
  Eigen::Matrix<double, 4, 4> Q;
  Q.setZero();
  Q.block<2, 2>(0, 0) = rotMat;
  Q.block<2, 2>(2, 2) = rotMat;
  return Q;
}

static auto vonMisesTrussTest() {
  TestSuite t("vonMisesTrussTest");
  std::cout << "Started: " << t.name() << std::endl;
  /// Construct grid
  constexpr double h = 1.0;
  constexpr double L = 10.0;
  auto grid          = createGrid<Grids::OneDFoamGridIn2D>();
  auto gridView      = grid->leafGridView();

  /// Construct basis
  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved{}));

  /// Create finite elements
  constexpr double E   = 30000;
  constexpr double A   = 0.1;
  constexpr double tol = 1e-10;
  auto sk              = skills(truss(E, A));
  using FEType         = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;
  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  /// Collect dirichlet nodes
  auto basisP = std::make_shared<const decltype(basis)>(basis);
  DirichletValues dirichletValues(basisP->flat());
  dirichletValues.fixBoundaryDOFs(
      [&](auto& dirichletFlags, auto&& globalIndex) { dirichletFlags[globalIndex] = true; });
  t.check(dirichletValues.fixedDOFsize() == 4, "Number of fixed DOFs is not equal to 4");

  /// Create assembler
  auto denseFlatAssembler = makeAssemblerManipulator(DenseFlatAssembler(fes, dirichletValues));

  /// Create non-linear operator
  double lambda = 0.0;
  Eigen::VectorXd d;
  d.setZero(basis.flat().size());

  auto req = FEType::Requirement(basis);
  denseFlatAssembler->bind(req, AffordanceCollections::elastoStatics, DBCOption::Full);

  Dune::FieldVector<double, 2> pointLoadPos{L, h};
  const auto globalIndicesAtPointLoadPos = utils::globalIndexFromGlobalPosition(basis.flat(), pointLoadPos);
  auto pointLoad                         = [&](const auto&, const auto&, auto, auto, Eigen::VectorXd& vec) -> void {
    vec[globalIndicesAtPointLoadPos[1]] -= -req.parameter();
  };
  denseFlatAssembler->bind(pointLoad);

  /// Test tangent stiffness matrix for element 0 for geometrically linear case
  const auto& ele0 = fes[0].gridElement();
  Eigen::MatrixXd ke0;
  ke0.setZero(4, 4);
  calculateMatrix(fes[0], req, MatrixAffordance::stiffness, ke0);
  const auto X1            = Dune::toEigen(ele0.geometry().corner(0));
  const auto X2            = Dune::toEigen(ele0.geometry().corner(1));
  const Eigen::Vector2d A1 = (X2 - X1).eval();
  const Eigen::Vector2d A0 = {1, 0}; // horizontal axis
  const double L0          = sqrt(A1.squaredNorm());
  const double theta0      = std::acos(A1.dot(A0) / (A0.norm() * A1.norm()));
  Eigen::Matrix<double, 4, 4> kbar;
  kbar.setZero();
  kbar(0, 0) = kbar(2, 2) = 1.0;
  kbar(0, 2) = kbar(2, 0) = -1.0;
  const auto Q            = rotationMatrix2D(theta0);
  const auto k0           = (Q.transpose() * kbar * Q) * E * A / L0;
  t.check(isApproxSame(k0, ke0, tol), "(k0 is not equal to ke0 for geometrically linear case.") << "Expected:\n"
                                                                                                << k0 << "\nActual:\n"
                                                                                                << ke0;

  /// Choose linear solver
  auto linSolver = LinearSolver(SolverTypeTag::d_LDLT);

  NewtonRaphsonConfig nrConfig({}, linSolver);
  NonlinearSolverFactory nrFactory(nrConfig);
  auto nr = nrFactory.create(denseFlatAssembler);

  /// Create Observer to write information of the non-linear solver
  auto f = DifferentiableFunctionFactory::op(denseFlatAssembler);

  t.check(utils::checkGradient(f, req, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << "Check gradient failed";
  t.check(utils::checkHessian(f, req, {.draw = false, .writeSlopeStatementIfFailed = true})) << "Check Hessian failed";

  constexpr int loadSteps = 10;

  Eigen::Matrix3Xd lambdaAndDisp;
  lambdaAndDisp.setZero(Eigen::NoChange, loadSteps + 1);

  /// Create Observer which writes vtk files when control routines messages
  auto vtkWriter = ControlSubsamplingVertexVTKWriter(basis.flat(), 2);
  vtkWriter.setFieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2);
  vtkWriter.setFileNamePrefix("vonMisesTruss");

  /// Create loadcontrol
  auto lc = ControlRoutineFactory::create(LoadControlConfig{loadSteps, 0.0, 0.5}, nr, denseFlatAssembler);

  auto nonLinearSolverLogger = NonLinearSolverLogger();
  auto controlLogger         = ControlLogger();

  nonLinearSolverLogger.subscribeTo(lc.nonLinearSolver());
  controlLogger.subscribeTo(lc);
  vtkWriter.subscribeTo(lc);

  auto lvkObserver = GenericListener(lc, ControlMessages::SOLUTION_CHANGED, [&](const auto& state) {
    const auto& d              = state.domain.globalSolution();
    const auto& lambda         = state.domain.parameter();
    int step                   = state.loadStep;
    lambdaAndDisp(0, step + 1) = lambda; // load factor
    lambdaAndDisp(1, step + 1) = d[2];   // horizontal displacement at center node
    lambdaAndDisp(2, step + 1) = d[3];   // vertical displacement at center node
  });

  /// Execute!
  auto controlInfo = lc.run(req);
  t.check(controlInfo.success == true) << "Load control failed to converge";

  Eigen::VectorXd lambdaVec = lambdaAndDisp.row(0);
  Eigen::VectorXd uVec      = lambdaAndDisp.row(1);
  Eigen::VectorXd vVec      = -lambdaAndDisp.row(2);

  t.check(Dune::FloatCmp::ne(vVec.tail(1)[0], 0.0, tol))
      << std::setprecision(16) << "Displacement value should not be zero. Actual:\t" << vVec.tail(1)[0];
  checkScalars(t, lambdaVec.tail(1)[0], 0.5, " Incorrect last load factor", tol);

  // return the load factor as a function of the vertical displacement
  auto analyticalLoadDisplacementCurve = [&](auto& v) {
    const double Ltruss = std::sqrt(h * h + L * L);
    return 2.0 * E * A * Dune::power(h, 3) / Dune::power(Ltruss, 3) *
           (v / h - 1.5 * Dune::power(v / h, 2) + 0.5 * Dune::power(v / h, 3));
  };

  for (std::size_t i = 0; i < lambdaAndDisp.cols(); ++i) {
    checkScalars(t, lambdaVec[i], analyticalLoadDisplacementCurve(vVec[i]), " Incorrect load factor", tol);
    checkScalars(t, uVec[i], 0.0, " Incorrect horizontal displacement", tol);
  }

  double expectedCauchyAxialForce = -2.785092479363964262;
  double expectedPK2AxialForce    = -2.787684078189828;

  // due to the symmetry of the problem, same axial forces are expected in both elements
  for (const auto& fe : fes) {
    const auto N    = fe.calculateAt<ResultTypes::cauchyAxialForce>(req, {0.5}).asVec()[0];
    const auto NPK2 = fe.calculateAt<ResultTypes::PK2AxialForce>(req, {0.5}).asVec()[0];
    checkScalars(t, N, expectedCauchyAxialForce, " Incorrect Cauchy Axial force", tol);
    checkScalars(t, NPK2, expectedPK2AxialForce, " Incorrect PK2 Axial force", tol);
  }

  for (auto i : Dune::range(0, loadSteps))
    fileExists(t, "vonMisesTruss" + std::to_string(i) + ".vtp");

  return t;
}

template <bool useArcLength = false>
static auto vonMisesTrussWithIDBCTest(const DBCOption dbcOption) {
  TestSuite t("vonMisesTrussTest with inhomogeneous Dirichlet BCs");
  std::cout << t.name() << " started with DBCOption = " << toString(dbcOption)
            << " and useArcLength = " << std::boolalpha << useArcLength << std::endl;

  /// Construct grid
  constexpr double h = 1.0;
  constexpr double L = 10.0;
  auto grid          = createGrid<Grids::OneDFoamGridIn2D>();
  auto gridView      = grid->leafGridView();

  /// Construct basis
  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved{}));

  /// Create finite elements
  constexpr double E   = 30000;
  constexpr double A   = 0.1;
  constexpr double tol = 1e-10;
  auto sk              = skills(truss(E, A));
  using FEType         = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;
  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  /// Collect dirichlet nodes
  auto basisP = std::make_shared<const decltype(basis)>(basis);
  DirichletValues dirichletValues(basisP->flat());
  dirichletValues.fixBoundaryDOFs(
      [&](auto& dirichletFlags, auto&& globalIndex) { dirichletFlags[globalIndex] = true; });
  t.check(dirichletValues.fixedDOFsize() == 4, "Number of constrained DOFs is not equal to 4");

  // Inhomogeneous Boundary Conditions
  auto inhomogeneousDisplacement = [L, h]<typename T>(const auto& globalCoord, const T& lambda) {
    Eigen::Vector<T, 2> localInhomogeneous;
    if (Dune::FloatCmp::eq(globalCoord[0], L) and Dune::FloatCmp::eq(globalCoord[1], h)) {
      localInhomogeneous[0] = 0;
      localInhomogeneous[1] = -lambda;
    } else
      localInhomogeneous.setZero();
    return localInhomogeneous;
  };

  dirichletValues.storeInhomogeneousBoundaryCondition(inhomogeneousDisplacement);
  t.check(dirichletValues.fixedDOFsize() == 5, "Number of constrained DOFs is not equal to 5");

  /// Create assembler
  auto denseFlatAssembler = makeDenseFlatAssembler(fes, dirichletValues);

  auto req = FEType::Requirement(basis);
  denseFlatAssembler->bind(req, AffordanceCollections::elastoStatics, dbcOption);
  auto linSolver          = LinearSolver(SolverTypeTag::d_LDLT);
  constexpr int loadSteps = 10;

  auto vtkWriter = ControlSubsamplingVertexVTKWriter(basis.flat(), 2);
  vtkWriter.setFieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2);
  vtkWriter.setFileNamePrefix("vonMisesTrussWithIDBC");

  auto controlRoutine = [&]() -> decltype(auto) {
    if constexpr (useArcLength) {
      NewtonRaphsonWithSubsidiaryFunctionConfig nrConfig({}, linSolver);
      auto nrFactory = NonlinearSolverFactory(nrConfig).withIDBCForceFunction(denseFlatAssembler);
      auto nr        = nrFactory.create(denseFlatAssembler);
      auto pft       = Ikarus::ArcLength{}; // Type of path following technique
      auto alc       = Ikarus::PathFollowing(nr, loadSteps, sqrt(0.5) / loadSteps, pft);
      return alc;
    } else { /// Create loadcontrol
      NewtonRaphsonConfig nrConfig({}, linSolver);
      auto nrFactory = NonlinearSolverFactory(nrConfig).withIDBCForceFunction(denseFlatAssembler);
      auto nr        = nrFactory.create(denseFlatAssembler);
      auto lc        = ControlRoutineFactory::create(LoadControlConfig{loadSteps, 0.0, 0.5}, nr, denseFlatAssembler);
      return lc;
    }
  };

  auto cr = controlRoutine();

  auto nonLinearSolverLogger = NonLinearSolverLogger();
  auto controlLogger         = ControlLogger();
  nonLinearSolverLogger.subscribeTo(cr.nonLinearSolver());
  controlLogger.subscribeTo(cr);
  vtkWriter.subscribeTo(cr);

  Eigen::Matrix3Xd lambdaAndDisp;
  lambdaAndDisp.setZero(Eigen::NoChange, loadSteps + 1);
  auto lvkObserver = GenericListener(cr, ControlMessages::SOLUTION_CHANGED, [&](const auto& state) {
    const auto& d              = state.domain.globalSolution();
    const auto& lambda         = state.domain.parameter();
    int step                   = state.loadStep;
    lambdaAndDisp(0, step + 1) = lambda; // load factor
    lambdaAndDisp(1, step + 1) = d[2];   // horizontal displacement at center node
    lambdaAndDisp(2, step + 1) = d[3];   // vertical displacement at center node
  });

  /// Execute!
  auto controlInfo = cr.run(req);
  t.check(controlInfo.success == true) << "Load control failed to converge";

  Eigen::VectorXd lambdaVec = lambdaAndDisp.row(0);
  Eigen::VectorXd uVec      = lambdaAndDisp.row(1);
  Eigen::VectorXd vVec      = lambdaAndDisp.row(2);

  std::vector<int> expectedIterations(loadSteps + 1, 0);

  t.check(Dune::FloatCmp::eq(vVec.tail(1)[0], -lambdaVec.tail(1)[0], tol))
      << std::setprecision(16) << "Displacement value should be equal to lambda. Actual:\t" << vVec.tail(1)[0];
  checkScalars(t, lambdaVec.tail(1)[0], 0.5, " Incorrect last load factor", tol);
  checkSolverInfos(t, expectedIterations, controlInfo, loadSteps + 1);

  // return the load factor as a function of the vertical displacement
  auto analyticalLoadDisplacementCurve = [&](auto& v) {
    return -v; // lambda = -v because an inhomogeneous Dirichlet BC is applied
  };

  for (std::size_t i = 0; i < lambdaAndDisp.cols(); ++i) {
    checkScalars(t, lambdaVec[i], analyticalLoadDisplacementCurve(vVec[i]), " Incorrect load factor", tol);
    checkScalars(t, uVec[i], 0.0, " Incorrect horizontal displacement", tol);
  }
  return t;
}

static auto truss3dTest() {
  TestSuite t("Truss3DTest");
  /// Construct grid
  auto grid     = createGrid<Grids::OneDFoamGridIn3D>();
  auto gridView = grid->leafGridView();

  /// Construct basis
  using namespace Dune::Functions::BasisFactory;
  auto basis      = Ikarus::makeBasis(gridView, power<3>(lagrange<1>(), FlatInterleaved{}));
  using LocalView = std::remove_cvref_t<decltype(basis.flat().localView())>;

  /// Create finite elements
  constexpr double E   = 2000;
  constexpr double A   = 1000;
  constexpr double tol = 1e-10;
  auto sk              = skills(truss(E, A));
  using FEType         = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;
  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  /// Collect dirichlet nodes
  auto basisP = std::make_shared<const decltype(basis)>(basis);
  DirichletValues dirichletValues(basisP->flat());
  dirichletValues.fixBoundaryDOFs(
      [&](auto& dirichletFlags, auto&& globalIndex) { dirichletFlags[globalIndex] = true; });
  std::array<LocalView::MultiIndex, 1> fixDOFs = {3}; // array containing the DOFs to be fixed at the center node
  for (const auto dof : fixDOFs)
    dirichletValues.setSingleDOF(dof, true);
  t.check(dirichletValues.fixedDOFsize() == 7, "Number of fixed DOFs is not equal to 7");

  /// Create assembler
  auto denseFlatAssembler = makeDenseFlatAssembler(fes, dirichletValues);

  double lambda = 1.0;
  Eigen::VectorXd dI;
  dI.setZero(basis.flat().size());

  auto req = FEType::Requirement();
  req.insertGlobalSolution(dI).insertParameter(lambda);
  auto& d = req.globalSolution();
  denseFlatAssembler->bind(req, AffordanceCollections::elastoStatics, DBCOption::Full);
  const auto& K = denseFlatAssembler->matrix();
  auto R        = denseFlatAssembler->vector();
  R[4]          = lambda;

  /// Choose linear solver
  auto linSolver = LinearSolver(SolverTypeTag::d_LDLT);
  linSolver.compute(K);
  linSolver.solve(d, -R);

  /// Compute eigenvalues of the stiffness matrix
  auto essaK = makeIdentitySymEigenSolver<EigenValueSolverType::Eigen>(K);
  essaK.compute();
  auto eigenValuesComputed = essaK.eigenvalues();

  Eigen::Vector2d axialForces = Eigen::Vector2d::Zero();
  for (int i = 0; const auto& fe : fes) {
    axialForces[i] = fe.calculateAt<ResultTypes::linearAxialForce>(req, {0.5}).asVec()[0];
    ++i;
  }

  Eigen::VectorXd expectedDisplacements, expectedEigenValues, expectedAxialForces;
  expectedDisplacements.setZero(basis.flat().size());
  expectedEigenValues.setZero(basis.flat().size());
  expectedAxialForces.setZero(2);
  expectedDisplacements << 0, 0, 0, 0, -0.0008521190304919996, 0.0002491175412207767, 0, 0, 0;
  expectedEigenValues << 1, 1, 1, 1, 1, 1, 1, 1081.004736091822, 734025.9250446925;
  expectedAxialForces << 13.78521130992862, -14.91038752769381;

  t.check(isApproxSame(expectedEigenValues, eigenValuesComputed, tol))
      << std::setprecision(16) << "Incorrect eigenvalues of K - Expected: " << expectedEigenValues.transpose()
      << " Actual: " << eigenValuesComputed.transpose();

  t.check(isApproxSame(expectedDisplacements, d, tol))
      << std::setprecision(16) << "Incorrect displacements - Expected: " << expectedDisplacements.transpose()
      << " Actual: " << d.transpose();

  t.check(isApproxSame(expectedAxialForces, axialForces, tol))
      << std::setprecision(16) << "Incorrect axial forces - Expected: " << expectedAxialForces.transpose()
      << " Actual: " << axialForces.transpose();

  return t;
}

template <int worldDim, typename GridType>
auto trussAutoDiffTest() {
  TestSuite t("Truss autodiff");
  using namespace Dune::Functions::BasisFactory;
  {
    auto grid     = createGrid<GridType>();
    auto gridView = grid->leafGridView();
    t.subTest(checkFESByAutoDiff(gridView, power<worldDim>(lagrange<1>(), FlatInterleaved{}),
                                 skills(truss(1000.0, 0.0956)), AffordanceCollections::elastoStatics));
  }
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t("TrussTest");
  t.subTest(trussAutoDiffTest<2, Grids::OneDFoamGridIn2D>());
  t.subTest(trussAutoDiffTest<3, Grids::OneDFoamGridIn3D>());
  t.subTest(vonMisesTrussTest());
  t.subTest(truss3dTest());
  t.subTest(vonMisesTrussWithIDBCTest(DBCOption::Full));
  t.subTest(vonMisesTrussWithIDBCTest(DBCOption::Reduced));
  t.subTest(vonMisesTrussWithIDBCTest<true>(DBCOption::Full));
  t.subTest(vonMisesTrussWithIDBCTest<true>(DBCOption::Reduced));
  return t.exit();
}
