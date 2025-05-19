// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "checkfebyautodiff.hh"
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
#include <ikarus/finiteelements/autodifffe.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/nonlinearelastic.hh>
#include <ikarus/solver/eigenvaluesolver/generalizedeigensolver.hh>
#include <ikarus/solver/eigenvaluesolver/generalizedeigensolverfactory.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverfactory.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/functionhelper.hh>
#include <ikarus/utils/init.hh>

using namespace Ikarus;
using Dune::TestSuite;

template <int gridDim>
auto assumedStressParametersFunc() {
  std::array<int, (gridDim == 2) ? 1 : 3> asParameters;
  if constexpr (gridDim == 2)
    asParameters = {5};
  else
    asParameters = {18, 24, 30};
  return asParameters;
}

template <int gridDim, typename TestSuitType, typename MAT>
void asAutoDiffTest(TestSuitType& t, const MAT& mat) {
  using namespace Dune::Functions::BasisFactory;

  auto grid               = createUGGridFromCorners<gridDim>(CornerDistortionFlag::unDistorted);
  auto gridView           = grid->leafGridView();
  const auto asParameters = assumedStressParametersFunc<gridDim>();

  for (const int numberOfASParameters : asParameters) {
    auto sk = skills(nonLinearElastic(mat), assumedStress<PS::PK2Stress>(numberOfASParameters));
    t.subTest(checkFESByAutoDiff(gridView, power<gridDim>(lagrange<1>(), FlatInterleaved{}), sk,
                                 AffordanceCollections::elastoStatics,
                                 " (numberOfInternalVariables = " + std::to_string(numberOfASParameters) + ")"));
  }
}

auto singleElementTest() {
  constexpr int gridDim                   = 3;
  constexpr int numberOfInternalVariables = 24;
  TestSuite t("Single element test for non-linear H1S24 element");
  using namespace Ikarus;

  auto grid     = createUGGridFromCorners<gridDim>(CornerDistortionFlag::fixedDistorted);
  auto gridView = grid->leafGridView();

  using namespace Ikarus;

  using namespace Dune::Functions::BasisFactory;
  auto basis        = Ikarus::makeBasis(gridView, power<gridDim>(lagrange<1>(), FlatInterleaved{}));
  auto element      = gridView.template begin<0>();
  auto nDOF         = basis.flat().size();
  const double tol  = 1e-10;
  auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = 100.0, .nu = 0.3});
  Materials::StVenantKirchhoff mat(matParameter);

  auto fe = makeFE(basis, skills(nonLinearElastic(mat), assumedStress<PS::PK2Stress>(numberOfInternalVariables)));
  fe.bind(*element);

  Eigen::VectorXd d;
  d.setZero(nDOF);
  double lambda = 0.0;

  auto req = typename decltype(fe)::Requirement(d, lambda);

  Eigen::MatrixXd K;
  Eigen::VectorXd R;
  R.setZero(nDOF);
  K.setZero(nDOF, nDOF);
  fe.updateState(req, d);
  const auto beta = fe.internalVariable();
  calculateVector(fe, req, Ikarus::VectorAffordance::forces, R);
  calculateMatrix(fe, req, Ikarus::MatrixAffordance::stiffness, K);
  auto essaK = Ikarus::makeIdentitySymEigenSolver<EigenValueSolverType::Eigen>(K);
  essaK.compute();
  auto eigenValuesComputed = essaK.eigenvalues();

  /// The expected values are applicable only if 2x2 Gauss integration points are used
  Eigen::VectorXd eigenValuesExpected, RExpected, betaExpected;
  eigenValuesExpected.setZero(nDOF);
  RExpected.setZero(nDOF);
  betaExpected.setZero(numberOfInternalVariables);

  eigenValuesExpected << 1e-16, 1e-16, 1e-16, 1e-16, 1e-16, 1e-16, 5.39121162342041660110, 7.05781676246691892760,
      9.36884335567226036850, 9.56659918651701746000, 10.40366553068070634300, 11.331509846741802377,
      13.19114818364121971800, 13.46399600344537256400, 21.12738335399774334400, 24.011092030411375791,
      24.89333220181543055200, 26.49727774428625148200, 35.08628709218272357900, 35.69141870944951953900,
      39.83534197522907002700, 41.80698841900341171400, 41.83687685676100413200, 126.61728448545971601000;

  checkApproxVectors(t, R, RExpected, testLocation() + "\nIncorrect residual vectors.", tol);
  checkApproxVectors(t, beta, betaExpected, testLocation() + "\nIncorrect internal parameter beta.", tol);

  for (size_t i = 0; i < nDOF; ++i) {
    if (abs(eigenValuesComputed[i]) > tol) {
      checkScalars(t, abs(eigenValuesComputed[i]), eigenValuesExpected[i],
                   " Mismatch in the " + std::to_string(i + 1) + "-th eigen value.", tol);
    }
  }
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t("Nonlinear Assumed Stress Test");
  auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = 100.0, .nu = 0.3});
  Materials::StVenantKirchhoff matSVK(matParameter);
  Materials::NeoHooke matNH(matParameter);
  auto matBK         = Materials::makeBlatzKo(40.0);
  auto reducedMatSVK = planeStrain(matSVK);
  auto reducedMatNH  = planeStrain(matNH);
  auto reducedMatBK  = planeStrain(matBK);

  asAutoDiffTest<2>(t, reducedMatSVK);
  asAutoDiffTest<3>(t, matSVK);

  asAutoDiffTest<2>(t, reducedMatNH);
  asAutoDiffTest<3>(t, matNH);

  asAutoDiffTest<2>(t, reducedMatBK);
  asAutoDiffTest<3>(t, matBK);

  t.subTest(singleElementTest());

  return t.exit();
}
