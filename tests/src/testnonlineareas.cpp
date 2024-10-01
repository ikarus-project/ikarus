// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcommon.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <ikarus/finiteelements/autodifffe.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/nonlinearelastic.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>

using namespace Ikarus;
using namespace Dune::Indices;
using Dune::TestSuite;

template <int gridDim, typename GV, typename MAT, typename BH, typename VectorType>
auto nonlinearEASAutoDiffTest(const GV& gridView, const MAT& mat, const BH& basis, VectorType& d,
                              const int numberOfEASParameters) {
  TestSuite t("Nonlinear EAS AutoDiff with gridDim = " + std::to_string(gridDim) +
              " and numberOfEASParameters = " + std::to_string(numberOfEASParameters));
  double lambda = 7.3;
  auto sk       = skills(nonLinearElastic(mat), eas(numberOfEASParameters));
  auto fe       = Ikarus::makeFE(basis, sk);
  using FE      = decltype(fe);
  auto req      = typename FE::Requirement();
  req.insertGlobalSolution(d).insertParameter(lambda);

  auto affordance = AffordanceCollections::elastoStatics;

  for (auto element : elements(gridView)) {
    auto localView = basis.flat().localView();
    localView.bind(element);
    auto nDOF        = localView.size();
    const double tol = 1e-10;

    fe.bind(element);
    updateState(fe, req, d); // here d = correction vector (DeltaD)

    const std::string feClassName = Dune::className(fe);

    using AutoDiffBasedFE = Ikarus::AutoDiffFE<decltype(fe), true>;
    AutoDiffBasedFE feAutoDiff(fe);

    Eigen::MatrixXd K, KAutoDiff;
    K.setZero(nDOF, nDOF);
    KAutoDiff.setZero(nDOF, nDOF);

    calculateMatrix(fe, req, affordance.matrixAffordance(), K);
    calculateMatrix(feAutoDiff, req, affordance.matrixAffordance(), KAutoDiff);

    t.check(fe.numberOfEASParameters() == feAutoDiff.realFE().numberOfEASParameters())
        << "Number of EAS parameters for FE(" << fe.numberOfEASParameters()
        << ") and number of EAS parameters for AutodiffFE(" << feAutoDiff.realFE().numberOfEASParameters()
        << ") are not equal";

    t.check(fe.alpha().isApprox(feAutoDiff.alpha(), tol), "Mismatch in alpha.")
        << "alpha is\n"
        << fe.alpha().transpose() << "\n alphaAutoDiff is\n"
        << feAutoDiff.alpha();

    t.check(K.isApprox(KAutoDiff, tol),
            "Mismatch between the stiffness matrices obtained from explicit implementation and the one based on "
            "automatic differentiation.")
        << "K is \n"
        << K << "\nKAutoDiff is \n"
        << KAutoDiff << "\nThe difference is\n"
        << (K - KAutoDiff);

    if (numberOfEASParameters != 0) {
      t.checkThrow<Dune::NotImplemented>(
          [&]() {
            Eigen::VectorXd RAutoDiff;
            RAutoDiff.setZero(nDOF);
            calculateVector(feAutoDiff, req, affordance.vectorAffordance(), RAutoDiff);
          },
          "calculateVector with AutoDiff should throw a Dune::NotImplemented");

      t.checkThrow<Dune::NotImplemented>(
          [&]() { auto energyAutoDiff = calculateScalar(feAutoDiff, req, affordance.scalarAffordance()); },
          "calculateScalar with AutoDiff should throw a Dune::NotImplemented");
    }
  }
  return t;
}

template <int gridDim, typename TestSuitType>
void easAutoDiffTest(const TestSuitType& t) {
  std::array<int, (gridDim == 2) ? 4 : 3> easParameters;
  if constexpr (gridDim == 2)
    easParameters = {0, 4, 5, 7};
  else
    easParameters = {0, 9, 21};

  auto grid = createUGGridFromCorners<gridDim>(CornerDistortionFlag::randomlyDistorted);
  grid->globalRefine(2);
  auto gridView     = grid->leafGridView();
  auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = 100.0, .nu = 0.3});
  StVenantKirchhoff matSVK(matParameter);
  auto reducedMat = planeStress(matSVK);

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<gridDim>(lagrange<1>()));
  Eigen::VectorXd d;

  auto autoDiffTestFunctor = [&](int numberOfEASParameters) {
    if constexpr (gridDim == 2)
      nonlinearEASAutoDiffTest<gridDim>(gridView, reducedMat, basis, d, numberOfEASParameters);
    else
      nonlinearEASAutoDiffTest<gridDim>(gridView, matSVK, basis, d, numberOfEASParameters);
  };

  for (const int numberOfEASParameters : easParameters) {
    d.setZero(basis.flat().size());
    autoDiffTestFunctor(numberOfEASParameters);
    d.setRandom(basis.flat().size());
    autoDiffTestFunctor(numberOfEASParameters);
  }
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t("Nonlinear EAS Test");
  easAutoDiffTest<2>(t);
  easAutoDiffTest<3>(t);
  return t.exit();
}