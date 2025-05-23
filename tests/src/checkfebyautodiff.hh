// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "testhelpers.hh"

#include <dune/common/float_cmp.hh>
#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/finiteelements/autodifffe.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverstate.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/broadcaster/broadcaster.hh>

template <typename GridView, typename BasisHandler, typename Skills, typename AffordanceColl, typename VectorType>
auto checkFESByAutoDiffImpl(const GridView& gridView, const BasisHandler& basis, Skills&& skills,
                            AffordanceColl affordance, VectorType& d, const std::string& testName = "",
                            double tol = 1e-10) {
  double lambda = 7.3;
  auto fe       = Ikarus::makeFE(basis, std::forward<Skills>(skills));
  using FE      = decltype(fe);

  auto req = typename FE::Requirement();
  VectorType dZero;
  dZero.setZero(d.size());

  req.insertGlobalSolution(dZero).insertParameter(lambda);
  Dune::TestSuite t("Check calculateMatrixImpl() and calculateVectorImpl() by Automatic Differentiation" + testName);
  for (auto element : elements(gridView)) {
    auto localView = basis.flat().localView();
    localView.bind(element);
    auto nDOF = localView.size();

    fe.bind(element);

    fe.updateState(req, d); // here d = correction vector (DeltaD)
    req.insertGlobalSolution(d).insertParameter(lambda);

    const std::string feClassName = Dune::className(fe);

    using AutoDiffBasedFE = Ikarus::AutoDiffFE<decltype(fe), true>;
    AutoDiffBasedFE feAutoDiff(fe);

    Eigen::MatrixXd K, KAutoDiff;
    K.setZero(nDOF, nDOF);
    KAutoDiff.setZero(nDOF, nDOF);

    Eigen::VectorXd R, RAutoDiff;
    R.setZero(nDOF);
    RAutoDiff.setZero(nDOF);

    calculateMatrix(fe, req, affordance.matrixAffordance(), K);
    calculateMatrix(feAutoDiff, req, affordance.matrixAffordance(), KAutoDiff);

    checkApproxMatrices(t, K, KAutoDiff, testLocation() + "\nIncorrect stiffness matrices." + feClassName, tol);
    checkSymmetricMatrix(t, K, tol, "K");
    checkSymmetricMatrix(t, KAutoDiff, tol, "KAutoDiff");

    if constexpr (fe.hasEAS() or fe.hasAssumedStress()) {
      checkScalars(t, fe.numberOfInternalVariables(), feAutoDiff.realFE().numberOfInternalVariables(),
                   "Incorrect number of EAS parameters");
      checkApproxVectors(t, fe.internalVariable(), feAutoDiff.internalVariable(),
                         testLocation() + "\nIncorrect internal variable alpha." + feClassName, tol);

      if constexpr (fe.hasEAS())
        if (fe.numberOfInternalVariables() != 0) {
          t.checkThrow<Dune::NotImplemented>(
              [&]() { calculateVector(feAutoDiff, req, affordance.vectorAffordance(), RAutoDiff); },
              "calculateVector with AutoDiff should throw a Dune::NotImplemented");

          t.checkThrow<Dune::NotImplemented>(
              [&]() { auto energyAutoDiff = calculateScalar(feAutoDiff, req, affordance.scalarAffordance()); },
              "calculateScalar with AutoDiff should throw a Dune::NotImplemented");
        }
    } else {
      calculateVector(fe, req, affordance.vectorAffordance(), R);
      calculateVector(feAutoDiff, req, affordance.vectorAffordance(), RAutoDiff);

      checkApproxVectors(t, R, RAutoDiff, testLocation() + "\nIncorrect residual vectors." + feClassName, tol);
      checkScalars(t, calculateScalar(fe, req, affordance.scalarAffordance()),
                   calculateScalar(feAutoDiff, req, affordance.scalarAffordance()),
                   testLocation() + "\nIncorrect Energies." + feClassName, tol);
    }

    using NRStateDummy = Ikarus::NonlinearSolverState<int, int>;
    struct DummyBroadcaster : public Ikarus::Broadcaster<Ikarus::ControlMessages, NRStateDummy>
    {
      using State = NRStateDummy;
    };
    auto broadcaster = DummyBroadcaster{};

    fe.template subscribeTo<Ikarus::ControlMessages>(broadcaster);
    feAutoDiff.template subscribeTo<Ikarus::ControlMessages>(broadcaster);
  }

  return t;
}

template <typename GridView, typename PreBasis, typename Skills, typename AffordanceColl>
auto checkFESByAutoDiff(const GridView& gridView, const PreBasis& pb, Skills&& skills, AffordanceColl affordance,
                        const std::string& testName = "", double tol = 1e-10) {
  Dune::TestSuite t("AutoDiff Test" + testName);
  auto basis = Ikarus::makeBasis(gridView, pb);
  Eigen::VectorXd d;
  d.setZero(basis.flat().dimension());
  t.subTest(checkFESByAutoDiffImpl(gridView, basis, skills, affordance, d, " Zero Displacements", tol));
  d.setRandom(basis.flat().dimension());
  d *= 0.01; // to avoid a negative determinant of deformation gradient
  t.subTest(checkFESByAutoDiffImpl(gridView, basis, skills, affordance, d, " Non-zero Displacements", tol));
  return t;
}
