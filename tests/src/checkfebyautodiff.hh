// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <dune/common/float_cmp.hh>
#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/finiteelements/autodifffe.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/utils/basis.hh>

template <typename GridView, typename PreBasis, typename Skills>
auto checkFESByAutoDiff(const GridView& gridView, const PreBasis& pb, Skills&& skills) {
  auto basis = Ikarus::makeBasis(gridView, pb);
  Eigen::VectorXd d;
  d.setRandom(basis.flat().dimension());
  double lambda = 7.3;

  auto req = Ikarus::FErequirements().addAffordance(Ikarus::AffordanceCollections::elastoStatics);
  req.insertGlobalSolution(Ikarus::FESolutions::displacement, d)
      .insertParameter(Ikarus::FEParameter::loadfactor, lambda);
  Dune::TestSuite t("Check calculateScalarImpl() and calculateVectorImpl() by Automatic Differentiation");
  for (auto element : elements(gridView)) {
    auto localView = basis.flat().localView();
    localView.bind(element);
    auto nDOF        = localView.size();
    const double tol = 1e-10;

    auto fe = Ikarus::makeFE(basis, std::forward<Skills>(skills));
    fe.bind(element);

    const std::string feClassName = Dune::className(fe);

    using AutoDiffBasedFE = Ikarus::AutoDiffFE<decltype(fe), true>;
    AutoDiffBasedFE feAutoDiff(fe);

    Eigen::MatrixXd K, KAutoDiff;
    K.setZero(nDOF, nDOF);
    KAutoDiff.setZero(nDOF, nDOF);

    Eigen::VectorXd R, RAutoDiff;
    R.setZero(nDOF);
    RAutoDiff.setZero(nDOF);

    calculateMatrix(fe, req, K);
    calculateMatrix(feAutoDiff, req, KAutoDiff);

    calculateVector(fe, req, R);
    calculateVector(feAutoDiff, req, RAutoDiff);

    t.check(K.isApprox(KAutoDiff, tol),
            "Mismatch between the stiffness matrices obtained from explicit implementation and the one based on "
            "automatic differentiation for " +
                feClassName)
        << "The difference is " << (K - KAutoDiff);

    t.check(R.isApprox(RAutoDiff, tol),
            "Mismatch between the residual vectors obtained from explicit implementation and the one based on "
            "automatic differentiation for " +
                feClassName)
        << "The difference is " << (R - RAutoDiff);

    t.check(Dune::FloatCmp::eq(calculateScalar(fe, req), calculateScalar(feAutoDiff, req), tol),
            "Mismatch between the energies obtained from explicit implementation and the one based on "
            "automatic differentiation for " +
                feClassName);
  }

  return t;
}
