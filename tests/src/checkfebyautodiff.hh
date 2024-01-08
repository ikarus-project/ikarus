// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later


#pragma once
#include <dune/common/test/testsuite.hh>

template <template<typename ,typename > class FE, typename GridView,typename PreBasis, typename... ElementArgsType>
auto checkFEByAutoDiff( const GridView& gridView, const PreBasis& pb,const ElementArgsType&... eleArgs) {

  using namespace Dune::Functions::BasisFactory;
  auto basis       = Ikarus::makeBasis(gridView, pb);
  auto element     = gridView.template begin<0>();
  auto nDOF        = basis.flat().size();
  const double tol = 1e-10;

  FE fe(basis,*element,eleArgs...);

  const std::string feClassName = Dune::className(fe);
  Dune::TestSuite t("Check calculateScalarImpl() and calculateVectorImpl() by Automatic Differentiation for gridDim = "
            + feClassName);
  using AutoDiffBasedFE = Ikarus::AutoDiffFE<decltype(fe)>;
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
          "automatic differentiation for "
              + feClassName);

  t.check(R.isApprox(RAutoDiff, tol),
          "Mismatch between the residual vectors obtained from explicit implementation and the one based on "
          "automatic differentiation for "
              + feClassName);

  t.check(Dune::FloatCmp::eq(fe.calculateScalar(req), feAutoDiff.calculateScalar(req), tol),
          "Mismatch between the energies obtained from explicit implementation and the one based on "
          "automatic differentiation for "
              + feClassName);

  return t;
}
