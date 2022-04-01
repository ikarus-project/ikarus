//
// Created by Alex on 21.04.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../../config.h"
#include "common.h"
#include "testHelpers.h"

#include <array>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <complex>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <Eigen/Core>

#include <ikarus/LinearAlgebra/NonLinearOperator.h>
#include <ikarus/LocalFunctions/ProjectionBasedLocalFunction.h>
#include <ikarus/LocalFunctions/StandardLocalFunction.h>
#include <ikarus/LocalFunctions/UnitNormalLocalFunction.h>
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/utils/functionSanityChecks.h>

using namespace Dune::Functions::BasisFactory;

TEST(UnitNormalLocalFunction, Test1) {
  using Manifold     = Ikarus::RealTuple<double, 3>;
  constexpr int size = Manifold::valueSize;
  using namespace Ikarus;
  using Grid        = Dune::YaspGrid<2>;
  const double L    = 2;
  const double h    = 1;
  const size_t elex = 1;
  const size_t eley = 1;

  Dune::FieldVector<double, 2> bbox = {L, h};
  std::array<int, 2> eles           = {elex, eley};
  auto grid                         = std::make_shared<Grid>(bbox, eles);

  auto gridView = grid->leafGridView();
  auto basis    = makeBasis(gridView, power<3>(lagrange<1>(), BlockedInterleaved()));

  Dune::BlockVector<Manifold> vBlocked(basis.size());
  for (auto& vsingle : vBlocked)
    vsingle.setValue(0.1 * Eigen::Vector<double, size>::Random() + Eigen::Vector<double, size>::UnitX());

  auto localView = basis.localView();
  for (auto& ele : elements(gridView)) {
    localView.bind(ele);
    const auto& fe  = localView.tree().child(0).finiteElement();
    auto localBasis = Ikarus::LocalBasis(localView.tree().child(0).finiteElement().localBasis());

    const auto& rule = Dune::QuadratureRules<double, 2>::rule(localView.element().type(), 3);
    localBasis.bind(rule, bindDerivatives(0, 1));

    Dune::BlockVector<Ikarus::RealTuple<double, 3>> vBlockedLocal(fe.size());

    for (size_t i = 0; i < fe.size(); ++i) {
      auto globalIndex = localView.index(localView.tree().child(0).localIndex(i));
      vBlockedLocal[i] = vBlocked[globalIndex[0]];
    }

    auto localF = Ikarus::UnitNormalLocalFunction(localBasis, vBlockedLocal);
    using namespace Ikarus::DerivativeDirections;
    for (int gpIndex = 0; auto& gp : rule) {
      for (size_t i = 0; i < fe.size(); ++i) {
        const auto N         = localF.evaluateFunction(gpIndex).getValue();
        const auto dNdC         = localF.evaluateDerivative(gpIndex, wrt(coeffs), coeffIndices(i));
        const auto NExptected = Eigen::Vector<double,size>::UnitZ();
        const auto dNdCExpected = Eigen::Matrix<double,size,size>::Identity();
        // Check untransformed derivatives
        EXPECT_THAT(N, EigenApproxEqual(NExptected, 1e-15));
        EXPECT_DOUBLE_EQ(N.norm(),1);
        EXPECT_THAT(dNdC, EigenApproxEqual(dNdCExpected, 1e-15));
      }
      ++gpIndex;
    }
  }
}
