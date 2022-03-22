//
// Created by Alex on 21.04.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "../../config.h"
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

#include "ikarus/LinearAlgebra/NonLinearOperator.h"
#include <ikarus/LocalFunctions/ProjectionBasedLocalFunction.h>
#include <ikarus/LocalFunctions/StandardLocalFunction.h>
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/utils/functionSanityChecks.h>
#include "common.h"
const double tol = 1e-15;
using namespace Dune::Functions::BasisFactory;




template <typename T>
class LocalFunctionProjectionBasedUnitVector : public testing::Test {
public:
  static constexpr int size = T::value;

};
using test_types = ::testing::Types<
                                                        std::integral_constant<std::size_t,2>,
                                                        std::integral_constant<std::size_t,3>,
                                                        std::integral_constant<std::size_t,4>,
    std::integral_constant<std::size_t,5>>;


TYPED_TEST_SUITE(LocalFunctionProjectionBasedUnitVector, test_types);

TYPED_TEST(LocalFunctionProjectionBasedUnitVector, ProjectionBasedUnitVector) {
  constexpr int size = TypeParam::value;
  using namespace Ikarus;
  auto grid = createGrid<Grids::Alu>();

  auto gridView = grid->leafGridView();
  auto basis    = makeBasis(gridView, power<3>(lagrange<1>(), BlockedInterleaved()));
  Dune::BlockVector<Ikarus::UnitVector<double, size>> vBlocked(basis.size());
  for (auto& vsingle : vBlocked) {
    vsingle.setValue(0.1 * Eigen::Vector<double, size>::Random() + Eigen::Vector<double, size>::UnitX());
  }
  auto localView = basis.localView();
  for (auto& ele : elements(gridView)) {
    localView.bind(ele);
    const auto& fe   = localView.tree().child(0).finiteElement();
    auto localBasis  = Ikarus::LocalBasis(localView.tree().child(0).finiteElement().localBasis());
    const auto& rule = Dune::QuadratureRules<double, 2>::rule(localView.element().type(), 3);
    localBasis.bind(rule, 0, 1);
    Dune::BlockVector<Ikarus::UnitVector<double, size>> vBlockedLocal(fe.size());
    Dune::BlockVector<Ikarus::UnitVector<autodiff::dual, size>> vBlockedLocalDual(fe.size());
    Dune::BlockVector<Ikarus::UnitVector<autodiff::dual2nd, size>> vBlockedLocalDual2nd(fe.size());
    for (size_t i = 0; i < fe.size(); ++i) {
      auto globalIndex = localView.index(localView.tree().child(0).localIndex(i));
      vBlockedLocal[i] = vBlocked[globalIndex[0]];
      vBlockedLocalDual[i].setValue(vBlocked[globalIndex[0]].getValue());
      vBlockedLocalDual2nd[i].setValue(vBlocked[globalIndex[0]].getValue());
    }
    auto localF     = Ikarus::ProjectionBasedLocalFunction(localBasis, vBlockedLocal);
    auto localFdual = [&](auto& x, int gpI) {
      Dune::BlockVector<Ikarus::UnitVector<autodiff::dual, size>> v = vBlockedLocalDual;

      Ikarus::LinearAlgebra::viewAsFlatEigenVector(v) += x;
      auto localF_ = Ikarus::ProjectionBasedLocalFunction(localBasis, v);
      return localF_.evaluateFunction(gpI).getValue();
    };

    const auto vasMat = Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(vBlockedLocal);
    using namespace Ikarus::DerivativeDirections;
    for (int gpIndex = 0; auto& gp : rule) {
      const auto& directorCached             = localF.evaluateFunction(gpIndex).getValue();
      const Eigen::Vector<double,size> directorEmbedded = vasMat * localBasis.getFunction(gpIndex);
      const auto& directoreval               = localF.evaluateFunction(gp.position());

      const auto J    = toEigenMatrix(ele.geometry().jacobianTransposed(gp.position())).transpose().eval();
      const auto Jinv = J.inverse().eval();
      const auto jaco2 = localF.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv));

      const auto jaco2e     = localF.evaluateDerivative(gp.position(), wrt(spatialall), transformWith(Jinv));
      const auto jaco2col0  = localF.evaluateDerivative(gpIndex, wrt(spatial<0>{}), transformWith(Jinv));
      const auto jaco2col0e = localF.evaluateDerivative(gp.position(), wrt(spatial<0>{}), transformWith(Jinv));
      const auto jaco2col1  = localF.evaluateDerivative(gpIndex, wrt(spatial1), transformWith(Jinv));
      const auto jaco2col1e = localF.evaluateDerivative(gp.position(), wrt(spatial1), transformWith(Jinv));
      EXPECT_THAT(jaco2.col(0), EigenApproxEqual(jaco2col0, 1e-15));
      EXPECT_THAT(jaco2.col(1), EigenApproxEqual(jaco2col1, 1e-15));
      EXPECT_THAT(jaco2, EigenApproxEqual(jaco2e, 1e-15));
      EXPECT_THAT(jaco2col0, EigenApproxEqual(jaco2col0e, 1e-15));
      EXPECT_THAT(jaco2col1, EigenApproxEqual(jaco2col1e, 1e-15));

      // Check untransformed derivatives
      const auto jaco2un     = localF.evaluateDerivative(gpIndex, wrt(spatialall));
      const auto jaco2col0un = localF.evaluateDerivative(gpIndex, wrt(spatial<0>{}));
      const auto jaco2col1un = localF.evaluateDerivative(gpIndex, wrt(spatial1));
      EXPECT_THAT(jaco2un.col(0), EigenApproxEqual(jaco2col0un, 1e-15));
      EXPECT_THAT(jaco2un.col(1), EigenApproxEqual(jaco2col1un, 1e-15));

      auto func = [&](auto& gpOffset_) { return localF.evaluateFunction(toFieldVector(gpOffset_)).getValue(); };
      auto deriv
          = [&](auto& gpOffset_) { return localF.evaluateDerivative(toFieldVector(gpOffset_), wrt(spatialall)); };
      Eigen::Vector<double, 2> gpOffset = toEigenVector(gp.position());
      auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(func, deriv), parameter(gpOffset));

      EXPECT_TRUE((checkJacobian<decltype(nonLinOp), Eigen::Vector<double, 2>>(nonLinOp, false)));

      auto localFdual_ = [&](auto& x) { return localFdual(x, gpIndex); };
      Eigen::VectorXdual xv(vasMat.cols() * vasMat.rows());
      xv.setZero();
      const Eigen::MatrixXd Jdual = jacobian(localFdual_, autodiff::wrt(xv), at(xv));

      const Eigen::Vector<double,size> testVec = Eigen::Vector<double,size> ::UnitX();
      for (size_t i = 0; i < fe.size(); ++i) {
        const auto jacobianWRTCoeffs = localF.evaluateDerivative(gpIndex, wrt(coeffs), coeffIndices(i));
        EXPECT_THAT(jacobianWRTCoeffs, EigenApproxEqual(Jdual.block<size, size>(0, i * size), 1e-15));

        const auto Warray
            = localF.evaluateDerivative(gpIndex, wrt(coeffs, spatialall), transformWith(Jinv), coeffIndices(i));
        const auto Warray2
            = localF.evaluateDerivative(gpIndex, wrt(spatialall, coeffs), transformWith(Jinv), coeffIndices(i));
        for (int j = 0; j < 2; ++j)
          EXPECT_THAT(Warray[j], EigenApproxEqual(Warray2[j], 1e-15));

        const auto W0 = localF.evaluateDerivative(gpIndex, wrt(coeffs, spatial0), transformWith(Jinv), coeffIndices(i));
        const auto W1 = localF.evaluateDerivative(gpIndex, wrt(coeffs, spatial1), transformWith(Jinv), coeffIndices(i));

        EXPECT_THAT(Warray[0], EigenApproxEqual(W0, 1e-15));
        EXPECT_THAT(Warray[1], EigenApproxEqual(W1, 1e-15));

        const auto Warrayun  = localF.evaluateDerivative(gpIndex, wrt(coeffs, spatialall), coeffIndices(i));
        const auto Warray2un = localF.evaluateDerivative(gpIndex, wrt(spatialall, coeffs), coeffIndices(i));
        for (int j = 0; j < 2; ++j)
          EXPECT_THAT(Warrayun[j], EigenApproxEqual(Warray2un[j], 1e-15));

        const auto W0un = localF.evaluateDerivative(gpIndex, wrt(coeffs, spatial0), coeffIndices(i));
        const auto W1un = localF.evaluateDerivative(gpIndex, wrt(coeffs, spatial1), coeffIndices(i));

        EXPECT_THAT(Warrayun[0], EigenApproxEqual(W0un, 1e-15));
        EXPECT_THAT(Warrayun[1], EigenApproxEqual(W1un, 1e-15));
        const auto Sun = localF.evaluateDerivative(gpIndex, wrt(coeffs, coeffs), along(testVec), coeffIndices(i, i));

        const auto S = localF.evaluateDerivative(gpIndex, wrt(coeffs, coeffs), along(testVec), transformWith(Jinv),
                                                 coeffIndices(i, i));

        const auto chi = localF.evaluateDerivative(gpIndex, wrt(spatialall, coeffs, coeffs), along(testVec),
                                                   transformWith(Jinv), coeffIndices(i, i));

        const auto chi0 = localF.evaluateDerivative(gpIndex, wrt(spatial0, coeffs, coeffs), along(testVec),
                                                    transformWith(Jinv), coeffIndices(i, i));

        const auto chi1 = localF.evaluateDerivative(gpIndex, wrt(spatial1, coeffs, coeffs), along(testVec),
                                                    transformWith(Jinv), coeffIndices(i, i));

        EXPECT_THAT(chi[0], EigenApproxEqual(chi0, 1e-15));
        EXPECT_THAT(chi[1], EigenApproxEqual(chi1, 1e-15));
      }

      EXPECT_DOUBLE_EQ(directorCached.norm(), 1.0);
      EXPECT_DOUBLE_EQ(directoreval.getValue().norm(), 1.0);
      EXPECT_THAT(directorCached, EigenApproxEqual(directoreval.getValue(), 1e-15));
      EXPECT_NEAR((directoreval.getValue().transpose() * jaco2).norm(), 0.0, 1e-15);
      EXPECT_NEAR((Ikarus::UnitVector<double, size>::derivativeOfProjectionWRTposition(directoreval.getValue())
                   * directoreval.getValue())
                      .norm(),
                  0.0, 1e-15);

      ++gpIndex;
    }
  }
}


template <typename T>
class LocalFunctionVector : public testing::Test {
public:
  using Manifold = T;

};
using test_types_1 = ::testing::Types<
    Ikarus::RealTuple<double, 1>,
    Ikarus::RealTuple<double, 2>,
    Ikarus::RealTuple<double, 3>>;


TYPED_TEST_SUITE(LocalFunctionVector, test_types_1);

TYPED_TEST(LocalFunctionVector,Test1) {
  using Manifold = typename TestFixture::Manifold;
  constexpr int size = Manifold::valueSize;
  using namespace Ikarus;
  auto grid = createGrid<Grids::Yasp>();

  auto gridView = grid->leafGridView();
  auto basis    = makeBasis(gridView, power<3>(lagrange<1>(), BlockedInterleaved()));
  Dune::BlockVector<Manifold> vBlocked(basis.size());
  for (auto& vsingle : vBlocked) {
    vsingle.setValue(0.1 * Eigen::Vector<double, size>::Random() + Eigen::Vector<double, size>::UnitX());
  }
  auto localView = basis.localView();
  for (auto& ele : elements(gridView)) {
    localView.bind(ele);
    const auto& fe   = localView.tree().child(0).finiteElement();
    auto localBasis  = Ikarus::LocalBasis(localView.tree().child(0).finiteElement().localBasis());
    const auto& rule = Dune::QuadratureRules<double, 2>::rule(localView.element().type(), 3);
    localBasis.bind(rule, 0, 1);
    Dune::BlockVector<Manifold> vBlockedLocal(fe.size());
    Dune::BlockVector<typename Manifold::template Rebind<autodiff::dual>::type> vBlockedLocalDual(fe.size());
    Dune::BlockVector<typename Manifold::template Rebind<autodiff::dual2nd>::type> vBlockedLocalDual2nd(fe.size());
    for (size_t i = 0; i < fe.size(); ++i) {
      auto globalIndex = localView.index(localView.tree().child(0).localIndex(i));
      vBlockedLocal[i] = vBlocked[globalIndex[0]];
      vBlockedLocalDual[i].setValue(vBlocked[globalIndex[0]].getValue());
      vBlockedLocalDual2nd[i].setValue(vBlocked[globalIndex[0]].getValue());
    }
    auto localF     = Ikarus::StandardLocalFunction(localBasis, vBlockedLocal);
    auto localFdual = [&](auto& x, int gpI) {
      Dune::BlockVector<typename Manifold::template Rebind<autodiff::dual>::type> v = vBlockedLocalDual;

      Ikarus::LinearAlgebra::viewAsFlatEigenVector(v) += x;
      auto localF_ = Ikarus::StandardLocalFunction(localBasis, v);
      return localF_.evaluateFunction(gpI).getValue();
    };

    const auto vasMat = Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(vBlockedLocal);
    using namespace Ikarus::DerivativeDirections;
    for (int gpIndex = 0; auto& gp : rule) {
      const auto& directorCached             = localF.evaluateFunction(gpIndex).getValue();
      const Eigen::Vector<double,size> directorEmbedded = vasMat * localBasis.getFunction(gpIndex);
      const auto& directoreval               = localF.evaluateFunction(gp.position());

      const auto J    = toEigenMatrix(ele.geometry().jacobianTransposed(gp.position())).transpose().eval();
      const auto Jinv = J.inverse().eval();
      const auto jaco2 = localF.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv));

      const auto jaco2e     = localF.evaluateDerivative(gp.position(), wrt(spatialall), transformWith(Jinv));
      const auto jaco2col0  = localF.evaluateDerivative(gpIndex, wrt(spatial<0>{}), transformWith(Jinv));
      const auto jaco2col0e = localF.evaluateDerivative(gp.position(), wrt(spatial<0>{}), transformWith(Jinv));
      const auto jaco2col1  = localF.evaluateDerivative(gpIndex, wrt(spatial1), transformWith(Jinv));
      const auto jaco2col1e = localF.evaluateDerivative(gp.position(), wrt(spatial1), transformWith(Jinv));
      EXPECT_THAT(jaco2.col(0), EigenApproxEqual(jaco2col0, 1e-15));
      EXPECT_THAT(jaco2.col(1), EigenApproxEqual(jaco2col1, 1e-15));
      EXPECT_THAT(jaco2, EigenApproxEqual(jaco2e, 1e-15));
      EXPECT_THAT(jaco2col0, EigenApproxEqual(jaco2col0e, 1e-15));
      EXPECT_THAT(jaco2col1, EigenApproxEqual(jaco2col1e, 1e-15));

      // Check untransformed derivatives
      const auto jaco2un     = localF.evaluateDerivative(gpIndex, wrt(spatialall));
      const auto jaco2col0un = localF.evaluateDerivative(gpIndex, wrt(spatial<0>{}));
      const auto jaco2col1un = localF.evaluateDerivative(gpIndex, wrt(spatial1));
      EXPECT_THAT(jaco2un.col(0), EigenApproxEqual(jaco2col0un, 1e-15));
      EXPECT_THAT(jaco2un.col(1), EigenApproxEqual(jaco2col1un, 1e-15));

      auto func = [&](auto& gpOffset_) { return localF.evaluateFunction(toFieldVector(gpOffset_)).getValue(); };
      auto deriv
          = [&](auto& gpOffset_) { return localF.evaluateDerivative(toFieldVector(gpOffset_), wrt(spatialall)); };
      Eigen::Vector<double, 2> gpOffset = toEigenVector(gp.position());
      auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(func, deriv), parameter(gpOffset));

      EXPECT_TRUE((checkJacobian<decltype(nonLinOp), Eigen::Vector<double, 2>>(nonLinOp, false)));

      auto localFdual_ = [&](auto& x) { return localFdual(x, gpIndex); };
      Eigen::VectorXdual xv(vasMat.cols() * vasMat.rows());
      xv.setZero();
      const Eigen::MatrixXd Jdual = jacobian(localFdual_, autodiff::wrt(xv), at(xv));

      const Eigen::Vector<double,size> testVec = Eigen::Vector<double,size> ::UnitX();
      for (size_t i = 0; i < fe.size(); ++i) {
        const auto jacobianWRTCoeffs = localF.evaluateDerivative(gpIndex, wrt(coeffs), coeffIndices(i));
        EXPECT_THAT(jacobianWRTCoeffs, EigenApproxEqual(Jdual.block<size, size>(0, i * size), 1e-15));

        const auto Warray
            = localF.evaluateDerivative(gpIndex, wrt(coeffs, spatialall), transformWith(Jinv), coeffIndices(i));
        const auto Warray2
            = localF.evaluateDerivative(gpIndex, wrt(spatialall, coeffs), transformWith(Jinv), coeffIndices(i));
        for (int j = 0; j < 2; ++j)
          EXPECT_THAT(Warray[j], EigenApproxEqual(Warray2[j], 1e-15));

        const auto W0 = localF.evaluateDerivative(gpIndex, wrt(coeffs, spatial0), transformWith(Jinv), coeffIndices(i));
        const auto W1 = localF.evaluateDerivative(gpIndex, wrt(coeffs, spatial1), transformWith(Jinv), coeffIndices(i));

        EXPECT_THAT(Warray[0], EigenApproxEqual(W0, 1e-15));
        EXPECT_THAT(Warray[1], EigenApproxEqual(W1, 1e-15));

        const auto Warrayun  = localF.evaluateDerivative(gpIndex, wrt(coeffs, spatialall), coeffIndices(i));
        const auto Warray2un = localF.evaluateDerivative(gpIndex, wrt(spatialall, coeffs), coeffIndices(i));
        for (int j = 0; j < 2; ++j)
          EXPECT_THAT(Warrayun[j], EigenApproxEqual(Warray2un[j], 1e-15));

        const auto W0un = localF.evaluateDerivative(gpIndex, wrt(coeffs, spatial0), coeffIndices(i));
        const auto W1un = localF.evaluateDerivative(gpIndex, wrt(coeffs, spatial1), coeffIndices(i));

        EXPECT_THAT(Warrayun[0], EigenApproxEqual(W0un, 1e-15));
        EXPECT_THAT(Warrayun[1], EigenApproxEqual(W1un, 1e-15));
      }


      EXPECT_THAT(directorCached, EigenApproxEqual(directoreval.getValue(), 1e-15));

      ++gpIndex;
    }
  }
}