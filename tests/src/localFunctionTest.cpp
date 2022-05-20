//
// Created by Alex on 21.04.2021.
//
#include <config.h>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "common.hh"
#include "testHelpers.hh"

#include <array>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <complex>
#include <vector>

#include <dune/common/classname.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <Eigen/Core>

#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/localFunctions/expressions.hh>
#include <ikarus/localFunctions/impl/projectionBasedLocalFunction.hh>
#include <ikarus/localFunctions/impl/standardLocalFunction.hh>
#include <ikarus/localFunctions/localFunctionName.hh>
#include <ikarus/manifolds/realTuple.hh>
#include <ikarus/manifolds/unitVector.hh>
#include <ikarus/utils/functionSanityChecks.hh>

using namespace Dune::Functions::BasisFactory;

template <typename T>
class LocalFunctionProjectionBasedUnitVector : public testing::Test {
public:
  static constexpr int size = T::value;
};
using test_types = ::testing::Types<std::integral_constant<std::size_t, 2>, std::integral_constant<std::size_t, 3>,
                                    std::integral_constant<std::size_t, 4>, std::integral_constant<std::size_t, 5>>;

TYPED_TEST_SUITE(LocalFunctionProjectionBasedUnitVector, test_types);

TYPED_TEST(LocalFunctionProjectionBasedUnitVector, ProjectionBasedUnitVector) {
  constexpr int size   = TypeParam::value;
  constexpr double tol = 1e-14;
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
    localBasis.bind(rule, bindDerivatives(0, 1));
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
    auto localFdual = [&](const auto& x, int gpI) {
      Dune::BlockVector<Ikarus::UnitVector<autodiff::dual, size>> v = vBlockedLocalDual;

      Ikarus::viewAsFlatEigenVector(v) += x;
      auto localF_ = Ikarus::ProjectionBasedLocalFunction(localBasis, v);
      return localF_.evaluateFunction(gpI);
    };

    const auto vasMat = Ikarus::viewAsEigenMatrixFixedDyn(vBlockedLocal);
    using namespace Ikarus::DerivativeDirections;
    for (int gpIndex = 0; auto& gp : rule) {
      const auto directorCached                          = localF.evaluateFunction(gpIndex);
      const Eigen::Vector<double, size> directorEmbedded = vasMat * localBasis.evaluateFunction(gpIndex);
      const auto directoreval                            = localF.evaluateFunction(gp.position());

      const auto J     = toEigenMatrix(ele.geometry().jacobianTransposed(gp.position())).transpose().eval();
      const auto Jinv  = J.inverse().eval();
      const auto jaco2 = localF.evaluateDerivative(gpIndex, wrt(spatialAll), transformWith(Jinv));

      const auto jaco2e     = localF.evaluateDerivative(gp.position(), wrt(spatialAll), transformWith(Jinv));
      const auto jaco2col0  = localF.evaluateDerivative(gpIndex, wrt(spatial(0)), transformWith(Jinv));
      const auto jaco2col0e = localF.evaluateDerivative(gp.position(), wrt(spatial(0)), transformWith(Jinv));
      const auto jaco2col1  = localF.evaluateDerivative(gpIndex, wrt(spatial(1)), transformWith(Jinv));
      const auto jaco2col1e = localF.evaluateDerivative(gp.position(), wrt(spatial(1)), transformWith(Jinv));
      EXPECT_THAT(jaco2.col(0), EigenApproxEqual(jaco2col0, tol));
      EXPECT_THAT(jaco2.col(1), EigenApproxEqual(jaco2col1, tol));
      EXPECT_THAT(jaco2, EigenApproxEqual(jaco2e, tol));
      EXPECT_THAT(jaco2col0, EigenApproxEqual(jaco2col0e, tol));
      EXPECT_THAT(jaco2col1, EigenApproxEqual(jaco2col1e, tol));

      // Check untransformed derivatives
      const auto jaco2un     = localF.evaluateDerivative(gpIndex, wrt(spatialAll));
      const auto jaco2col0un = localF.evaluateDerivative(gpIndex, wrt(spatial(0)));
      const auto jaco2col1un = localF.evaluateDerivative(gpIndex, wrt(spatial(1)));
      EXPECT_THAT(jaco2un.col(0), EigenApproxEqual(jaco2col0un, tol));
      EXPECT_THAT(jaco2un.col(1), EigenApproxEqual(jaco2col1un, tol));

      auto func = [&](auto& gpOffset_) { return localF.evaluateFunction(toFieldVector(gpOffset_)); };
      auto deriv
          = [&](auto& gpOffset_) { return localF.evaluateDerivative(toFieldVector(gpOffset_), wrt(spatialAll)); };
      Eigen::Vector<double, 2> gpOffset = toEigenVector(gp.position());
      auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(func, deriv), parameter(gpOffset));

      EXPECT_TRUE((checkJacobian<decltype(nonLinOp), Eigen::Vector<double, 2>>(nonLinOp, {.draw=false,.writeSlopeStatement=false})));

      auto localFdual_ = [&](auto& x) { return localFdual(x, gpIndex); };
      Eigen::VectorXdual xv(vasMat.cols() * vasMat.rows());
      xv.setZero();
      const Eigen::MatrixXd Jdual = jacobian(localFdual_, autodiff::wrt(xv), at(xv));

      const Eigen::Vector<double, size> testVec       = Eigen::Vector<double, size>::UnitX();
      const Eigen::Matrix<double, size, size> testMat = Eigen::Matrix<double, size, size>::Random();
      for (size_t i = 0; i < fe.size(); ++i) {
        const auto jacobianWRTCoeffs = localF.evaluateDerivative(gpIndex, wrt(coeff(i)));
        EXPECT_THAT(jacobianWRTCoeffs,
                    EigenApproxEqual(Jdual.block<size, size>(0, i * size) * vBlockedLocal[i].orthonormalFrame(), tol));

        const auto Warray  = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatialAll), transformWith(Jinv));
        const auto Warray2 = localF.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i)), transformWith(Jinv));
        for (int j = 0; j < 2; ++j)
          EXPECT_THAT(Warray[j], EigenApproxEqual(Warray2[j], tol));

        const auto W0 = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(0)), transformWith(Jinv));
        const auto W1 = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(1)), transformWith(Jinv));

        EXPECT_THAT(Warray[0], EigenApproxEqual(W0, tol));
        EXPECT_THAT(Warray[1], EigenApproxEqual(W1, tol));

        const auto Warrayun  = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatialAll));
        const auto Warray2un = localF.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i)));
        for (int j = 0; j < 2; ++j)
          EXPECT_THAT(Warrayun[j], EigenApproxEqual(Warray2un[j], tol));

        const auto W0un = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(0)));
        const auto W1un = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(1)));

        EXPECT_THAT(Warrayun[0], EigenApproxEqual(W0un, tol));
        EXPECT_THAT(Warrayun[1], EigenApproxEqual(W1un, tol));
        const auto Sun = localF.evaluateDerivative(gpIndex, wrt(coeff(i, i)), along(testVec));

        const auto S = localF.evaluateDerivative(gpIndex, wrt(coeff(i, i)), along(testVec), transformWith(Jinv));

        const auto chi
            = localF.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i, i)), along(testMat), transformWith(Jinv));

        const auto chi0 = localF.evaluateDerivative(gpIndex, wrt(spatial(0), coeff(i, i)), along(testMat.col(0)),
                                                    transformWith(Jinv));

        const auto chi1 = localF.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(i, i)), along(testMat.col(1)),
                                                    transformWith(Jinv));

        EXPECT_THAT(chi, EigenApproxEqual(chi0 + chi1, tol));
      }

      EXPECT_DOUBLE_EQ(directorCached.norm(), 1.0);
      EXPECT_DOUBLE_EQ(directoreval.norm(), 1.0);
      EXPECT_THAT(directorCached, EigenApproxEqual(directoreval, tol));
      EXPECT_NEAR((directoreval.transpose() * jaco2).norm(), 0.0, tol);
      EXPECT_NEAR(
          (Ikarus::UnitVector<double, size>::derivativeOfProjectionWRTposition(directoreval) * directoreval).norm(),
          0.0, tol);

      ++gpIndex;
    }
  }
}

template <typename T>
class LocalFunctionVector : public testing::Test {
public:
  using Manifold = T;
};
using test_types_1
    = ::testing::Types<Ikarus::RealTuple<double, 1>, Ikarus::RealTuple<double, 2>, Ikarus::RealTuple<double, 3>>;

TYPED_TEST_SUITE(LocalFunctionVector, test_types_1);

TYPED_TEST(LocalFunctionVector, Test1) {
  using Manifold     = typename TestFixture::Manifold;
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
    localBasis.bind(rule, bindDerivatives(0, 1));
    Dune::BlockVector<Manifold> vBlockedLocal(fe.size());
    Dune::BlockVector<typename Manifold::template Rebind<autodiff::dual>::other> vBlockedLocalDual(fe.size());
    Dune::BlockVector<typename Manifold::template Rebind<autodiff::dual2nd>::other> vBlockedLocalDual2nd(fe.size());
    for (size_t i = 0; i < fe.size(); ++i) {
      auto globalIndex = localView.index(localView.tree().child(0).localIndex(i));
      vBlockedLocal[i] = vBlocked[globalIndex[0]];
      vBlockedLocalDual[i].setValue(vBlocked[globalIndex[0]].getValue());
      vBlockedLocalDual2nd[i].setValue(vBlocked[globalIndex[0]].getValue());
    }
    auto localF     = Ikarus::StandardLocalFunction(localBasis, vBlockedLocal);
    auto localFdual = [&](const auto& x, int gpI) {
      Dune::BlockVector<typename Manifold::template Rebind<autodiff::dual>::other> v = vBlockedLocalDual;

      Ikarus::viewAsFlatEigenVector(v) += x;
      auto localF_ = Ikarus::StandardLocalFunction(localBasis, v);
      return localF_.evaluateFunction(gpI);
    };

    const auto vasMat = Ikarus::viewAsEigenMatrixFixedDyn(vBlockedLocal);
    using namespace Ikarus::DerivativeDirections;
    for (int gpIndex = 0; auto& gp : rule) {
      const auto& directorCached                         = localF.evaluateFunction(gpIndex);
      const Eigen::Vector<double, size> directorEmbedded = vasMat * localBasis.evaluateFunction(gpIndex);
      const auto& directoreval                           = localF.evaluateFunction(gp.position());

      const auto J     = toEigenMatrix(ele.geometry().jacobianTransposed(gp.position())).transpose().eval();
      const auto Jinv  = J.inverse().eval();
      const auto jaco2 = localF.evaluateDerivative(gpIndex, wrt(spatialAll), transformWith(Jinv));

      const auto jaco2e     = localF.evaluateDerivative(gp.position(), wrt(spatialAll), transformWith(Jinv));
      const auto jaco2col0  = localF.evaluateDerivative(gpIndex, wrt(spatial(0)), transformWith(Jinv));
      const auto jaco2col0e = localF.evaluateDerivative(gp.position(), wrt(spatial(0)), transformWith(Jinv));
      const auto jaco2col1  = localF.evaluateDerivative(gpIndex, wrt(spatial(1)), transformWith(Jinv));
      const auto jaco2col1e = localF.evaluateDerivative(gp.position(), wrt(spatial(1)), transformWith(Jinv));
      EXPECT_THAT(jaco2.col(0), EigenApproxEqual(jaco2col0, 1e-15));
      EXPECT_THAT(jaco2.col(1), EigenApproxEqual(jaco2col1, 1e-15));
      EXPECT_THAT(jaco2, EigenApproxEqual(jaco2e, 1e-15));
      EXPECT_THAT(jaco2col0, EigenApproxEqual(jaco2col0e, 1e-15));
      EXPECT_THAT(jaco2col1, EigenApproxEqual(jaco2col1e, 1e-15));

      // Check untransformed derivatives
      const auto jaco2un     = localF.evaluateDerivative(gpIndex, wrt(spatialAll));
      const auto jaco2col0un = localF.evaluateDerivative(gpIndex, wrt(spatial(0)));
      const auto jaco2col1un = localF.evaluateDerivative(gpIndex, wrt(spatial(1)));
      EXPECT_THAT(jaco2un.col(0), EigenApproxEqual(jaco2col0un, 1e-15));
      EXPECT_THAT(jaco2un.col(1), EigenApproxEqual(jaco2col1un, 1e-15));

      auto func = [&](auto& gpOffset_) { return localF.evaluateFunction(toFieldVector(gpOffset_)); };
      auto deriv
          = [&](auto& gpOffset_) { return localF.evaluateDerivative(toFieldVector(gpOffset_), wrt(spatialAll)); };
      Eigen::Vector<double, 2> gpOffset = toEigenVector(gp.position());
      auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(func, deriv), parameter(gpOffset));

      EXPECT_TRUE((checkJacobian<decltype(nonLinOp), Eigen::Vector<double, 2>>(nonLinOp, {.draw=false,.writeSlopeStatement=false})));

      auto localFdual_ = [&](auto& x) { return localFdual(x, gpIndex); };
      Eigen::VectorXdual xv(vasMat.cols() * vasMat.rows());
      xv.setZero();
      const Eigen::MatrixXd Jdual = jacobian(localFdual_, autodiff::wrt(xv), at(xv));

      const Eigen::Vector<double, size> testVec = Eigen::Vector<double, size>::UnitX();
      for (size_t i = 0; i < fe.size(); ++i) {
        const auto jacobianWRTCoeffs = localF.evaluateDerivative(gpIndex, wrt(coeff(i)));
        EXPECT_THAT(jacobianWRTCoeffs, EigenApproxEqual(Jdual.block<size, size>(0, i * size), 1e-15));

        const auto Warray  = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatialAll), transformWith(Jinv));
        const auto Warray2 = localF.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i)), transformWith(Jinv));
        for (int j = 0; j < 2; ++j)
          EXPECT_THAT(Warray[j], EigenApproxEqual(Warray2[j], 1e-15));

        const auto W0 = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(0)), transformWith(Jinv));
        const auto W1 = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(1)), transformWith(Jinv));

        EXPECT_THAT(Warray[0], EigenApproxEqual(W0, 1e-15));
        EXPECT_THAT(Warray[1], EigenApproxEqual(W1, 1e-15));

        const auto Warrayun  = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatialAll));
        const auto Warray2un = localF.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i)));
        for (int j = 0; j < 2; ++j)
          EXPECT_THAT(Warrayun[j], EigenApproxEqual(Warray2un[j], 1e-15));

        const auto W0un = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(0)));
        const auto W1un = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(1)));

        EXPECT_THAT(Warrayun[0], EigenApproxEqual(W0un, 1e-15));
        EXPECT_THAT(Warrayun[1], EigenApproxEqual(W1un, 1e-15));
      }

      EXPECT_THAT(directorCached, EigenApproxEqual(directoreval, 1e-15));

      ++gpIndex;
    }
  }
}
