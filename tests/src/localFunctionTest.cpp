//
#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
using Dune::TestSuite;

#include "common.hh"
#include "factories.hh"
#include "testHelpers.hh"

#include <array>
#include <complex>
#include <vector>

#include <dune/common/classname.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <Eigen/Core>

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/localFunctions/expressions.hh>
#include <ikarus/localFunctions/impl/projectionBasedLocalFunction.hh>
#include <ikarus/localFunctions/impl/standardLocalFunction.hh>
#include <ikarus/localFunctions/localFunctionName.hh>
#include <ikarus/manifolds/realTuple.hh>
#include <ikarus/manifolds/unitVector.hh>
#include <ikarus/utils/functionSanityChecks.hh>
using namespace Dune::Functions::BasisFactory;

template <int Dim>
auto projectionBasedLocalFunctionTest() {
  TestSuite t("projectionBasedLocalFunctionTest" + std::to_string(Dim));
  constexpr int size   = Dim;
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
      const auto directorCached = localF.evaluateFunction(gpIndex);
      const auto directoreval   = localF.evaluateFunction(gp.position());

      const auto J     = toEigen(ele.geometry().jacobianTransposed(gp.position())).transpose().eval();
      const auto Jinv  = J.inverse().eval();
      const auto jaco2 = localF.evaluateDerivative(gpIndex, wrt(spatialAll), transformWith(Jinv));

      const auto jaco2e     = localF.evaluateDerivative(gp.position(), wrt(spatialAll), transformWith(Jinv));
      const auto jaco2col0  = localF.evaluateDerivative(gpIndex, wrt(spatial(0)), transformWith(Jinv));
      const auto jaco2col0e = localF.evaluateDerivative(gp.position(), wrt(spatial(0)), transformWith(Jinv));
      const auto jaco2col1  = localF.evaluateDerivative(gpIndex, wrt(spatial(1)), transformWith(Jinv));
      const auto jaco2col1e = localF.evaluateDerivative(gp.position(), wrt(spatial(1)), transformWith(Jinv));
      t.check(isApproxSame(jaco2.col(0), jaco2col0, tol));
      t.check(isApproxSame(jaco2.col(1), jaco2col1, tol));
      t.check(isApproxSame(jaco2, jaco2e, tol));
      t.check(isApproxSame(jaco2col0, jaco2col0e, tol));
      t.check(isApproxSame(jaco2col1, jaco2col1e, tol));

      // Check untransformed derivatives
      const auto jaco2un     = localF.evaluateDerivative(gpIndex, wrt(spatialAll));
      const auto jaco2col0un = localF.evaluateDerivative(gpIndex, wrt(spatial(0)));
      const auto jaco2col1un = localF.evaluateDerivative(gpIndex, wrt(spatial(1)));
      t.check(isApproxSame(jaco2un.col(0), jaco2col0un, tol), "spatialAll[0] == spatial(0)");
      t.check(isApproxSame(jaco2un.col(1), jaco2col1un, tol), "spatialAll[1] == spatial(1)");

      auto func  = [&](auto& gpOffset_) { return localF.evaluateFunction(toDune(gpOffset_)); };
      auto deriv = [&](auto& gpOffset_) { return localF.evaluateDerivative(toDune(gpOffset_), wrt(spatialAll)); };
      Eigen::Vector<double, 2> gpOffset = toEigen(gp.position());
      auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(func, deriv), parameter(gpOffset));

      t.check((checkJacobian<decltype(nonLinOp), Eigen::Vector<double, 2>>(
          nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = false})));

      auto localFdual_ = [&](auto& x) { return localFdual(x, gpIndex); };
      Eigen::VectorXdual xv(vasMat.cols() * vasMat.rows());
      xv.setZero();
      const Eigen::MatrixXd Jdual = jacobian(localFdual_, autodiff::wrt(xv), at(xv));

      const Eigen::Matrix<double, size, size> testMat = Eigen::Matrix<double, size, size>::Random();
      for (size_t i = 0; i < fe.size(); ++i) {
        const auto jacobianWRTCoeffs = localF.evaluateDerivative(gpIndex, wrt(coeff(i)));
        t.check(isApproxSame(jacobianWRTCoeffs,
                             Jdual.block<size, size>(0, i * size) * vBlockedLocal[i].orthonormalFrame(), tol));

        const auto Warray  = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatialAll), transformWith(Jinv));
        const auto Warray2 = localF.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i)), transformWith(Jinv));
        for (int j = 0; j < 2; ++j)
          t.check(isApproxSame(Warray[j], Warray2[j], tol));

        const auto W0 = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(0)), transformWith(Jinv));
        const auto W1 = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(1)), transformWith(Jinv));

        t.check(isApproxSame(Warray[0], W0, tol));
        t.check(isApproxSame(Warray[1], W1, tol));

        const auto Warrayun  = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatialAll));
        const auto Warray2un = localF.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i)));
        for (int j = 0; j < 2; ++j)
          t.check(isApproxSame(Warrayun[j], Warray2un[j], tol));

        const auto W0un = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(0)));
        const auto W1un = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(1)));

        t.check(isApproxSame(Warrayun[0], W0un, tol));
        t.check(isApproxSame(Warrayun[1], W1un, tol));
        //        const auto Sun = localF.evaluateDerivative(gpIndex, wrt(coeff(i, i)), along(testVec));

        //        const auto S = localF.evaluateDerivative(gpIndex, wrt(coeff(i, i)), along(testVec),
        //        transformWith(Jinv));

        const auto chi
            = localF.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i, i)), along(testMat), transformWith(Jinv));

        const auto chi0 = localF.evaluateDerivative(gpIndex, wrt(spatial(0), coeff(i, i)), along(testMat.col(0)),
                                                    transformWith(Jinv));

        const auto chi1 = localF.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(i, i)), along(testMat.col(1)),
                                                    transformWith(Jinv));

        t.check(isApproxSame(chi, chi0 + chi1, tol));
      }

      t.check(Dune::FloatCmp::eq(1.0, directorCached.norm(), tol), "DirectorLength");
      t.check(Dune::FloatCmp::eq(1.0, directoreval.norm(), tol), "DirectorLength2");
      t.check(isApproxSame(directorCached, directoreval, tol), "directorCached==directoreval");
      t.check(std::abs((directoreval.transpose() * jaco2).norm()) < tol, "Director is normal to gradDirector ");
      t.check(

          std::abs(
              (Ikarus::UnitVector<double, size>::derivativeOfProjectionWRTposition(directoreval) * directoreval).norm())
              < tol,
          "Director is in kernel of tangent base projector ");

      ++gpIndex;
    }
  }
  return t;
}

template <int Dim>
auto standardLocalFunctionTest() {
  TestSuite t("standardLocalFunctionTest" + std::to_string(Dim));
  using Manifold     = Ikarus::RealTuple<double, Dim>;
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
      const auto& directorCached = localF.evaluateFunction(gpIndex);
      const auto& directoreval   = localF.evaluateFunction(gp.position());

      const auto J     = toEigen(ele.geometry().jacobianTransposed(gp.position())).transpose().eval();
      const auto Jinv  = J.inverse().eval();
      const auto jaco2 = localF.evaluateDerivative(gpIndex, wrt(spatialAll), transformWith(Jinv));

      const auto jaco2e     = localF.evaluateDerivative(gp.position(), wrt(spatialAll), transformWith(Jinv));
      const auto jaco2col0  = localF.evaluateDerivative(gpIndex, wrt(spatial(0)), transformWith(Jinv));
      const auto jaco2col0e = localF.evaluateDerivative(gp.position(), wrt(spatial(0)), transformWith(Jinv));
      const auto jaco2col1  = localF.evaluateDerivative(gpIndex, wrt(spatial(1)), transformWith(Jinv));
      const auto jaco2col1e = localF.evaluateDerivative(gp.position(), wrt(spatial(1)), transformWith(Jinv));
      t.check(isApproxSame(jaco2.col(0), jaco2col0, 1e-15));
      t.check(isApproxSame(jaco2.col(1), jaco2col1, 1e-15));
      t.check(isApproxSame(jaco2, jaco2e, 1e-15));
      t.check(isApproxSame(jaco2col0, jaco2col0e, 1e-15));
      t.check(isApproxSame(jaco2col1, jaco2col1e, 1e-15));

      // Check untransformed derivatives
      const auto jaco2un     = localF.evaluateDerivative(gpIndex, wrt(spatialAll));
      const auto jaco2col0un = localF.evaluateDerivative(gpIndex, wrt(spatial(0)));
      const auto jaco2col1un = localF.evaluateDerivative(gpIndex, wrt(spatial(1)));
      t.check(isApproxSame(jaco2un.col(0), jaco2col0un, 1e-15));
      t.check(isApproxSame(jaco2un.col(1), jaco2col1un, 1e-15));

      auto func  = [&](auto& gpOffset_) { return localF.evaluateFunction(toDune(gpOffset_)); };
      auto deriv = [&](auto& gpOffset_) { return localF.evaluateDerivative(toDune(gpOffset_), wrt(spatialAll)); };
      Eigen::Vector<double, 2> gpOffset = toEigen(gp.position());
      auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(func, deriv), parameter(gpOffset));

      t.check((checkJacobian<decltype(nonLinOp), Eigen::Vector<double, 2>>(
          nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = false})));

      auto localFdual_ = [&](auto& x) { return localFdual(x, gpIndex); };
      Eigen::VectorXdual xv(vasMat.cols() * vasMat.rows());
      xv.setZero();
      const Eigen::MatrixXd Jdual = jacobian(localFdual_, autodiff::wrt(xv), at(xv));

      for (size_t i = 0; i < fe.size(); ++i) {
        const auto jacobianWRTCoeffs = localF.evaluateDerivative(gpIndex, wrt(coeff(i)));
        t.check(isApproxSame(jacobianWRTCoeffs, Jdual.block<size, size>(0, i * size), 1e-15));

        const auto Warray  = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatialAll), transformWith(Jinv));
        const auto Warray2 = localF.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i)), transformWith(Jinv));
        for (int j = 0; j < 2; ++j)
          t.check(isApproxSame(Warray[j], Warray2[j], 1e-15));

        const auto W0 = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(0)), transformWith(Jinv));
        const auto W1 = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(1)), transformWith(Jinv));

        t.check(isApproxSame(Warray[0], W0, 1e-15));
        t.check(isApproxSame(Warray[1], W1, 1e-15));

        const auto Warrayun  = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatialAll));
        const auto Warray2un = localF.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i)));
        for (int j = 0; j < 2; ++j)
          t.check(isApproxSame(Warrayun[j], Warray2un[j], 1e-15));

        const auto W0un = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(0)));
        const auto W1un = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(1)));

        t.check(isApproxSame(Warrayun[0], W0un, 1e-15));
        t.check(isApproxSame(Warrayun[1], W1un, 1e-15));
      }

      t.check(isApproxSame(directorCached, directoreval, 1e-15));

      ++gpIndex;
    }
  }
  return t;
}

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t;

  t.subTest(standardLocalFunctionTest<1>());
  t.subTest(standardLocalFunctionTest<2>());
  t.subTest(standardLocalFunctionTest<3>());
  t.subTest(projectionBasedLocalFunctionTest<2>());
  t.subTest(projectionBasedLocalFunctionTest<3>());
  t.subTest(projectionBasedLocalFunctionTest<4>());
  t.subTest(projectionBasedLocalFunctionTest<5>());

  return t.exit();
}