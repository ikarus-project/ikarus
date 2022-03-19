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
#include "ikarus/LocalBasis/localBasis.h"
#include "ikarus/LocalFunctions/GeometryWithExternalInput.h"
#include "ikarus/utils/utils/algorithms.h"
#include <ikarus/LocalFunctions/ProjectionBasedLocalFunction.h>
#include <ikarus/LocalFunctions/SimpleLocalFunction.h>
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/utils/functionSanityChecks.h>

const double tol = 1e-15;
using namespace Dune::Functions::BasisFactory;
TEST(LocalFunction, ProjectionBasedUnitVector) {
  using namespace Ikarus;
  using Grid        = Dune::YaspGrid<2>;
  const double Lx   = 7;
  const double Ly   = 5;
  const size_t elex = 1;
  const size_t eley = 1;

  Dune::FieldVector<double, 2> bbox = {Lx, Ly};
  std::array<int, 2> eles           = {elex, eley};
  auto grid                         = std::make_shared<Grid>(bbox, eles);

  auto gridView = grid->leafGridView();
  auto basis    = makeBasis(gridView, power<3>(lagrange<1>(), BlockedInterleaved()));
  Dune::BlockVector<Ikarus::UnitVector<double, 2>> vBlocked(basis.size());
  for (auto& vsingle : vBlocked) {
    vsingle.setValue(0.1 * Eigen::Vector<double, 2>::Random() + Eigen::Vector<double, 2>::UnitX());
  }
  auto localView = basis.localView();
  for (auto& ele : elements(gridView)) {
    localView.bind(ele);
    const auto& fe   = localView.tree().child(0).finiteElement();
    auto localBasis  = Ikarus::LocalBasis(localView.tree().child(0).finiteElement().localBasis());
    const auto& rule = Dune::QuadratureRules<double, 2>::rule(localView.element().type(), 3);
    localBasis.bind(rule, 0, 1);
    Dune::BlockVector<Ikarus::UnitVector<double, 2>> vBlockedLocal(fe.size());
    Dune::BlockVector<Ikarus::UnitVector<autodiff::dual, 2>> vBlockedLocalDual(fe.size());
    Dune::BlockVector<Ikarus::UnitVector<autodiff::dual2nd, 2>> vBlockedLocalDual2nd(fe.size());
    for (size_t i = 0; i < fe.size(); ++i) {
      auto globalIndex = localView.index(localView.tree().child(0).localIndex(i));
      vBlockedLocal[i] = vBlocked[globalIndex[0]];
      vBlockedLocalDual[i].setValue(vBlocked[globalIndex[0]].getValue());
      vBlockedLocalDual2nd[i].setValue(vBlocked[globalIndex[0]].getValue());
    }
    auto localF     = Ikarus::ProjectionBasedLocalFunction(localBasis, vBlockedLocal);
    auto localFdual = [&](auto& x, int gpI) {
      Dune::BlockVector<Ikarus::UnitVector<autodiff::dual, 2>> v = vBlockedLocalDual;

      Ikarus::LinearAlgebra::viewAsFlatEigenVector(v) += x;
      auto localF_ = Ikarus::ProjectionBasedLocalFunction(localBasis, v);
      return localF_.evaluateFunction(gpI).getValue();
    };

//    auto localFdual2nd = [&](auto& x, int gpI) {
//      Dune::BlockVector<Ikarus::UnitVector<autodiff::dual2nd, 2>> v = vBlockedLocalDual2nd;
//
//      Ikarus::LinearAlgebra::viewAsFlatEigenVector(v) += x;
//      auto localF_ = Ikarus::ProjectionBasedLocalFunction(localBasis, v);
//      return localF_.evaluateFunction(gpI).getValue();
//    };

    const auto vasMat = Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(vBlockedLocal);
    using namespace Ikarus::DerivativeDirections;
    for (int gpIndex = 0; auto& gp : rule) {
      const auto& directorCached             = localF.evaluateFunction(gpIndex).getValue();
      const Eigen::Vector2d directorEmbedded = vasMat * localBasis.getFunction(gpIndex);
      const auto& directoreval               = localF.evaluateFunction(gp.position());

      const auto J    = toEigenMatrix(ele.geometry().jacobianTransposed(gp.position())).transpose().eval();
      const auto Jinv = J.inverse().eval();
      //      const auto& jaco = localF.evaluateDerivative(gpIndex, directoreval, wrt(spatial<all>),
      //      transformWith(Jinv));
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
      const Eigen::MatrixXd Jdual   = jacobian(localFdual_, autodiff::wrt(xv), at(xv));

//      auto localFdual2nd_ = [&](auto& x) { return localFdual(x, gpIndex)[0]; };
//      Eigen::VectorXdual2nd xv2nd(vasMat.cols() * vasMat.rows());
//      xv.setZero();
//      const Eigen::MatrixXd Jdual2nd   = hessian(localFdual2nd_, autodiff::wrt(xv2nd), at(xv2nd));


      const Eigen::Vector2d testVec = Eigen::Vector2d::UnitX();
      for (size_t i = 0; i < fe.size(); ++i) {
        const auto jacobianWRTCoeffs = localF.evaluateDerivative(gpIndex, wrt(coeffs), coeffIndices(i));
        EXPECT_THAT(jacobianWRTCoeffs, EigenApproxEqual(Jdual.block<2, 2>(0, i * 2), 1e-15));

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
        //        const auto Warray2 = localF.evaluateDerivative(gpIndex, wrt(spatial<all>, coeffs), coeffIndices(i));
        const auto S
            = localF.evaluateDerivative(gpIndex, wrt(coeffs, coeffs), along(testVec), coeffIndices(i, i));

//        EXPECT_THAT(S, EigenApproxEqual(Jdual2nd.block<2, 2>(0, 0), 1e-15));
        ////        const auto Warray4
        ////            = localF.evaluateDerivative(gpIndex, wrt(spatial<all>, coeffs, coeffs), along(testVec),
        /// coeffIndices(i, i));
        //        std::cout << Warray4[0] << std::endl;
        //        std::cout << Warray4[1] << std::endl;
        //        for (int j = 0; j < 2; ++j) {
        //          EXPECT_THAT(Warray[j], EigenApproxEqual(Warray2[j], 1e-15));
        //        }
        //
        //        const auto& dN = localBasis.getJacobian(gpIndex);
        //        auto JacoEmbedded     = (Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(vBlockedLocal) * dN).eval();
        //        for (int dir = 0; dir < 2; ++dir) {
        //          EXPECT_THAT(Warray2[dir],
        //                      EigenApproxEqual(Ikarus::UnitVector<double,
        //                      2>::secondDerivativeOfProjectionWRTposition(directorEmbedded, JacoEmbedded.col(dir)) *
        //                      localBasis.getFunction(gpIndex)[i]
        //                                           + Ikarus::UnitVector<double,
        //                                           2>::derivativeOfProjectionWRTposition(directorEmbedded)
        //                                                 * localBasis.getJacobian(gpIndex)(i, dir)
        //                                           , 1e-15));
        //        }

        //        EXPECT_THAT(jacobianWRTCoeffs,
        //                    EigenApproxEqual(Ikarus::UnitVector<double,
        //                    2>::derivativeOfProjectionWRTposition(directorEmbedded)
        //                                         * localBasis.getFunction(gpIndex)[i],
        //                                     1e-15));
      }

      EXPECT_DOUBLE_EQ(directorCached.norm(), 1.0);
      EXPECT_DOUBLE_EQ(directoreval.getValue().norm(), 1.0);
      EXPECT_THAT(directorCached, EigenApproxEqual(directoreval.getValue(), 1e-15));
      //      EXPECT_NEAR((directoreval.getValue().transpose() * jaco).norm(), 0.0, 1e-15);
      EXPECT_NEAR((Ikarus::UnitVector<double, 2>::derivativeOfProjectionWRTposition(directoreval.getValue())
                   * directoreval.getValue())
                      .norm(),
                  0.0, 1e-15);

      ++gpIndex;
    }
  }
}
