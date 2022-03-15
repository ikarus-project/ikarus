//
// Created by Alex on 21.04.2021.
//
#include "../../config.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"

#include <array>
#include <complex>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <Eigen/Core>

#include "ikarus/Geometries/GeometryWithExternalInput.h"
#include "ikarus/LocalBasis/localBasis.h"
#include "ikarus/utils/utils/algorithms.h"
#include <ikarus/Geometries/ProjectionBasedLocalFunction.h>
#include <ikarus/Geometries/SimpleLocalFunction.h>
#include <ikarus/Variables/VariableDefinitions.h>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

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
    for (size_t i = 0; i < fe.size(); ++i) {
      auto globalIndex = localView.index(localView.tree().child(0).localIndex(i));
      vBlockedLocal[i] = vBlocked[globalIndex[0]];
      vBlockedLocalDual[i].setValue(vBlocked[globalIndex[0]].getValue());
    }
    auto localF = Ikarus::ProjectionBasedLocalFunction(localBasis, vBlockedLocal);
    auto localFdual = [&](auto&x,int gpI){
      Dune::BlockVector<Ikarus::UnitVector<autodiff::dual, 2>>  v=vBlockedLocalDual;

      Ikarus::LinearAlgebra::viewAsFlatEigenVector(v) += x;
      auto localF_ = Ikarus::ProjectionBasedLocalFunction(localBasis, v) ;
      return localF_.evaluateFunction(gpI).getValue();
    };

    const auto vasMat = Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(vBlockedLocal);
    using namespace Ikarus::DerivativeDirections;
    for (int gpIndex = 0; auto& gp : rule) {

      const auto& directorCached = localF.evaluateFunction(gpIndex).getValue();
      const Eigen::Vector2d directorEmbedded = vasMat*localBasis.getFunction(gpIndex);
      const auto& directoreval   = localF.evaluateFunction(gp.position());

      const auto& jaco           = localF.evaluateDerivative(gpIndex, directoreval);
      const auto jaco2           = localF.evaluateDerivative(gpIndex, wrt(spatial));
      auto localFdual_ = [&](auto&x){
        return localFdual(x,gpIndex);
      };
      Eigen::VectorXdual xv(vasMat.cols()*vasMat.rows());
      xv.setZero();
      const Eigen::MatrixXd Jdual = jacobian(localFdual_, autodiff::wrt(xv), at(xv));
      for (size_t i = 0; i < fe.size(); ++i) {
        const auto jacobianWRTCoeffs           = localF.evaluateDerivative(gpIndex, wrt(coeffs), coeffIndices(i));
        EXPECT_THAT(jacobianWRTCoeffs, EigenApproxEqual(Jdual.block<2,2>(0,i*2), 1e-15));

        const auto Warray           = localF.evaluateDerivative(gpIndex, wrt(coeffs,spatial), coeffIndices(i));
        const auto Warray2           = localF.evaluateDerivative(gpIndex, wrt(spatial,coeffs), coeffIndices(i));
        for (int j = 0; j < 2; ++j) {
          EXPECT_THAT(Warray[j], EigenApproxEqual(Warray2[j], 1e-15));
        }

        EXPECT_THAT(jacobianWRTCoeffs, EigenApproxEqual(Ikarus::UnitVector<double, 2>::derivativeOfProjectionWRTposition(directorEmbedded)*localBasis.getFunction(gpIndex)[i], 1e-15));
      }

      EXPECT_DOUBLE_EQ(directorCached.norm(), 1.0);
      EXPECT_DOUBLE_EQ(directoreval.getValue().norm(), 1.0);
      EXPECT_THAT(directorCached, EigenApproxEqual(directoreval.getValue(), 1e-15));
      EXPECT_NEAR((directoreval.getValue().transpose() * jaco).norm(), 0.0, 1e-15);
      EXPECT_NEAR(
          (Ikarus::UnitVector<double, 2>::derivativeOfProjectionWRTposition(directoreval.getValue()) * directoreval.getValue())
              .norm(),
          0.0, 1e-15);

      ++gpIndex;
    }
  }
}
