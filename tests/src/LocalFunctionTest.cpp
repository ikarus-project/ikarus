//
// Created by Alex on 21.04.2021.
//

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
#include <ikarus/Geometries/ProjectionBasedLocalFunction.h>
#include <ikarus/Geometries/SimpleLocalFunction.h>
#include <ikarus/Variables/VariableDefinitions.h>

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
    for (size_t i = 0; i < fe.size(); ++i) {
      auto globalIndex = localView.index(localView.tree().child(0).localIndex(i));
      vBlockedLocal[i] = vBlocked[globalIndex[0]];
    }
    auto localF = Ikarus::ProjectionBasedLocalFunction(localBasis, vBlockedLocal);
    const auto vasMat = Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(vBlockedLocal);
    for (int gpIndex = 0; auto& gp : rule) {

      const auto& directorCached = localF.evaluateFunction(gpIndex).getValue();
      const Eigen::Vector2d directorEmbedded = vasMat*localBasis.getFunction(gpIndex);
      const auto& directoreval   = localF.evaluateFunction(gp.position());
      const auto& jaco           = localF.evaluateDerivative(gpIndex, directoreval);
      const auto jaco2           = localF.evaluateDerivative(gpIndex, wrt(spatial{}));
      for (size_t i = 0; i < fe.size(); ++i) {
        const auto jaco3           = localF.evaluateDerivative(gpIndex, wrt(coeffs{}), {i});
        EXPECT_THAT(jaco3, EigenApproxEqual(Ikarus::UnitVector<double, 2>::derivativeOfProjection(directorEmbedded)*localBasis.getFunction(gpIndex)[i], 1e-15));
      }

      EXPECT_DOUBLE_EQ(directorCached.norm(), 1.0);
      EXPECT_DOUBLE_EQ(directoreval.getValue().norm(), 1.0);
      EXPECT_THAT(directorCached, EigenApproxEqual(directoreval.getValue(), 1e-15));
      //      std::cout<<jaco*directoreval.getValue()<<std::endl;
      EXPECT_NEAR((directoreval.getValue().transpose() * jaco).norm(), 0.0, 1e-15);
      EXPECT_NEAR(
          (Ikarus::UnitVector<double, 2>::derivativeOfProjection(directoreval.getValue()) * directoreval.getValue())
              .norm(),
          0.0, 1e-15);
      //      std::cout<<vBlocked<<std::endl;
      //      std::cout<<Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(vBlocked)<<std::endl;
      ++gpIndex;
    }
  }
}
