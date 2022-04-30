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
#include <ikarus/manifolds/realTuple.hh>
#include <ikarus/manifolds/unitVector.hh>
#include <ikarus/utils/functionSanityChecks.hh>

using namespace Dune::Functions::BasisFactory;

template <typename LF,std::size_t... I>
auto getCoeffRefHelper(const LF& lf)
{
  if constexpr(LF::isLeaf)
    return lf.coefficientsRef()
  else
    return collectLeafNodeLocalFunctions(lf).coefficientsRef();
}

template <typename LF>
void testLocalFunction(const LF& lf, int ipIndex) {
  using namespace Ikarus::DerivativeDirections;
  using namespace autodiff;
  using namespace Ikarus;

  const auto& coeffs  = getCoeffRefHelper(lf);
  const auto& basis  = getBasisHelper(lf);
  const size_t coeffSize = coeffs.size();
  const int gridDim = LF::gridDim;
  using Manifold = typename std::remove_cvref_t<decltype(coeffs)>::value_type;
  constexpr int coeffValueSize = Manifold::valueSize;
  constexpr int coeffCorrectionSize = Manifold::correctionSize;

  typename LF::template Rebind<dual>::other lfDual(basis(), convertUnderlying<dual>(coeffs));

  auto localFdual = [&](auto& x) {
    lfDual.coefficientsRef() += x;
    auto value =  lfDual.evaluateFunction(ipIndex).getValue();
    lfDual.coefficientsRef() -= x;
    return value;
  };

  auto localFdualSpatialSingle = [&](auto& x,int i) {
    lfDual.coefficientsRef() += x;
    auto value =  lfDual.evaluateDerivative(ipIndex,Ikarus::wrt(spatial(i)));
    lfDual.coefficientsRef() -= x;
    return value;
  };

  Eigen::VectorXdual xv(correctionSize(lfDual.coefficientsRef()));
  xv.setZero();
  const Eigen::MatrixXd Jcoeff = jacobian(localFdual, autodiff::wrt(xv), at(xv));
  const Eigen::MatrixXd JspatialSingleCoeff = jacobian(localFdual, autodiff::wrt(xv), at(xv));

  std::array<Eigen::MatrixXd,gridDim> JdualSpatialWRTCoeffs;
  for (int d = 0; d < gridDim; ++d)
   JdualSpatialWRTCoeffs[d] = jacobian(localFdualSpatialSingle, autodiff::wrt(xv), at(xv,d));

  //  const Eigen::Vector<double, size> testVec = Eigen::Vector<double, size>::UnitX();
  for (size_t i = 0; i < coeffSize; ++i) {
    const auto jacobianWRTCoeffs = lf.evaluateDerivative(ipIndex, Ikarus::wrt(coeff(i)));
    EXPECT_THAT(
        jacobianWRTCoeffs,
        EigenApproxEqual(Jcoeff.block<coeffValueSize, coeffCorrectionSize>(0, i * coeffValueSize), 1e-15));

    for (int d = 0; d < gridDim; ++d) {
      const auto jacoWrtCoeffAndSpatial = lf.evaluateDerivative(ipIndex, Ikarus::wrt(coeff(i),spatial(d)));
      const auto jacoWrtSpatialAndCoeff = lf.evaluateDerivative(ipIndex, Ikarus::wrt(spatial(d),coeff(i)));
      EXPECT_THAT(  jacoWrtCoeffAndSpatial,  EigenApproxEqual(jacoWrtSpatialAndCoeff, 1e-15));

      EXPECT_THAT(jacoWrtCoeffAndSpatial,
          EigenApproxEqual(JdualSpatialWRTCoeffs[d].template block<coeffValueSize, coeffCorrectionSize>(0, i * coeffValueSize), 1e-15));
    }
    }
}

TEST(LocalFunctionTests, TestExpressions) {
  constexpr int sizeD = 3;
  using Manifold      = Ikarus::RealTuple<double, sizeD>;
  using Manifold2     = Ikarus::UnitVector<double, sizeD>;
  using VectorType    = Eigen::Vector<double, sizeD>;
  using MatrixType    = Eigen::Matrix<double, sizeD, sizeD>;
  constexpr int size  = Manifold::valueSize;
  using namespace Ikarus;
  using namespace Dune::Indices;
  auto grid = createGrid<Grids::Yasp>();

  auto gridView = grid->leafGridView();
  auto basis    = makeBasis(gridView, lagrange<2>());
  Dune::BlockVector<Manifold> vBlocked(basis.size());
  Dune::BlockVector<Manifold> vBlocked2(basis.size());
  Dune::BlockVector<Manifold2> vBlocked3(basis.size());
  for (auto& vsingle : vBlocked)
    vsingle.setValue(0.1 * Eigen::Vector<double, size>::Random() + Eigen::Vector<double, size>::UnitX());

  for (auto& vsingle : vBlocked2)
    vsingle.setValue(0.5 * Eigen::Vector<double, size>::Random() + 7 * Eigen::Vector<double, size>::UnitX());

  for (auto& vsingle : vBlocked3)
    vsingle.setValue(0.25 * Eigen::Vector<double, size>::Random() + 37 * Eigen::Vector<double, size>::UnitX());

  auto localView = basis.localView();
  for (auto& ele : elements(gridView)) {
    localView.bind(ele);
    const auto& fe   = localView.tree().finiteElement();
    auto localBasis  = Ikarus::LocalBasis(localView.tree().finiteElement().localBasis());
    const auto& rule = Dune::QuadratureRules<double, 2>::rule(localView.element().type(), 3);
    localBasis.bind(rule, bindDerivatives(0, 1));
    Dune::BlockVector<Manifold> vBlockedLocal(fe.size());
    Dune::BlockVector<Manifold> vBlockedLocal2(fe.size());
    Dune::BlockVector<Manifold2> vBlockedLocal3(fe.size());

    for (size_t i = 0; i < fe.size(); ++i) {
      auto globalIndex  = localView.index(localView.tree().localIndex(i));
      vBlockedLocal[i]  = vBlocked[globalIndex[0]];
      vBlockedLocal2[i] = vBlocked2[globalIndex[0]];
      vBlockedLocal3[i] = vBlocked3[globalIndex[0]];
    }
    auto f = Ikarus::StandardLocalFunction(localBasis, vBlockedLocal);
    auto g = Ikarus::StandardLocalFunction(localBasis, vBlockedLocal);


    static_assert(countUniqueNonArithmeticLeafNodes(f) == 1);
    static_assert(countUniqueNonArithmeticLeafNodes(g) == 1);
    using namespace Ikarus::DerivativeDirections;
    auto h = f + g;

    auto a     = collectNonArithmeticLeafNodes(h);
    auto hLeaf = collectLeafNodeLocalFunctions(h);
    static_assert(std::tuple_size_v<decltype(a)> == 2);

    static_assert(countUniqueNonArithmeticLeafNodes(h) == 1);
    static_assert(std::is_same_v<decltype(h)::Ids, std::tuple<Dune::index_constant<0>, Dune::index_constant<0>>>);

    for (size_t i = 0; i < fe.size(); ++i) {
      EXPECT_TRUE(collectLeafNodeLocalFunctions(h).getCoeffs(_0)[i] == vBlockedLocal[i]);
      EXPECT_TRUE(collectLeafNodeLocalFunctions(h).getCoeffs(_1)[i] == vBlockedLocal[i]);
    }

    for (int gpIndex = 0; auto& gp : rule) {
      testLocalFunction(f,gpIndex);
      testLocalFunction(h,gpIndex);

      EXPECT_THAT(f.evaluateFunction(gpIndex).getValue() + g.evaluateFunction(gpIndex).getValue(),
                  EigenApproxEqual(h.evaluateFunction(gpIndex).getValue(), 1e-15));
      EXPECT_THAT(f.evaluateDerivative(gpIndex, wrt(spatial(0))) + g.evaluateDerivative(gpIndex, wrt(spatial(0))),
                  EigenApproxEqual(h.evaluateDerivative(gpIndex, wrt(spatial(0))), 1e-15));
      EXPECT_THAT(f.evaluateDerivative(gpIndex, wrt(spatialAll)) + g.evaluateDerivative(gpIndex, wrt(spatialAll)),
                  EigenApproxEqual(h.evaluateDerivative(gpIndex, wrt(spatialAll)), 1e-15));
      for (size_t i = 0; i < fe.size(); ++i) {
        const MatrixType dfgdi
            = f.evaluateDerivative(gpIndex, wrt(coeff(i))) + g.evaluateDerivative(gpIndex, wrt(coeff(i)));
        const MatrixType dhdi    = h.evaluateDerivative(gpIndex, wrt(coeff(i)));
        const MatrixType dhdSdi  = h.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(i)));
        const MatrixType dfgdSdi = f.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(i)))
                                   + g.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(i)));
//        EXPECT_THAT(dfgdi, EigenApproxEqual(dhdi, 1e-14));
        EXPECT_THAT(dhdSdi, EigenApproxEqual(dfgdSdi, 1e-14));
        for (size_t j = 0; j < fe.size(); ++j) {
          const MatrixType dfgdj
              = f.evaluateDerivative(gpIndex, wrt(coeff(j))) + g.evaluateDerivative(gpIndex, wrt(coeff(j)));
          const MatrixType dhdj    = h.evaluateDerivative(gpIndex, wrt(coeff(j)));
          const MatrixType dhdSdj  = h.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(j)));
          const MatrixType dfgdSdj = f.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(j)))
                                     + g.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(j)));
          if (i == j) {
            EXPECT_THAT(dfgdj, EigenApproxEqual(dhdj, 1e-14));
            EXPECT_THAT(dfgdSdj, EigenApproxEqual(dhdSdj, 1e-14));
          }
        }
      }
      ++gpIndex;
    }

    auto k = -dot(f + f, 3.0 * (g / 5.0) * 5.0);
    auto b = collectNonArithmeticLeafNodes(k);
    static_assert(std::tuple_size_v<decltype(b)> == 3);

    //    std::cout<<Dune::className(a)<<std::endl;
//    std::cout << Dune::className(b) << std::endl;
    static_assert(countUniqueNonArithmeticLeafNodes(k) == 1);
    static_assert(std::is_same_v<decltype(k)::Ids, std::tuple<Dune::index_constant<0>, Dune::index_constant<0>,
                                                              Ikarus::Arithmetic, Dune::index_constant<0>>>);

    const double tol = 1e-13;
    for (int gpIndex = 0; auto& gp : rule) {
      const auto& N  = localBasis.evaluateFunction(gpIndex);
      const auto& dN = localBasis.evaluateJacobian(gpIndex);
      EXPECT_DOUBLE_EQ((-2 * 3) * f.evaluateFunction(gpIndex).getValue().dot(g.evaluateFunction(gpIndex).getValue()),
                       k.evaluateFunction(gpIndex));
      auto resSingleSpatial = ((-2 * 3) * f.evaluateDerivative(gpIndex, wrt(spatial(0))).transpose()
                                   * g.evaluateFunction(gpIndex).getValue()
                               + (-2 * 3) * f.evaluateFunction(gpIndex).getValue().transpose()
                                     * g.evaluateDerivative(gpIndex, wrt(spatial(0))))
                                  .eval();
      EXPECT_THAT(resSingleSpatial, EigenApproxEqual(k.evaluateDerivative(gpIndex, wrt(spatial(0))), tol));
      auto resSpatialAll = (((-2 * 3) * f.evaluateDerivative(gpIndex, wrt(spatialAll)).transpose()
                             * g.evaluateFunction(gpIndex).getValue())
                                .transpose()
                            + (-2 * 3) * f.evaluateFunction(gpIndex).getValue().transpose()
                                  * g.evaluateDerivative(gpIndex, wrt(spatialAll)))
                               .eval();

      static_assert(resSpatialAll.cols() == 2);
      static_assert(resSpatialAll.rows() == 1);

      EXPECT_THAT(resSpatialAll, EigenApproxEqual(k.evaluateDerivative(gpIndex, wrt(spatialAll)), tol));
      for (size_t i = 0; i < fe.size(); ++i) {
        const VectorType dfgdi = ((transpose((-2 * 3) * f.evaluateDerivative(gpIndex, wrt(coeff(i))))
                                   * g.evaluateFunction(gpIndex).getValue())
                                      .transpose()
                                  + (-2 * 3) * f.evaluateFunction(gpIndex).getValue().transpose()
                                        * g.evaluateDerivative(gpIndex, wrt(coeff(i))))
                                     .eval();

        const VectorType dkdi = k.evaluateDerivative(gpIndex, wrt(coeff(i)));

        EXPECT_THAT(dfgdi, EigenApproxEqual(dkdi, tol));

        const VectorType dkdSdi = k.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(i)));

        const VectorType dkdSdi2 = k.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(1)));

        EXPECT_THAT(dkdSdi, EigenApproxEqual(dkdSdi2, tol));

        for (size_t j = 0; j < fe.size(); ++j) {
          const MatrixType dkdSdi3         = k.evaluateDerivative(gpIndex, wrt(coeff(i, j)));
          const MatrixType dkdSdi3Expected = -12 * N[i] * N[j] * MatrixType::Identity();
          EXPECT_THAT(dkdSdi3, EigenApproxEqual(dkdSdi3Expected, tol));

          const MatrixType dkdSdi4 = k.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(i, j)));
          const MatrixType dkdSdi5 = k.evaluateDerivative(gpIndex, wrt(coeff(i, j), spatial(1)));
          const MatrixType dkdSdi6 = k.evaluateDerivative(gpIndex, wrt(coeff(j, i), spatial(1)));
          EXPECT_THAT(dkdSdi4, EigenApproxEqual(dkdSdi5, tol));
          EXPECT_THAT(dkdSdi5, EigenApproxEqual(dkdSdi6, tol));
          const MatrixType dkdSdi4Expected = -12 * (dN(i, 1) * N[j] + dN(j, 1) * N[i]) * MatrixType::Identity();
          EXPECT_THAT(dkdSdi4, EigenApproxEqual(dkdSdi4Expected, tol));
        }
      }
      ++gpIndex;
    }

    auto f2 = Ikarus::StandardLocalFunction(localBasis, vBlockedLocal, _0);
    auto g2 = Ikarus::StandardLocalFunction(localBasis, vBlockedLocal2, _1);
    static_assert(countUniqueNonArithmeticLeafNodes(f2) == 1);
    static_assert(countUniqueNonArithmeticLeafNodes(g2) == 1);

    auto k2 = dot(f2 + g2, g2);
    static_assert(countUniqueNonArithmeticLeafNodes(k2) == 2);
    static_assert(
        std::is_same_v<decltype(k2)::Ids,
                       std::tuple<Dune::index_constant<0>, Dune::index_constant<1>, Dune::index_constant<1>>>);

    auto b2 = collectNonArithmeticLeafNodes(k2);
    static_assert(std::tuple_size_v<decltype(b2)> == 3);

    for (int gpIndex = 0; auto& gp : rule) {
      const auto& N  = localBasis.evaluateFunction(gpIndex);
      const auto& dN = localBasis.evaluateJacobian(gpIndex);
      EXPECT_DOUBLE_EQ((f2.evaluateFunction(gpIndex).getValue() + g2.evaluateFunction(gpIndex).getValue())
                           .dot(g2.evaluateFunction(gpIndex).getValue()),
                       k2.evaluateFunction(gpIndex));
      auto resSingleSpatial
          = ((f2.evaluateDerivative(gpIndex, wrt(spatial(0))) + g2.evaluateDerivative(gpIndex, wrt(spatial(0))))
                     .transpose()
                 * g2.evaluateFunction(gpIndex).getValue()
             + (f2.evaluateFunction(gpIndex).getValue() + g2.evaluateFunction(gpIndex).getValue()).transpose()
                   * g2.evaluateDerivative(gpIndex, wrt(spatial(0))))
                .eval();
      EXPECT_THAT(resSingleSpatial, EigenApproxEqual(k2.evaluateDerivative(gpIndex, wrt(spatial(0))), tol));
      auto resSpatialAll
          = (((f2.evaluateDerivative(gpIndex, wrt(spatialAll)) + g2.evaluateDerivative(gpIndex, wrt(spatialAll)))
                  .transpose()
              * g2.evaluateFunction(gpIndex).getValue())
                 .transpose()
             + (f2.evaluateFunction(gpIndex).getValue() + g2.evaluateFunction(gpIndex).getValue()).transpose()
                   * g2.evaluateDerivative(gpIndex, wrt(spatialAll)))
                .eval();
      static_assert(resSpatialAll.cols() == 2);
      static_assert(resSpatialAll.rows() == 1);

      EXPECT_THAT(resSpatialAll, EigenApproxEqual(k2.evaluateDerivative(gpIndex, wrt(spatialAll)), tol));

      for (size_t i = 0; i < fe.size(); ++i) {
        const VectorType dfdi = g2.evaluateFunction(gpIndex).getValue() * N[i];

        const VectorType dkdi = k2.evaluateDerivative(gpIndex, wrt(coeff(_0, i)));

        EXPECT_THAT(dfdi, EigenApproxEqual(dkdi, tol));

        //        const VectorType dkdSdi  = k.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(i)));
        //        const VectorType dkdSdi2 = k.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(1)));
        //
        //        EXPECT_THAT(dkdSdi, EigenApproxEqual(dkdSdi2, tol));
        //
        for (size_t j = 0; j < fe.size(); ++j) {
          const MatrixType dkdij         = k2.evaluateDerivative(gpIndex, wrt(coeff(_0, i, _1, j)));
          const MatrixType dkdijExpected = N[j] * N[i] * MatrixType::Identity();
          EXPECT_THAT(dkdijExpected, EigenApproxEqual(dkdij, tol));

          const MatrixType dkdij2         = k2.evaluateDerivative(gpIndex, wrt(coeff(_0, i, _0, j)));
          const MatrixType dkdijExpected2 = MatrixType::Zero();
          EXPECT_THAT(dkdijExpected2, EigenApproxEqual(dkdij2, tol));
          const MatrixType dkdij3         = k2.evaluateDerivative(gpIndex, wrt(coeff(_1, i, _1, j)));
          const MatrixType dkdijExpected3 = 2 * N[i] * N[j] * MatrixType::Identity();
          EXPECT_THAT(dkdijExpected3, EigenApproxEqual(dkdij3, tol));
          //          const MatrixType dkdSdi3 = k.evaluateDerivative(gpIndex, wrt(coeff(i, j)));
          //
          //          EXPECT_DOUBLE_EQ(dkdSdi3(0, 0), -12 * N[i] * N[j]);
          //          EXPECT_DOUBLE_EQ(dkdSdi3(1, 1), -12 * N[i] * N[j]);
          //          EXPECT_DOUBLE_EQ(dkdSdi3(0, 1), 0);
          //          EXPECT_DOUBLE_EQ(dkdSdi3(1, 0), 0);

          const MatrixType dkdSij         = k2.evaluateDerivative(gpIndex, wrt(spatial(0), coeff(_0, i, _1, j)));
          const MatrixType dkdSijR        = k2.evaluateDerivative(gpIndex, wrt(coeff(_0, i, _1, j), spatial(0)));
          const MatrixType dkdSijExpected = (dN(j, 0) * N[i] + N[j] * dN(i, 0)) * MatrixType::Identity();
          EXPECT_THAT(dkdSijR, EigenApproxEqual(dkdSij, tol));
          EXPECT_THAT(dkdSijExpected, EigenApproxEqual(dkdSij, tol));
          //          const MatrixType dkdSdi4 = k.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(i, j)));
          //          const MatrixType dkdSdi5 = k.evaluateDerivative(gpIndex, wrt(coeff(i, j), spatial(1)));
          //          const MatrixType dkdSdi6 = k.evaluateDerivative(gpIndex, wrt(coeff(j, i), spatial(1)));
          //          EXPECT_THAT(dkdSdi4, EigenApproxEqual(dkdSdi5, tol));
          //          EXPECT_THAT(dkdSdi5, EigenApproxEqual(dkdSdi6, tol));
          //          EXPECT_NEAR(dkdSdi4(0, 0), -12 * (dN(i, 1) * N[j]+dN(j, 1) * N[i]),tol);
          //          EXPECT_NEAR(dkdSdi4(1, 1), -12 * (dN(i, 1) * N[j]+dN(j, 1) * N[i]),tol);
          //          EXPECT_DOUBLE_EQ(dkdSdi4(0, 1), 0);
          //          EXPECT_DOUBLE_EQ(dkdSdi4(1, 0), 0);
        }
      }
      ++gpIndex;
    }

    auto gP = Ikarus::ProjectionBasedLocalFunction(localBasis, vBlockedLocal3);
  }
}

template <typename T>
class LocalFunctionProjectionBasedUnitVector : public testing::Test {
public:
  static constexpr int size = T::value;
};
using test_types = ::testing::Types<std::integral_constant<std::size_t, 2>, std::integral_constant<std::size_t, 3>,
                                    std::integral_constant<std::size_t, 4>, std::integral_constant<std::size_t, 5>>;

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
    auto localFdual = [&](auto& x, int gpI) {
      Dune::BlockVector<Ikarus::UnitVector<autodiff::dual, size>> v = vBlockedLocalDual;

      Ikarus::viewAsFlatEigenVector(v) += x;
      auto localF_ = Ikarus::ProjectionBasedLocalFunction(localBasis, v);
      return localF_.evaluateFunction(gpI).getValue();
    };

    const auto vasMat = Ikarus::viewAsEigenMatrixFixedDyn(vBlockedLocal);
    using namespace Ikarus::DerivativeDirections;
    for (int gpIndex = 0; auto& gp : rule) {
      const auto directorCached                          = localF.evaluateFunction(gpIndex).getValue();
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

      auto func = [&](auto& gpOffset_) { return localF.evaluateFunction(toFieldVector(gpOffset_)).getValue(); };
      auto deriv
          = [&](auto& gpOffset_) { return localF.evaluateDerivative(toFieldVector(gpOffset_), wrt(spatialAll)); };
      Eigen::Vector<double, 2> gpOffset = toEigenVector(gp.position());
      auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(func, deriv), parameter(gpOffset));

      EXPECT_TRUE((checkJacobian<decltype(nonLinOp), Eigen::Vector<double, 2>>(nonLinOp, false)));

      auto localFdual_ = [&](auto& x) { return localFdual(x, gpIndex); };
      Eigen::VectorXdual xv(vasMat.cols() * vasMat.rows());
      xv.setZero();
      const Eigen::MatrixXd Jdual = jacobian(localFdual_, autodiff::wrt(xv), at(xv));

      const Eigen::Vector<double, size> testVec = Eigen::Vector<double, size>::UnitX();
      for (size_t i = 0; i < fe.size(); ++i) {
        const auto jacobianWRTCoeffs = localF.evaluateDerivative(gpIndex, wrt(coeff(i)));
        EXPECT_THAT(
            jacobianWRTCoeffs,
            EigenApproxEqual(Jdual.block<size, size>(0, i * size) * vBlockedLocal[i].orthonormalFrame(), 1e-15));

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
        const auto Sun = localF.evaluateDerivative(gpIndex, wrt(coeff(i, i)), along(testVec));

        const auto S = localF.evaluateDerivative(gpIndex, wrt(coeff(i, i)), along(testVec), transformWith(Jinv));

        const auto chi
            = localF.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i, i)), along(testVec), transformWith(Jinv));

        const auto chi0
            = localF.evaluateDerivative(gpIndex, wrt(spatial(0), coeff(i, i)), along(testVec), transformWith(Jinv));

        const auto chi1
            = localF.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(i, i)), along(testVec), transformWith(Jinv));

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
    auto localFdual = [&](auto& x, int gpI) {
      Dune::BlockVector<typename Manifold::template Rebind<autodiff::dual>::other> v = vBlockedLocalDual;

      Ikarus::viewAsFlatEigenVector(v) += x;
      auto localF_ = Ikarus::StandardLocalFunction(localBasis, v);
      return localF_.evaluateFunction(gpI).getValue();
    };

    const auto vasMat = Ikarus::viewAsEigenMatrixFixedDyn(vBlockedLocal);
    using namespace Ikarus::DerivativeDirections;
    for (int gpIndex = 0; auto& gp : rule) {
      const auto& directorCached                         = localF.evaluateFunction(gpIndex).getValue();
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

      auto func = [&](auto& gpOffset_) { return localF.evaluateFunction(toFieldVector(gpOffset_)).getValue(); };
      auto deriv
          = [&](auto& gpOffset_) { return localF.evaluateDerivative(toFieldVector(gpOffset_), wrt(spatialAll)); };
      Eigen::Vector<double, 2> gpOffset = toEigenVector(gp.position());
      auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(func, deriv), parameter(gpOffset));

      EXPECT_TRUE((checkJacobian<decltype(nonLinOp), Eigen::Vector<double, 2>>(nonLinOp, false)));

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

      EXPECT_THAT(directorCached, EigenApproxEqual(directoreval.getValue(), 1e-15));

      ++gpIndex;
    }
  }
}
