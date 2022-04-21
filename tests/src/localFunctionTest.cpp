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

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <Eigen/Core>

#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/localFunctions/projectionBasedLocalFunction.hh>
#include <ikarus/localFunctions/standardLocalFunction.hh>
#include <ikarus/utils/functionSanityChecks.hh>
#include <ikarus/manifolds/unitVector.hh>
#include <ikarus/manifolds/realTuple.hh>
#include <dune/common/classname.hh>

using namespace Dune::Functions::BasisFactory;

TEST(LocalFunctionTests,TestInterface)
{
  using namespace Ikarus::DerivativeDirections;
  auto wrt1 = Ikarus::wrt(spatial(0),coeff(7)) ;

  auto counter = countDerivativesType<decltype(wrt1)>();

  EXPECT_EQ(counter.singleCoeffDerivs,1);
  EXPECT_EQ(counter.twoCoeffDerivs,0);
  EXPECT_EQ(counter.spatialDerivs,1);
  EXPECT_EQ(counter.spatialAll,0);

  auto wrt2 = Ikarus::wrt(spatialall,coeff(7,1)) ;

  counter = countDerivativesType<decltype(wrt2)>();

  EXPECT_EQ(counter.singleCoeffDerivs,0);
  EXPECT_EQ(counter.twoCoeffDerivs,1);
  EXPECT_EQ(counter.spatialDerivs,0);
  EXPECT_EQ(counter.spatialAll,1);
}

//
//
//
//template <typename T>
//class LocalFunctionProjectionBasedUnitVector : public testing::Test {
//public:
//  static constexpr int size = T::value;
//};
//using test_types = ::testing::Types<std::integral_constant<std::size_t, 2>, std::integral_constant<std::size_t, 3>,
//                                    std::integral_constant<std::size_t, 4>, std::integral_constant<std::size_t, 5>>;
//
//TYPED_TEST_SUITE(LocalFunctionProjectionBasedUnitVector, test_types);
//
//TYPED_TEST(LocalFunctionProjectionBasedUnitVector, ProjectionBasedUnitVector) {
//  constexpr int size = TypeParam::value;
//  using namespace Ikarus;
//  auto grid = createGrid<Grids::Alu>();
//
//  auto gridView = grid->leafGridView();
//  auto basis    = makeBasis(gridView, power<3>(lagrange<1>(), BlockedInterleaved()));
//  Dune::BlockVector<Ikarus::UnitVector<double, size>> vBlocked(basis.size());
//  for (auto& vsingle : vBlocked) {
//    vsingle.setValue(0.1 * Eigen::Vector<double, size>::Random() + Eigen::Vector<double, size>::UnitX());
//  }
//  auto localView = basis.localView();
//  for (auto& ele : elements(gridView)) {
//    localView.bind(ele);
//    const auto& fe   = localView.tree().child(0).finiteElement();
//    auto localBasis  = Ikarus::LocalBasis(localView.tree().child(0).finiteElement().localBasis());
//    const auto& rule = Dune::QuadratureRules<double, 2>::rule(localView.element().type(), 3);
//    localBasis.bind(rule, bindDerivatives(0, 1));
//    Dune::BlockVector<Ikarus::UnitVector<double, size>> vBlockedLocal(fe.size());
//    Dune::BlockVector<Ikarus::UnitVector<autodiff::dual, size>> vBlockedLocalDual(fe.size());
//    Dune::BlockVector<Ikarus::UnitVector<autodiff::dual2nd, size>> vBlockedLocalDual2nd(fe.size());
//    for (size_t i = 0; i < fe.size(); ++i) {
//      auto globalIndex = localView.index(localView.tree().child(0).localIndex(i));
//      vBlockedLocal[i] = vBlocked[globalIndex[0]];
//      vBlockedLocalDual[i].setValue(vBlocked[globalIndex[0]].getValue());
//      vBlockedLocalDual2nd[i].setValue(vBlocked[globalIndex[0]].getValue());
//    }
//    auto localF     = Ikarus::ProjectionBasedLocalFunction(localBasis, vBlockedLocal);
//    auto localFdual = [&](auto& x, int gpI) {
//      Dune::BlockVector<Ikarus::UnitVector<autodiff::dual, size>> v = vBlockedLocalDual;
//
//      Ikarus::LinearAlgebra::viewAsFlatEigenVector(v) += x;
//      auto localF_ = Ikarus::ProjectionBasedLocalFunction(localBasis, v);
//      return localF_.evaluateFunction(gpI).getValue();
//    };
//
//    const auto vasMat = Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(vBlockedLocal);
//    using namespace Ikarus::DerivativeDirections;
//    for (int gpIndex = 0; auto& gp : rule) {
//      const auto& directorCached                         = localF.evaluateFunction(gpIndex).getValue();
//      const Eigen::Vector<double, size> directorEmbedded = vasMat * localBasis.evaluateFunction(gpIndex);
//      const auto& directoreval                           = localF.evaluateFunction(gp.position());
//
//      const auto J     = toEigenMatrix(ele.geometry().jacobianTransposed(gp.position())).transpose().eval();
//      const auto Jinv  = J.inverse().eval();
//      const auto jaco2 = localF.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv));
//
//      const auto jaco2e     = localF.evaluateDerivative(gp.position(), wrt(spatialall), transformWith(Jinv));
//      const auto jaco2col0  = localF.evaluateDerivative(gpIndex, wrt(spatial(0)), transformWith(Jinv));
//      const auto jaco2col0e = localF.evaluateDerivative(gp.position(), wrt(spatial(0)), transformWith(Jinv));
//      const auto jaco2col1  = localF.evaluateDerivative(gpIndex, wrt(spatial(1)), transformWith(Jinv));
//      const auto jaco2col1e = localF.evaluateDerivative(gp.position(), wrt(spatial(1)), transformWith(Jinv));
//      EXPECT_THAT(jaco2.col(0), EigenApproxEqual(jaco2col0, 1e-15));
//      EXPECT_THAT(jaco2.col(1), EigenApproxEqual(jaco2col1, 1e-15));
//      EXPECT_THAT(jaco2, EigenApproxEqual(jaco2e, 1e-15));
//      EXPECT_THAT(jaco2col0, EigenApproxEqual(jaco2col0e, 1e-15));
//      EXPECT_THAT(jaco2col1, EigenApproxEqual(jaco2col1e, 1e-15));
//
//      // Check untransformed derivatives
//      const auto jaco2un     = localF.evaluateDerivative(gpIndex, wrt(spatialall));
//      const auto jaco2col0un = localF.evaluateDerivative(gpIndex, wrt(spatial(0)));
//      const auto jaco2col1un = localF.evaluateDerivative(gpIndex, wrt(spatial(1)));
//      EXPECT_THAT(jaco2un.col(0), EigenApproxEqual(jaco2col0un, 1e-15));
//      EXPECT_THAT(jaco2un.col(1), EigenApproxEqual(jaco2col1un, 1e-15));
//
//      auto func = [&](auto& gpOffset_) { return localF.evaluateFunction(toFieldVector(gpOffset_)).getValue(); };
//      auto deriv
//          = [&](auto& gpOffset_) { return localF.evaluateDerivative(toFieldVector(gpOffset_), wrt(spatialall)); };
//      Eigen::Vector<double, 2> gpOffset = toEigenVector(gp.position());
//      auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(func, deriv), parameter(gpOffset));
//
//      EXPECT_TRUE((checkJacobian<decltype(nonLinOp), Eigen::Vector<double, 2>>(nonLinOp, false)));
//
//      auto localFdual_ = [&](auto& x) { return localFdual(x, gpIndex); };
//      Eigen::VectorXdual xv(vasMat.cols() * vasMat.rows());
//      xv.setZero();
//      const Eigen::MatrixXd Jdual = jacobian(localFdual_, autodiff::wrt(xv), at(xv));
//
//      const Eigen::Vector<double, size> testVec = Eigen::Vector<double, size>::UnitX();
//      for (size_t i = 0; i < fe.size(); ++i) {
//        const auto jacobianWRTCoeffs = localF.evaluateDerivative(gpIndex, wrt(coeff(i)));
//        EXPECT_THAT(jacobianWRTCoeffs, EigenApproxEqual((vBlockedLocal[i].orthonormalFrame()).transpose()
//                                                            * Jdual.block<size, size>(0, i * size),
//                                                        1e-15));
//
//        const auto Warray
//            = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatialall), transformWith(Jinv));
//        const auto Warray2
//            = localF.evaluateDerivative(gpIndex, wrt(spatialall, coeff(i)), transformWith(Jinv));
//        for (int j = 0; j < 2; ++j)
//          EXPECT_THAT(Warray[j], EigenApproxEqual(Warray2[j], 1e-15));
//
//        const auto W0
//            = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(0)), transformWith(Jinv));
//        const auto W1
//            = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(1)), transformWith(Jinv));
//
//        EXPECT_THAT(Warray[0], EigenApproxEqual(W0, 1e-15));
//        EXPECT_THAT(Warray[1], EigenApproxEqual(W1, 1e-15));
//
//        const auto Warrayun  = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatialall));
//        const auto Warray2un = localF.evaluateDerivative(gpIndex, wrt(spatialall, coeff(i)));
//        for (int j = 0; j < 2; ++j)
//          EXPECT_THAT(Warrayun[j], EigenApproxEqual(Warray2un[j], 1e-15));
//
//        const auto W0un = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(0)));
//        const auto W1un = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(1)));
//
//        EXPECT_THAT(Warrayun[0], EigenApproxEqual(W0un, 1e-15));
//        EXPECT_THAT(Warrayun[1], EigenApproxEqual(W1un, 1e-15));
//        const auto Sun = localF.evaluateDerivative(gpIndex, wrt(coeff(i,i)), along(testVec));
//
//        const auto S = localF.evaluateDerivative(gpIndex, wrt(coeff(i,i)), along(testVec), transformWith(Jinv));
//
//        const auto chi = localF.evaluateDerivative(gpIndex, wrt(spatialall, coeff(i,i)), along(testVec),
//                                                   transformWith(Jinv));
//
//        const auto chi0 = localF.evaluateDerivative(gpIndex, wrt(spatial(0), coeff(i,i)), along(testVec),
//                                                    transformWith(Jinv));
//
//        const auto chi1 = localF.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(i,i)), along(testVec),
//                                                    transformWith(Jinv));
//
//        EXPECT_THAT(chi[0], EigenApproxEqual(chi0, 1e-15));
//        EXPECT_THAT(chi[1], EigenApproxEqual(chi1, 1e-15));
//      }
//
//      EXPECT_DOUBLE_EQ(directorCached.norm(), 1.0);
//      EXPECT_DOUBLE_EQ(directoreval.getValue().norm(), 1.0);
//      EXPECT_THAT(directorCached, EigenApproxEqual(directoreval.getValue(), 1e-15));
//      EXPECT_NEAR((directoreval.getValue().transpose() * jaco2).norm(), 0.0, 1e-15);
//      EXPECT_NEAR((Ikarus::UnitVector<double, size>::derivativeOfProjectionWRTposition(directoreval.getValue())
//                   * directoreval.getValue())
//                      .norm(),
//                  0.0, 1e-15);
//
//      ++gpIndex;
//    }
//  }
//}

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
      const auto& directorCached                         = localF.evaluateFunction(gpIndex).getValue();
      const Eigen::Vector<double, size> directorEmbedded = vasMat * localBasis.evaluateFunction(gpIndex);
      const auto& directoreval                           = localF.evaluateFunction(gp.position());

      const auto J     = toEigenMatrix(ele.geometry().jacobianTransposed(gp.position())).transpose().eval();
      const auto Jinv  = J.inverse().eval();
      const auto jaco2 = localF.evaluateDerivative(gpIndex, wrt(spatialall), transformWith(Jinv));

      const auto jaco2e     = localF.evaluateDerivative(gp.position(), wrt(spatialall), transformWith(Jinv));
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
      const auto jaco2un     = localF.evaluateDerivative(gpIndex, wrt(spatialall));
      const auto jaco2col0un = localF.evaluateDerivative(gpIndex, wrt(spatial(0)));
      const auto jaco2col1un = localF.evaluateDerivative(gpIndex, wrt(spatial(1)));
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

      const Eigen::Vector<double, size> testVec = Eigen::Vector<double, size>::UnitX();
      for (size_t i = 0; i < fe.size(); ++i) {
        const auto jacobianWRTCoeffs = localF.evaluateDerivative(gpIndex, wrt(coeff(i)));
        EXPECT_THAT(jacobianWRTCoeffs, EigenApproxEqual(Jdual.block<size, size>(0, i * size), 1e-15));

        const auto Warray
            = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatialall), transformWith(Jinv));
        const auto Warray2
            = localF.evaluateDerivative(gpIndex, wrt(spatialall, coeff(i)), transformWith(Jinv));
        for (int j = 0; j < 2; ++j)
          EXPECT_THAT(Warray[j], EigenApproxEqual(Warray2[j], 1e-15));

        const auto W0
            = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(0)), transformWith(Jinv));
        const auto W1
            = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatial(1)), transformWith(Jinv));

        EXPECT_THAT(Warray[0], EigenApproxEqual(W0, 1e-15));
        EXPECT_THAT(Warray[1], EigenApproxEqual(W1, 1e-15));

        const auto Warrayun  = localF.evaluateDerivative(gpIndex, wrt(coeff(i), spatialall));
        const auto Warray2un = localF.evaluateDerivative(gpIndex, wrt(spatialall, coeff(i)));
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










#include <ikarus/localFunctions/expressions/binaryExpr.hh>




template <typename E1, typename E2>
 class LocalFunctionSum : public Ikarus::BinaryLocalFunctionExpression<LocalFunctionSum<E1,E2>,E1,E2> {

 public:
  using  Base= Ikarus::BinaryLocalFunctionExpression<LocalFunctionSum<E1,E2>,E1,E2>;
  using  Base::BinaryLocalFunctionExpression;
  using Traits = Ikarus::LocalFunctionTraits<LocalFunctionSum>;



  /** \brief Type used for coordinates */
  using ctype = typename Traits::ctype;
  //    /** \brief Dimension of the coeffs */
  static constexpr int valueSize = Traits::valueSize;
  /** \brief Dimension of the grid */
  static constexpr int gridDim = Traits::gridDim;
  /** \brief Type for coordinate vector in world space */
  using FunctionReturnType = typename Traits::FunctionReturnType;
  /** \brief Type for the directional derivatives */
  using AlongType = typename Traits::AlongType;
  /** \brief Type for the coordinates to store the return value */
  using GlobalE = typename FunctionReturnType::CoordinateType;
  /** \brief Type for the Jacobian matrix */
  using Jacobian = typename Traits::Jacobian;
  /** \brief Type for a column of the Jacobian matrix */
  using JacobianColType = typename Traits::JacobianColType;
  /** \brief Type for the derivatives wrT the coeffiecients */
  using CoeffDerivMatrix = typename Traits::CoeffDerivMatrix;
  /** \brief Type for ansatz function values */
  using AnsatzFunctionType = typename Traits::AnsatzFunctionType;
  /** \brief Type for the Jacobian of the ansatz function values */
  using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;

  template<typename DomainTypeOrIntegrationPointIndex,typename... TransformArgs>
  auto evaluateFunctionImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const Ikarus::TransformWith<TransformArgs...>& transArgs) const {
    return eval(this->_u.evaluateFunction(ipIndexOrPosition,transArgs)+this->_v.evaluateFunction(ipIndexOrPosition,transArgs));
  }


   template<int DerivativeOrder,typename LocalFunctionEvaluationArgs_>
   auto evaluateDerivativeOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {
      return eval(evaluateDerivativeImpl(this->_u,localFunctionArgs) + evaluateDerivativeImpl(this->_v,localFunctionArgs));  }
};

template <typename E1, typename E2>
struct Ikarus::LocalFunctionTraits<LocalFunctionSum<E1,E2>> : public Ikarus::LocalFunctionTraits<E1>{
using Base =  Ikarus::LocalFunctionTraits<E1>;

};

template <typename E1, typename E2>
LocalFunctionSum<E1, E2>
operator+(Ikarus::LocalFunctionInterface<E1> const& u, Ikarus::LocalFunctionInterface<E2> const& v) {
  static_assert(Ikarus::Concepts::AddAble<typename E1::FunctionReturnType,typename E2::FunctionReturnType>, "The function values of your local functions are not addable!");

  return LocalFunctionSum<E1,E2>(u,v);
}


template <typename E1, typename E2>
class LocalFunctionDot : public Ikarus::BinaryLocalFunctionExpression<LocalFunctionDot<E1,E2>,E1,E2> {

 public:
  using  Base= Ikarus::BinaryLocalFunctionExpression<LocalFunctionDot<E1,E2>,E1,E2>;
  using  Base::BinaryLocalFunctionExpression;
  using Traits = Ikarus::LocalFunctionTraits<LocalFunctionDot>;

  /** \brief Type used for coordinates */
  using ctype = typename Traits::ctype;
  //    /** \brief Dimension of the coeffs */
  static constexpr int valueSize = Traits::valueSize;
  /** \brief Dimension of the grid */
  static constexpr int gridDim = Traits::gridDim;
  /** \brief Type for coordinate vector in world space */
  using FunctionReturnType = typename Traits::FunctionReturnType;
  /** \brief Type for the directional derivatives */
  using AlongType = typename Traits::AlongType;
  /** \brief Type for the coordinates to store the return value */
  using GlobalE = typename FunctionReturnType::CoordinateType;
  /** \brief Type for the Jacobian matrix */
  using Jacobian = typename Traits::Jacobian;
  /** \brief Type for a column of the Jacobian matrix */
  using JacobianColType = typename Traits::JacobianColType;
  /** \brief Type for the derivatives wrT the coeffiecients */
  using CoeffDerivMatrix = typename Traits::CoeffDerivMatrix;
  /** \brief Type for ansatz function values */
  using AnsatzFunctionType = typename Traits::AnsatzFunctionType;
  /** \brief Type for the Jacobian of the ansatz function values */
  using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;

  template<typename DomainTypeOrIntegrationPointIndex,typename... TransformArgs>
  auto evaluateFunctionImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const Ikarus::TransformWith<TransformArgs...>& transArgs) const {
    return eval(this->_u.evaluateFunction(ipIndexOrPosition,transArgs).dot(this->_v.evaluateFunction(ipIndexOrPosition,transArgs)));
  }


  template<int DerivativeOrder,typename LocalFunctionEvaluationArgs_>
  auto evaluateDerivativeOfExpression(const LocalFunctionEvaluationArgs_& localFunctionArgs) const {


    return eval(evaluateDerivativeImpl(this->_u,localFunctionArgs) *this->_v.evaluateFunctionImpl(localFunctionArgs.ipIndexOrPosition, localFunctionArgs.transArgs) +  this->_u.evaluateFunctionImpl(localFunctionArgs.ipIndexOrPosition, localFunctionArgs.transArgs) * evaluateDerivativeImpl(this->_v,localFunctionArgs) );  }
};

template <typename E1, typename E2>
struct Ikarus::LocalFunctionTraits<LocalFunctionDot<E1,E2>> {
  /** \brief Type used for coordinates */
  using ctype = typename E1::ctype;
  static constexpr int valueSize = 1;
  /** \brief Dimension of the grid */
  static constexpr int gridDim = E1::gridDim;
  /** \brief Type for the return value */
  using FunctionReturnType = ctype;
  /** \brief Type for the Jacobian matrix */
  using Jacobian = Eigen::Matrix<ctype, valueSize, gridDim>;
  /** \brief Type for the derivatives wrt. the coeffiecients */
  using CoeffDerivMatrix = Eigen::DiagonalMatrix<ctype, valueSize>;

  /** \brief Type for the points for evaluation, usually the integration points */
  using DomainType = typename E1::DomainType;
  /** \brief Type for a column of the Jacobian matrix */
  using JacobianColType = typename Eigen::internal::plain_col_type<Jacobian>::type;
  /** \brief Type for the directional derivatives */
  using AlongType = Eigen::Vector<ctype, valueSize>;

};



template <typename E1, typename E2>
LocalFunctionDot<E1, E2>
dot(Ikarus::LocalFunctionInterface<E1> const& u, Ikarus::LocalFunctionInterface<E2> const& v) {

  return LocalFunctionDot<E1,E2>(u,v);
}



TEST(LocalFunctionTests,TestExpressions) {
  using Manifold = Ikarus::RealTuple<double, 2>;
  constexpr int size = Manifold::valueSize;
  using namespace Ikarus;
  auto grid = createGrid<Grids::Yasp>();

  auto gridView = grid->leafGridView();
  auto basis = makeBasis(gridView, power<3>(lagrange<2>(), BlockedInterleaved()));
  Dune::BlockVector<Manifold> vBlocked(basis.size());
  Dune::BlockVector<Manifold> vBlocked2(basis.size());
  for (auto &vsingle: vBlocked)
    vsingle.setValue(0.1*Eigen::Vector<double, size>::Random() + Eigen::Vector<double, size>::UnitX());

  for (auto &vsingle: vBlocked2)
    vsingle.setValue(0.5*Eigen::Vector<double, size>::Random() + 7*Eigen::Vector<double, size>::UnitX());

  auto localView = basis.localView();
  for (auto &ele: elements(gridView)) {
    localView.bind(ele);
    const auto &fe = localView.tree().child(0).finiteElement();
    auto localBasis = Ikarus::LocalBasis(localView.tree().child(0).finiteElement().localBasis());
    const auto &rule = Dune::QuadratureRules<double, 2>::rule(localView.element().type(), 3);
    localBasis.bind(rule, bindDerivatives(0, 1));
    Dune::BlockVector<Manifold> vBlockedLocal(fe.size());
    Dune::BlockVector<Manifold> vBlockedLocal2(fe.size());
    Dune::BlockVector<typename Manifold::template Rebind<autodiff::dual>::type> vBlockedLocalDual(fe.size());
    Dune::BlockVector<typename Manifold::template Rebind<autodiff::dual2nd>::type> vBlockedLocalDual2nd(fe.size());
    for (size_t i = 0; i < fe.size(); ++i) {
      auto globalIndex = localView.index(localView.tree().child(0).localIndex(i));
      vBlockedLocal[i] = vBlocked[globalIndex[0]];
      vBlockedLocal2[i] = vBlocked2[globalIndex[0]];
    }
    auto f = Ikarus::StandardLocalFunction(localBasis, vBlockedLocal);
    auto g = Ikarus::StandardLocalFunction(localBasis, vBlockedLocal2);
    using namespace Ikarus::DerivativeDirections;
    auto h = f + g;

    for (int gpIndex = 0; auto &gp: rule) {
      EXPECT_THAT(f.evaluateFunction(gpIndex).getValue()+g.evaluateFunction(gpIndex).getValue(), EigenApproxEqual(h.evaluateFunction(gpIndex).getValue(), 1e-15));
      EXPECT_THAT(f.evaluateDerivative(gpIndex, wrt(spatial(0)))+g.evaluateDerivative(gpIndex, wrt(spatial(0))), EigenApproxEqual(h.evaluateDerivative(gpIndex, wrt(spatial(0))), 1e-15));
      EXPECT_THAT(f.evaluateDerivative(gpIndex, wrt(spatialall))+g.evaluateDerivative(gpIndex, wrt(spatialall)), EigenApproxEqual(h.evaluateDerivative(gpIndex, wrt(spatialall)), 1e-15));
      for (size_t i = 0; i < fe.size(); ++i) {
        const Eigen::Matrix2d dfgdi= f.evaluateDerivative(gpIndex, wrt(coeff(i)))+g.evaluateDerivative(gpIndex, wrt(coeff(i)));
        const Eigen::Matrix2d dhdi= h.evaluateDerivative(gpIndex, wrt(coeff(i)));
        const Eigen::Matrix2d dhdSdi= h.evaluateDerivative(gpIndex, wrt(spatial(1),coeff(i)));
        const Eigen::Matrix2d dfgdSdi= f.evaluateDerivative(gpIndex, wrt(spatial(1),coeff(i)))+g.evaluateDerivative(gpIndex, wrt(spatial(1),coeff(i)));
        EXPECT_THAT(dfgdi, EigenApproxEqual(dhdi, 1e-14));
        EXPECT_THAT(dhdSdi, EigenApproxEqual(dfgdSdi, 1e-14));
        for (size_t j = 0; j < fe.size(); ++j) {
          const Eigen::Matrix2d dfgdj= f.evaluateDerivative(gpIndex, wrt(coeff(j)))+g.evaluateDerivative(gpIndex, wrt(coeff(j)));
          const Eigen::Matrix2d dhdj= h.evaluateDerivative(gpIndex, wrt(coeff(j)));
          const Eigen::Matrix2d dhdSdj= h.evaluateDerivative(gpIndex, wrt(spatial(1),coeff(j)));
          const Eigen::Matrix2d dfgdSdj= f.evaluateDerivative(gpIndex, wrt(spatial(1),coeff(j)))+g.evaluateDerivative(gpIndex, wrt(spatial(1),coeff(j)));
          if (i==j) {
            EXPECT_THAT(dfgdj, EigenApproxEqual(dhdj, 1e-14));
            EXPECT_THAT(dfgdSdj, EigenApproxEqual(dhdSdj, 1e-14));

          }
        }
      }
    }

//    auto k = dot(f,g);
  }
}