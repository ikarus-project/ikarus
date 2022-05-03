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

template <typename LF,std::size_t... I>
auto& getCoeffRefHelper(LF& lf)
{
  if constexpr(LF::isLeaf)
    return lf.coefficientsRef();
  else
    return collectLeafNodeLocalFunctions(lf).coefficientsRef();
}




template <typename LF,bool isCopy=false>
void testLocalFunction(const LF& lf, int ipIndex) {
  spdlog::info("Testing: " + std::string(isCopy==true? "Copy ": "") + Ikarus::localFunctionName(lf));
  const double tol = 1e-13;
  using namespace Ikarus::DerivativeDirections;
  using namespace autodiff;
  using namespace Ikarus;
  const auto& coeffs  = getCoeffRefHelper(lf);
  const size_t coeffSize = coeffs.size();

  const auto lfCopy = lf.clone();
  const auto& coeffCopy = getCoeffRefHelper(lfCopy);
  for (size_t i = 0; i < coeffSize; ++i)
    EXPECT_EQ(coeffs[i],coeffCopy[i]); //since the coeffs are copied the values should be the same

  EXPECT_NE(&coeffs,&coeffCopy); //since the coeffs are copied the adresses should differ

  constexpr int gridDim = LF::gridDim;
  using Manifold = typename std::remove_cvref_t<decltype(coeffs)>::value_type;
  constexpr int localFunctionValueSize = LF::valueSize;
  constexpr int coeffValueSize = Manifold::valueSize;
  constexpr int coeffCorrectionSize = Manifold::correctionSize;

  /// Check spatial derivatives
  /// Check spatial derivatives return sizes
  {
    const auto spatialAllDerivative = lf.evaluateDerivative(ipIndex, Ikarus::wrt(spatialAll));
    static_assert(spatialAllDerivative.cols()==gridDim);
    static_assert(spatialAllDerivative.rows()==localFunctionValueSize);
    const auto spatialSingleDerivative = lf.evaluateDerivative(ipIndex, Ikarus::wrt(spatial(0)));
    static_assert(spatialSingleDerivative.cols()==1);
    static_assert(spatialSingleDerivative.rows()==localFunctionValueSize);

    /// Check if spatial derivatives are really derivatives
    /// Perturb in a random direction in the elements parameter space and check spatial derivative
    auto func = [&](auto &gpOffset_) { return lf.evaluateFunction(toFieldVector(gpOffset_)).getValue(); };
    auto spatialDerivAll =
        [&](auto &gpOffset_) { return lf.evaluateDerivative(toFieldVector(gpOffset_), Ikarus::wrt(spatialAll)); };
    auto derivDerivSingle = [&](auto gpOffset_, int spatialIndex) {
      return lf.evaluateDerivative(toFieldVector(gpOffset_),
                                   Ikarus::wrt(spatial(spatialIndex)));
    };
    Eigen::Vector<double, gridDim> ipOffset = (Eigen::Vector<double, gridDim>::Random()).normalized();
    auto nonLinOpSpatialAll =
        Ikarus::NonLinearOperator(linearAlgebraFunctions(func, spatialDerivAll), parameter(ipOffset));
    EXPECT_TRUE((checkJacobian<decltype(nonLinOpSpatialAll), Eigen::Vector<double, gridDim>>(nonLinOpSpatialAll,
                                                                                             false)));

    /// Perturb each spatial direction and check with derivative value
    for (int i = 0; i < gridDim; ++i) {
      Eigen::Vector<double, 1> ipOffsetSingle ( ipOffset[i]);
      auto derivDerivSingleI = [&](auto gpOffset_) {
        auto offSetSingle = ipOffset;
        offSetSingle[i]+=gpOffset_[0];
        return lf.evaluateDerivative(toFieldVector(offSetSingle),
                                     Ikarus::wrt(spatial(i)));
      };

      auto funcSingle = [&](auto &gpOffset_) {
        auto offSetSingle = ipOffset;
        offSetSingle[i]+=gpOffset_[0];
        return lf.evaluateFunction(toFieldVector(offSetSingle)).getValue(); };

      auto nonLinOpSpatialSingle =
          Ikarus::NonLinearOperator(linearAlgebraFunctions(funcSingle, derivDerivSingleI), parameter(ipOffsetSingle));
      EXPECT_TRUE((checkJacobian<decltype(nonLinOpSpatialSingle), Eigen::Vector<double, 1>>(nonLinOpSpatialSingle,
                                                                                            false)));
    }
  }

  /// Check coeff and spatial derivatives
  Eigen::VectorXdual xv(correctionSize(coeffs));
  xv.setZero();

  const Eigen::Vector<double, localFunctionValueSize> alongVec = localFunctionValueSize==1 ? Eigen::Vector<double, localFunctionValueSize>::Ones().eval() : Eigen::Vector<double, localFunctionValueSize>::Random().eval();
  const Eigen::Matrix<double, localFunctionValueSize, gridDim> alongMat = localFunctionValueSize==1 ? Eigen::Matrix<double, localFunctionValueSize, gridDim>::Ones().eval() : Eigen::Matrix<double, localFunctionValueSize, gridDim>::Random().eval();
  /// Rebind local function to second order dual number
  auto lfDual2nd = lf.rebindClone(dual2nd ());
  auto lfDual2ndLeafNodeCollection = collectLeafNodeLocalFunctions(lfDual2nd);

  auto localFdual2nd = [&]  (const auto& x)  {
    lfDual2ndLeafNodeCollection.addToCoeffs(x);
    auto value =  (lfDual2nd.evaluateFunction(ipIndex).getValue().transpose()*alongVec).trace();
    lfDual2ndLeafNodeCollection.addToCoeffs(-x);
    return value;
  };

  auto localFdual2ndSpatialSingle = [&](const auto& x,int i) {
    lfDual2ndLeafNodeCollection.addToCoeffs(x);
    auto value =  (lfDual2nd.evaluateDerivative(ipIndex,Ikarus::wrt(spatial(i))).transpose()*alongVec).trace();
    lfDual2ndLeafNodeCollection.addToCoeffs(-x);
    return value;
  };

  auto localFdual2ndSpatialAll = [&](const auto& x) {
    lfDual2ndLeafNodeCollection.addToCoeffs(x);
    auto value =  (lfDual2nd.evaluateDerivative(ipIndex,Ikarus::wrt(spatialAll)).transpose()*alongMat).trace();
    lfDual2ndLeafNodeCollection.addToCoeffs(-x);
    return value;
  };

  Eigen::VectorXdual2nd  xvr(correctionSize(coeffs));
  xvr.setZero();

  dual2nd u;
  Eigen::VectorXd gradienWRTCoeffs;
  Eigen::MatrixXd hessianWRTCoeffs;
  autodiff::hessian(localFdual2nd, autodiff::wrt(xvr), autodiff::at(xvr), u, gradienWRTCoeffs, hessianWRTCoeffs);

  Eigen::VectorXd gradienWRTCoeffsSpatialAll;
  Eigen::MatrixXd hessianWRTCoeffsSpatialAll;
  autodiff::hessian(localFdual2ndSpatialAll, autodiff::wrt(xvr), autodiff::at(xvr), u, gradienWRTCoeffsSpatialAll, hessianWRTCoeffsSpatialAll);

  std::array<Eigen::MatrixXd,gridDim>  hessianWRTCoeffsTwoTimesSingleSpatial;
  std::array<Eigen::VectorXd,gridDim> gradientWRTCoeffsTwoTimesSingleSpatial;
  for (int d = 0; d < gridDim; ++d) {
    autodiff::hessian(localFdual2ndSpatialSingle,
                      autodiff::wrt(xvr),
                      autodiff::at(xvr,d),
                      u,
                      gradientWRTCoeffsTwoTimesSingleSpatial[d],
                      hessianWRTCoeffsTwoTimesSingleSpatial[d]);
  }

  for (size_t i = 0; i < coeffSize; ++i) {
    const auto jacobianWRTCoeffslf= lf.evaluateDerivative(ipIndex, Ikarus::wrt(coeff(i)));
    static_assert(jacobianWRTCoeffslf.ColsAtCompileTime == coeffCorrectionSize );
    static_assert(jacobianWRTCoeffslf.RowsAtCompileTime == localFunctionValueSize);
    const auto jacobianWRTCoeffs = ((alongVec.transpose()*jacobianWRTCoeffslf).transpose()).eval();
    static_assert(jacobianWRTCoeffs.cols() == 1);
    static_assert(jacobianWRTCoeffs.rows() == coeffCorrectionSize);
    EXPECT_THAT(
        jacobianWRTCoeffs,
        EigenApproxEqual(gradienWRTCoeffs.template segment<coeffCorrectionSize>(i*coeffCorrectionSize), tol));

    for (int d = 0; d < gridDim; ++d) {
      const auto jacoWrtCoeffAndSpatiallf= lf.evaluateDerivative(ipIndex, Ikarus::wrt(coeff(i), spatial(d)));
      static_assert(jacoWrtCoeffAndSpatiallf.ColsAtCompileTime == coeffCorrectionSize);
      static_assert(jacoWrtCoeffAndSpatiallf.RowsAtCompileTime == localFunctionValueSize );
      const auto jacoWrtCoeffAndSpatial = ((alongVec.transpose()*jacoWrtCoeffAndSpatiallf).transpose()).eval();
      const auto jacoWrtSpatialAndCoeff = ((alongVec.transpose()*lf.evaluateDerivative(ipIndex, Ikarus::wrt( spatial(d),coeff(i)))).transpose()).eval();

      static_assert(jacoWrtSpatialAndCoeff.cols() == 1);
      static_assert(jacoWrtSpatialAndCoeff.rows() == coeffCorrectionSize);
      static_assert(jacoWrtCoeffAndSpatial.rows() == jacoWrtCoeffAndSpatial.rows()  and jacoWrtCoeffAndSpatial.cols() == jacoWrtCoeffAndSpatial.cols());

      EXPECT_THAT(jacoWrtCoeffAndSpatial, EigenApproxEqual(jacoWrtSpatialAndCoeff, tol));

      EXPECT_THAT(jacoWrtCoeffAndSpatial,
                  EigenApproxEqual(gradientWRTCoeffsTwoTimesSingleSpatial[d].template segment<coeffCorrectionSize>(i
                                                                                              *coeffCorrectionSize
                                                                                                                ),
                                   tol));
    }

    const auto jacoWrtSpatialAllAndCoeff = lf.evaluateDerivative(ipIndex, Ikarus::wrt( spatialAll,coeff(i)));

    static_assert(jacoWrtSpatialAllAndCoeff[0].ColsAtCompileTime  == coeffCorrectionSize);
    static_assert(jacoWrtSpatialAllAndCoeff[0].RowsAtCompileTime  == localFunctionValueSize );

    Eigen::Matrix<double,1,coeffCorrectionSize> jacoWrtSpatialAllAndCoeffProd;
    jacoWrtSpatialAllAndCoeffProd.setZero();
//    std::cout<<"alongMat"<<std::endl;
//    std::cout<<alongMat<<std::endl;

    for (int d = 0; d < gridDim; ++d)
        jacoWrtSpatialAllAndCoeffProd+= (alongMat.col(d).transpose()*jacoWrtSpatialAllAndCoeff[d]).eval();

    EXPECT_THAT(jacoWrtSpatialAllAndCoeffProd,
                EigenApproxEqual(gradienWRTCoeffsSpatialAll.template segment<coeffCorrectionSize>(
                                     i*coeffCorrectionSize
                                 ).transpose(),
                                 tol));

    for (size_t j = 0; j < coeffSize; ++j) {
      const auto
          jacobianWRTCoeffsTwoTimes = lf.evaluateDerivative(ipIndex, Ikarus::wrt(coeff(i, j)), Ikarus::along(alongVec));
      static_assert(jacobianWRTCoeffsTwoTimes.cols() == coeffCorrectionSize );
      static_assert(jacobianWRTCoeffsTwoTimes.rows() == coeffCorrectionSize);
      const auto jacobianWRTCoeffsTwoTimesExpected =
          hessianWRTCoeffs.template block<coeffCorrectionSize, coeffCorrectionSize>(i*coeffCorrectionSize,
                                                                                    j*coeffCorrectionSize);
      EXPECT_THAT(
          jacobianWRTCoeffsTwoTimes,
          EigenApproxEqual(jacobianWRTCoeffsTwoTimesExpected, tol));

      const auto
          jacobianWRTCoeffsTwoTimesSpatialAll = lf.evaluateDerivative(ipIndex, Ikarus::wrt(coeff(i, j),spatialAll), Ikarus::along(alongMat));
      static_assert(jacobianWRTCoeffsTwoTimesSpatialAll.cols() == coeffCorrectionSize );
      static_assert(jacobianWRTCoeffsTwoTimesSpatialAll.rows() == coeffCorrectionSize);
      const auto jacobianWRTCoeffsTwoTimesSpatialAllExpected =
          hessianWRTCoeffsSpatialAll.template block<coeffCorrectionSize, coeffCorrectionSize>(i*coeffCorrectionSize,
                                                                                    j*coeffCorrectionSize);

      /// if the order of the function value is less then quadratic then this should yield a vanishing derivative
      if constexpr(lf.template order<> <quadratic) {
        EXPECT_TRUE(jacobianWRTCoeffsTwoTimesSpatialAll.norm()<tol);
        EXPECT_TRUE(jacobianWRTCoeffsTwoTimesSpatialAllExpected.norm()<tol);
      } else
      {
//        std::cout<<Dune::className(lf)<<std::endl;
        EXPECT_THAT(
            jacobianWRTCoeffsTwoTimesSpatialAll,
            EigenApproxEqual(jacobianWRTCoeffsTwoTimesSpatialAllExpected, tol));
      }

      for (int d = 0; d < gridDim; ++d) {
        const auto jacobianWRTCoeffsTwoTimesSingleSpatial =
            lf.evaluateDerivative(ipIndex, Ikarus::wrt(coeff(i, j), spatial(d)), Ikarus::along(alongVec));
        static_assert(jacobianWRTCoeffsTwoTimesSingleSpatial.cols() == coeffCorrectionSize );
        static_assert(jacobianWRTCoeffsTwoTimesSingleSpatial.rows() == coeffCorrectionSize);
        const auto jacobianWRTCoeffsTwoTimesSingleSpatialExpected =
            hessianWRTCoeffsTwoTimesSingleSpatial[d].template block<coeffCorrectionSize, coeffCorrectionSize>(
                i*coeffCorrectionSize, j*coeffCorrectionSize);
        EXPECT_THAT(
            jacobianWRTCoeffsTwoTimesSingleSpatial,
            EigenApproxEqual(jacobianWRTCoeffsTwoTimesSingleSpatialExpected, tol));
      }
    }
  }
  spdlog::info("done. ");
  if constexpr(not isCopy)
    testLocalFunction<std::remove_cvref_t<decltype(lfCopy)>,true>(lfCopy,ipIndex); //also test the cloned local function
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
  constexpr int gridDim =std::remove_reference_t<decltype(*grid)>::dimension;

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
    static_assert(f.template order<> == linear);
    auto g = Ikarus::StandardLocalFunction(localBasis, vBlockedLocal);
    static_assert(g.template order<> == linear);

    static_assert(countNonArithmeticLeafNodes(f) == 1);
    static_assert(countNonArithmeticLeafNodes(g) == 1);
    using namespace Ikarus::DerivativeDirections;
    auto h = f + g;
    static_assert(h.template order<> == linear);

    auto a     = collectNonArithmeticLeafNodes(h);
    auto hLeaf = collectLeafNodeLocalFunctions(h);
    static_assert(std::tuple_size_v<decltype(a)> == 2);

    static_assert(countNonArithmeticLeafNodes(h) == 2);
    static_assert(std::is_same_v<decltype(h)::Ids, std::tuple<Dune::index_constant<0>, Dune::index_constant<0>>>);

    auto mf = -f;
    static_assert(f.template order<> == mf.template order<>);

    for (size_t i = 0; i < fe.size(); ++i) {
      EXPECT_TRUE(collectLeafNodeLocalFunctions(h).coefficientsRef(_0)[i] == vBlockedLocal[i]);
      EXPECT_TRUE(collectLeafNodeLocalFunctions(h).coefficientsRef(_1)[i] == vBlockedLocal[i]);
    }

    for (int gpIndex = 0; auto& gp : rule) {
      testLocalFunction(f,gpIndex);
      testLocalFunction(2*f,gpIndex);
      testLocalFunction(2*f*3,gpIndex);
      testLocalFunction(mf,gpIndex);
      testLocalFunction(h,gpIndex);
//      testLocalFunction(h,gpIndex);

      EXPECT_THAT(f.evaluateFunction(gpIndex).getValue() + g.evaluateFunction(gpIndex).getValue(),
                  EigenApproxEqual(h.evaluateFunction(gpIndex).getValue(), 1e-15));

      for (size_t i = 0; i < fe.size(); ++i) {
        const auto dhdSAdi  = h.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i)));
        const auto dfgdSAdi = (f.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i)))
            + g.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i))));
        for (int d = 0; d < gridDim; ++d)
          EXPECT_THAT(dfgdSAdi[d], EigenApproxEqual(dhdSAdi[d], 1e-14));

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
    static_assert(k.template order<> == quadratic);
    auto b = collectNonArithmeticLeafNodes(k);
    static_assert(std::tuple_size_v<decltype(b)> == 3);

    //    std::cout<<Dune::className(a)<<std::endl;
//    std::cout << Dune::className(b) << std::endl;
    static_assert(countNonArithmeticLeafNodes(k) == 3);
    static_assert(std::is_same_v<decltype(k)::Ids, std::tuple<Dune::index_constant<0>, Dune::index_constant<0>,
                                                              Ikarus::Arithmetic, Dune::index_constant<0>>>);

    const double tol = 1e-13;

    auto dotff = dot(f,g);

    static_assert(countNonArithmeticLeafNodes(dotff) == 2);
    static_assert(dotff.template order<> == quadratic);
    static_assert(std::is_same_v<decltype(dotff)::Ids, std::tuple<Dune::index_constant<0>, Dune::index_constant<0>     >>);
    for (int gpIndex = 0; auto& gp : rule) {
      testLocalFunction(dotff,gpIndex);
      testLocalFunction(k,gpIndex);


      const auto& N  = localBasis.evaluateFunction(gpIndex);
      const auto& dN = localBasis.evaluateJacobian(gpIndex);
      EXPECT_DOUBLE_EQ((-2 * 3) * f.evaluateFunction(gpIndex).getValue().dot(g.evaluateFunction(gpIndex).getValue()),
                       k.evaluateFunction(gpIndex).getValue()[0]);
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
    static_assert(countNonArithmeticLeafNodes(f2) == 1);
    static_assert(countNonArithmeticLeafNodes(g2) == 1);

    auto k2 = dot(f2 + g2, g2);
    static_assert(countNonArithmeticLeafNodes(k2) == 3);
    static_assert(
        std::is_same_v<decltype(k2)::Ids,
                       std::tuple<Dune::index_constant<0>, Dune::index_constant<1>, Dune::index_constant<1>>>);

    auto b2 = collectNonArithmeticLeafNodes(k2);
    static_assert(std::tuple_size_v<decltype(b2)> == 3);

    for (int gpIndex = 0; auto& gp : rule) {
//      testLocalFunction(k2,gpIndex);
      const auto& N  = localBasis.evaluateFunction(gpIndex);
      const auto& dN = localBasis.evaluateJacobian(gpIndex);
      EXPECT_DOUBLE_EQ((f2.evaluateFunction(gpIndex).getValue() + g2.evaluateFunction(gpIndex).getValue())
                           .dot(g2.evaluateFunction(gpIndex).getValue()),
                       k2.evaluateFunction(gpIndex).getValue()[0]);
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
      const Eigen::Matrix<double, size, size> testMat= Eigen::Matrix<double, size, size>::Random();
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
            = localF.evaluateDerivative(gpIndex, wrt(spatialAll, coeff(i, i)), along(testMat), transformWith(Jinv));

        const auto chi0
            = localF.evaluateDerivative(gpIndex, wrt(spatial(0), coeff(i, i)), along(testMat.col(0)), transformWith(Jinv));

        const auto chi1
            = localF.evaluateDerivative(gpIndex, wrt(spatial(1), coeff(i, i)), along(testMat.col(1)), transformWith(Jinv));

        EXPECT_THAT(chi, EigenApproxEqual(chi0+chi1, 1e-15));

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
