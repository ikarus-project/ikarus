
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
#include <ikarus/utils/linearAlgebraHelper.hh>

using namespace Dune::Functions::BasisFactory;

template <typename LF, bool isCopy = false>
auto testLocalFunction(const LF &lf) {
  TestSuite t(Ikarus::localFunctionName(lf));
  spdlog::info("Testing: " + std::string(isCopy ? "Copy " : "") + Ikarus::localFunctionName(lf));

  const double tol = 1e-12;
  using namespace Ikarus::DerivativeDirections;
  using namespace autodiff;
  using namespace Ikarus;
  const auto &coeffs     = lf.coefficientsRef();
  const size_t coeffSize = coeffs.size();

  constexpr int gridDim                = LF::gridDim;
  using Manifold                       = typename std::remove_cvref_t<decltype(coeffs)>::value_type;
  constexpr int localFunctionValueSize = LF::Traits::valueSize;
  constexpr int coeffValueSize         = Manifold::valueSize;
  using ctype                          = typename Manifold::ctype;
  constexpr int coeffCorrectionSize    = Manifold::correctionSize;

  // dynamic sized vectors before the loop
  Eigen::VectorXdual2nd xvr(valueSize(coeffs));
  xvr.setZero();
  Eigen::VectorXd gradienWRTCoeffs;
  Eigen::MatrixXd hessianWRTCoeffs;
  Eigen::VectorXd gradienWRTCoeffsSpatialAll;
  Eigen::MatrixXd hessianWRTCoeffsSpatialAll;
  std::array<Eigen::MatrixXd, gridDim> hessianWRTCoeffsTwoTimesSingleSpatial;
  std::array<Eigen::VectorXd, gridDim> gradientWRTCoeffsTwoTimesSingleSpatial;
  for (const auto &[ipIndex, ip] : lf.viewOverIntegrationPoints()) {
    /// Check spatial derivatives
    /// Check spatial derivatives return sizes
    {
      if constexpr (requires { lf.evaluateDerivative(ipIndex, Ikarus::wrt(spatialAll)); }) {
        const decltype(lf.evaluateDerivative(ipIndex, Ikarus::wrt(spatialAll))) spatialAllDerivative;
        static_assert(spatialAllDerivative.cols() == gridDim);
        static_assert(spatialAllDerivative.rows() == localFunctionValueSize);
      }
      if constexpr (requires { lf.evaluateDerivative(ipIndex, Ikarus::wrt(spatial(0))); }) {
        const decltype(lf.evaluateDerivative(ipIndex, Ikarus::wrt(spatial(0)))) spatialSingleDerivative;
        static_assert(spatialSingleDerivative.cols() == 1);
        static_assert(spatialSingleDerivative.rows() == localFunctionValueSize);
      }

      /// Check if spatial derivatives are really derivatives
      /// Perturb in a random direction in the elements parameter space and check spatial derivative
      auto func = [&](auto &gpOffset_) { return lf.evaluateFunction(toDune(gpOffset_)); };
      auto spatialDerivAll
          = [&](auto &gpOffset_) { return lf.evaluateDerivative(toDune(gpOffset_), Ikarus::wrt(spatialAll)); };

      Eigen::Vector<double, gridDim> ipOffset = (Eigen::Vector<double, gridDim>::Random()).normalized() / 16;
      try {
        auto nonLinOpSpatialAll
            = Ikarus::NonLinearOperator(linearAlgebraFunctions(func, spatialDerivAll), parameter(ipOffset));

        t.check((checkJacobian<decltype(nonLinOpSpatialAll), Eigen::Vector<double, gridDim>>(
            nonLinOpSpatialAll, {.draw = false, .writeSlopeStatementIfFailed = true, .tolerance = 1e-2})));
      } catch (const Dune::NotImplemented &exception) {
        spdlog::info(
            "SpatialDerivative in all directions not tested, since it is not implemented by the local function "
            "(expression)");
      }

      /// Perturb each spatial direction and check with derivative value
      for (int i = 0; i < gridDim; ++i) {
        Eigen::Vector<double, 1> ipOffsetSingle(ipOffset[i]);
        auto derivDerivSingleI = [&](auto gpOffset_) {
          auto offSetSingle = ipOffset;
          offSetSingle[i] += gpOffset_[0];
          return lf.evaluateDerivative(toDune(offSetSingle), Ikarus::wrt(spatial(i)));
        };

        auto funcSingle = [&](const auto &gpOffset_) {
          auto offSetSingle = ipOffset;
          offSetSingle[i] += gpOffset_[0];
          return Eigen::Vector<ctype, localFunctionValueSize>(lf.evaluateFunction(toDune(offSetSingle)));
        };

        try {
          auto nonLinOpSpatialSingle = Ikarus::NonLinearOperator(linearAlgebraFunctions(funcSingle, derivDerivSingleI),
                                                                 parameter(ipOffsetSingle));
          t.check((checkJacobian<decltype(nonLinOpSpatialSingle), Eigen::Vector<double, 1>>(
              nonLinOpSpatialSingle, {.draw = false, .writeSlopeStatementIfFailed = true, .tolerance = 1e-2})));
        } catch (const Dune::NotImplemented &exception) {
          spdlog::info(
              "Single SpatialDerivative not tested, since it is not implemented by the local function (expression)");
        }
      }
    }

    /// Check coeff and spatial derivatives

    const Eigen::Vector<double, localFunctionValueSize> alongVec
        = localFunctionValueSize == 1 ? Eigen::Vector<double, localFunctionValueSize>::Ones().eval()
                                      : Eigen::Vector<double, localFunctionValueSize>::Random().eval();
    const Eigen::Matrix<double, localFunctionValueSize, gridDim> alongMat
        = localFunctionValueSize == 1 ? Eigen::Matrix<double, localFunctionValueSize, gridDim>::Ones().eval()
                                      : Eigen::Matrix<double, localFunctionValueSize, gridDim>::Random().eval();
    /// Rebind local function to second order dual number
    auto lfDual2nd                   = lf.rebindClone(dual2nd());
    auto lfDual2ndLeafNodeCollection = collectLeafNodeLocalFunctions(lfDual2nd);

    auto localFdual2nd = [&](const auto &x) {
      lfDual2ndLeafNodeCollection.addToCoeffsInEmbedding(x);
      auto value = (transpose(lfDual2nd.evaluateFunction(ipIndex)) * alongVec).trace();
      lfDual2ndLeafNodeCollection.addToCoeffsInEmbedding(-x);
      return value;
    };

    auto localFdual2ndSpatialSingle = [&](const auto &x, int i) {
      lfDual2ndLeafNodeCollection.addToCoeffsInEmbedding(x);
      auto value = (lfDual2nd.evaluateDerivative(ipIndex, Ikarus::wrt(spatial(i))).transpose() * alongVec).trace();
      lfDual2ndLeafNodeCollection.addToCoeffsInEmbedding(-x);
      return value;
    };

    auto localFdual2ndSpatialAll = [&](const auto &x) {
      lfDual2ndLeafNodeCollection.addToCoeffsInEmbedding(x);
      auto value = (lfDual2nd.evaluateDerivative(ipIndex, Ikarus::wrt(spatialAll)).transpose() * alongMat).trace();
      lfDual2ndLeafNodeCollection.addToCoeffsInEmbedding(-x);
      return value;
    };

    dual2nd u;

    autodiff::hessian(localFdual2nd, autodiff::wrt(xvr), autodiff::at(xvr), u, gradienWRTCoeffs, hessianWRTCoeffs);
    bool spatialSingleImplemented{true}, spatialAllImplemented{true};
    try {
      autodiff::hessian(localFdual2ndSpatialAll, autodiff::wrt(xvr), autodiff::at(xvr), u, gradienWRTCoeffsSpatialAll,
                        hessianWRTCoeffsSpatialAll);
    } catch (const Dune::NotImplemented &exception) {
      spdlog::info(
          "SpatialDerivative in all directions not tested, since it is not implemented by the local function "
          "(expression)");
      spatialSingleImplemented = false;
    }
    try {
      for (int d = 0; d < gridDim; ++d)
        autodiff::hessian(localFdual2ndSpatialSingle, autodiff::wrt(xvr), autodiff::at(xvr, d), u,
                          gradientWRTCoeffsTwoTimesSingleSpatial[d], hessianWRTCoeffsTwoTimesSingleSpatial[d]);

    } catch (const Dune::NotImplemented &exception) {
      spdlog::info(
          "Single SpatialDerivative not tested, since it is not implemented by the local function (expression)");
      spatialAllImplemented = false;
    }

    for (size_t i = 0; i < coeffSize; ++i) {
      const auto BLAi                = coeffs[i].orthonormalFrame();
      const auto jacobianWRTCoeffslf = lf.evaluateDerivative(ipIndex, Ikarus::wrt(coeff(i)));
      static_assert(jacobianWRTCoeffslf.ColsAtCompileTime == coeffCorrectionSize);
      static_assert(jacobianWRTCoeffslf.RowsAtCompileTime == localFunctionValueSize);
      const auto jacobianWRTCoeffs = ((alongVec.transpose() * jacobianWRTCoeffslf).transpose()).eval();
      static_assert(jacobianWRTCoeffs.cols() == 1);
      static_assert(jacobianWRTCoeffs.rows() == coeffCorrectionSize);
      t.check(isApproxSame(jacobianWRTCoeffs,
                           BLAi.transpose() * gradienWRTCoeffs.template segment<coeffValueSize>(i * coeffValueSize),
                           tol));

      if (spatialSingleImplemented)
        for (int d = 0; d < gridDim; ++d) {
          const auto jacoWrtCoeffAndSpatiallf = lf.evaluateDerivative(ipIndex, Ikarus::wrt(coeff(i), spatial(d)));
          static_assert(jacoWrtCoeffAndSpatiallf.ColsAtCompileTime == coeffCorrectionSize);
          static_assert(jacoWrtCoeffAndSpatiallf.RowsAtCompileTime == localFunctionValueSize);
          const auto jacoWrtCoeffAndSpatial = ((alongVec.transpose() * jacoWrtCoeffAndSpatiallf).transpose()).eval();
          const auto jacoWrtSpatialAndCoeff
              = ((alongVec.transpose() * lf.evaluateDerivative(ipIndex, Ikarus::wrt(spatial(d), coeff(i)))).transpose())
                    .eval();

          static_assert(jacoWrtSpatialAndCoeff.cols() == 1);
          static_assert(jacoWrtSpatialAndCoeff.rows() == coeffCorrectionSize);
          static_assert(jacoWrtCoeffAndSpatial.rows() == jacoWrtCoeffAndSpatial.rows()
                        and jacoWrtCoeffAndSpatial.cols() == jacoWrtCoeffAndSpatial.cols());

          t.check(isApproxSame(jacoWrtCoeffAndSpatial, jacoWrtSpatialAndCoeff, tol));

          t.check(isApproxSame(
              jacoWrtCoeffAndSpatial,
              BLAi.transpose()
                  * gradientWRTCoeffsTwoTimesSingleSpatial[d].template segment<coeffValueSize>(i * coeffValueSize),
              tol));
        }

      if (spatialAllImplemented and spatialSingleImplemented) {
        const auto jacoWrtSpatialAllAndCoeff = lf.evaluateDerivative(ipIndex, Ikarus::wrt(spatialAll, coeff(i)));

        static_assert(jacoWrtSpatialAllAndCoeff[0].ColsAtCompileTime == coeffCorrectionSize);
        static_assert(jacoWrtSpatialAllAndCoeff[0].RowsAtCompileTime == localFunctionValueSize);

        Eigen::Matrix<double, 1, coeffCorrectionSize> jacoWrtSpatialAllAndCoeffProd;
        jacoWrtSpatialAllAndCoeffProd.setZero();

        for (int d = 0; d < gridDim; ++d)
          jacoWrtSpatialAllAndCoeffProd += (alongMat.col(d).transpose() * jacoWrtSpatialAllAndCoeff[d]).eval();

        t.check(isApproxSame(
            jacoWrtSpatialAllAndCoeffProd,
            (BLAi.transpose() * gradienWRTCoeffsSpatialAll.template segment<coeffValueSize>(i * coeffValueSize))
                .transpose(),
            tol));

        // Check if spatialAll returns the same as the single spatial derivatives
        const auto Warray  = lf.evaluateDerivative(ipIndex, Ikarus::wrt(coeff(i), spatialAll));
        const auto Warray2 = lf.evaluateDerivative(ipIndex, Ikarus::wrt(spatialAll, coeff(i)));
        for (int j = 0; j < gridDim; ++j)
          t.check(isApproxSame(Warray[j], Warray2[j], tol));

        std::array<std::remove_cvref_t<decltype(Warray[0])>, gridDim> WarraySingle;
        for (int s = 0; s < gridDim; ++s)
          WarraySingle[s] = lf.evaluateDerivative(ipIndex, Ikarus::wrt(coeff(i), spatial(s)));

        for (int j = 0; j < gridDim; ++j)
          t.check(isApproxSame(Warray[j], WarraySingle[j], tol));
      }
      for (size_t j = 0; j < coeffSize; ++j) {
        const auto BLAj = coeffs[j].orthonormalFrame();
        const auto jacobianWRTCoeffsTwoTimes
            = lf.evaluateDerivative(ipIndex, Ikarus::wrt(coeff(i, j)), Ikarus::along(alongVec));
        static_assert(jacobianWRTCoeffsTwoTimes.cols() == coeffCorrectionSize);
        static_assert(jacobianWRTCoeffsTwoTimes.rows() == coeffCorrectionSize);
        const auto jacobianWRTCoeffsTwoTimesExpected
            = (BLAi.transpose()
                   * hessianWRTCoeffs.template block<coeffValueSize, coeffValueSize>(i * coeffValueSize,
                                                                                     j * coeffValueSize)
                   * BLAj
               + (i == j)
                     * coeffs[j].weingartenMap(gradienWRTCoeffs.template segment<coeffValueSize>(i * coeffValueSize)))
                  .eval();
        t.check(isApproxSame(jacobianWRTCoeffsTwoTimes, jacobianWRTCoeffsTwoTimesExpected, tol));

        if (spatialAllImplemented) {
          const auto jacobianWRTCoeffsTwoTimesSpatialAll
              = lf.evaluateDerivative(ipIndex, Ikarus::wrt(coeff(i, j), spatialAll), Ikarus::along(alongMat));
          static_assert(jacobianWRTCoeffsTwoTimesSpatialAll.cols() == coeffCorrectionSize);
          static_assert(jacobianWRTCoeffsTwoTimesSpatialAll.rows() == coeffCorrectionSize);
          const auto jacobianWRTCoeffsTwoTimesSpatialAllExpected
              = (BLAi.transpose()
                     * hessianWRTCoeffsSpatialAll.template block<coeffValueSize, coeffValueSize>(i * coeffValueSize,
                                                                                                 j * coeffValueSize)
                     * BLAj
                 + (i == j)
                       * coeffs[j].weingartenMap(
                           gradienWRTCoeffsSpatialAll.template segment<coeffValueSize>(i * coeffValueSize)))
                    .eval();

          /// if the order of the function value is less then quadratic then this should yield a vanishing derivative
          if constexpr (lf.order() < quadratic) {
            t.check(jacobianWRTCoeffsTwoTimesSpatialAll.norm() < tol);
            t.check(jacobianWRTCoeffsTwoTimesSpatialAllExpected.norm() < tol);
          } else {
            t.check(
                isApproxSame(jacobianWRTCoeffsTwoTimesSpatialAll, jacobianWRTCoeffsTwoTimesSpatialAllExpected, tol));
          }
        }
        if (spatialSingleImplemented)
          for (int d = 0; d < gridDim; ++d) {
            const auto jacobianWRTCoeffsTwoTimesSingleSpatial
                = lf.evaluateDerivative(ipIndex, Ikarus::wrt(coeff(i, j), spatial(d)), Ikarus::along(alongVec));
            static_assert(jacobianWRTCoeffsTwoTimesSingleSpatial.cols() == coeffCorrectionSize);
            static_assert(jacobianWRTCoeffsTwoTimesSingleSpatial.rows() == coeffCorrectionSize);
            const auto jacobianWRTCoeffsTwoTimesSingleSpatialExpected
                = (BLAi.transpose()
                       * hessianWRTCoeffsTwoTimesSingleSpatial[d].template block<coeffValueSize, coeffValueSize>(
                           i * coeffValueSize, j * coeffValueSize)
                       * BLAj
                   + (i == j)
                         * coeffs[j].weingartenMap(
                             gradientWRTCoeffsTwoTimesSingleSpatial[d].template segment<coeffValueSize>(
                                 i * coeffValueSize)))
                      .eval();
            t.check(isApproxSame(jacobianWRTCoeffsTwoTimesSingleSpatial, jacobianWRTCoeffsTwoTimesSingleSpatialExpected,
                                 tol));
          }
      }
    }
  }
  spdlog::info("done. ");
  if constexpr (not isCopy) {  //  test the cloned local function
    const auto lfCopy     = lf.clone();
    const auto &coeffCopy = lfCopy.coefficientsRef();
    for (size_t i = 0; i < coeffSize; ++i)
      t.check(coeffCopy[i] == coeffs[i]);

    t.check(&coeffCopy != &coeffs);

    t.subTest(testLocalFunction<std::remove_cvref_t<decltype(lfCopy)>, true>(lfCopy));
  }
  return t;
}

template <int domainDim, int worldDim, int order>
auto localFunctionTestConstructor(const Dune::GeometryType &geometryType, size_t nNodalTestPointsI = 1) {
  TestSuite t;
  using namespace Ikarus;
  using namespace Dune::Indices;
  using Manifold     = Ikarus::RealTuple<double, worldDim>;
  using Manifold2    = Ikarus::UnitVector<double, worldDim>;
  using VectorType   = Eigen::Vector<double, worldDim>;
  using MatrixType   = Eigen::Matrix<double, worldDim, worldDim>;
  constexpr int size = Manifold::valueSize;
  Dune::BlockVector<Manifold> testNodalPoints;
  using FECache = Dune::PQkLocalFiniteElementCache<double, double, domainDim, order>;
  FECache feCache;
  const auto &fe      = feCache.get(geometryType);
  auto localBasis     = Ikarus::LocalBasis(fe.localBasis());
  const size_t nNodes = fe.size();
  Dune::BlockVector<Manifold> testNodalPoints1;
  const int nNodalTestPoints = std::max(nNodalTestPointsI, nNodes);
  Ikarus::ValueFactory<Manifold>::construct(testNodalPoints1, nNodalTestPoints);

  Dune::BlockVector<Manifold2> testNodalPoints2;
  Ikarus::ValueFactory<Manifold2>::construct(testNodalPoints2, nNodalTestPoints);

  Dune::BlockVector<Manifold> vBlockedLocal(nNodes);
  Dune::BlockVector<Manifold> vBlockedLocal2(nNodes);
  Dune::BlockVector<Manifold2> vBlockedLocal3(nNodes);

  const auto &rule = Dune::QuadratureRules<double, domainDim>::rule(fe.type(), 2);
  localBasis.bind(rule, bindDerivatives(0, 1));

  for (size_t j = 0; j < fe.size(); ++j) {
    vBlockedLocal[j]  = testNodalPoints1[j];
    vBlockedLocal2[j] = testNodalPoints1[j];
    vBlockedLocal3[j] = testNodalPoints2[j];
  }

  auto f = Ikarus::StandardLocalFunction(localBasis, vBlockedLocal);
  t.subTest(testLocalFunction(f));
  static_assert(f.order() == linear);
  static_assert(countNonArithmeticLeafNodes(f) == 1);

  auto g = Ikarus::StandardLocalFunction(localBasis, vBlockedLocal);
  static_assert(countNonArithmeticLeafNodes(g) == 1);
  static_assert(g.order() == linear);

  auto h = f + g;
  t.subTest(testLocalFunction(h));
  static_assert(h.order() == linear);
  for (size_t k = 0; k < fe.size(); ++k) {
    t.check(h.coefficientsRef(_0)[k] == vBlockedLocal[k]);
    t.check(h.coefficientsRef(_1)[k] == vBlockedLocal[k]);
  }
  static_assert(std::tuple_size_v<decltype(collectNonArithmeticLeafNodes(h))> == 2);
  static_assert(countNonArithmeticLeafNodes(h) == 2);
  static_assert(
      std::is_same_v<typename decltype(h)::Ids, std::tuple<Dune::index_constant<0>, Dune::index_constant<0>>>);

  auto ft2 = 2 * f;
  t.subTest(testLocalFunction(ft2));
  static_assert(ft2.order() == linear);

  auto f23 = 2 * f * 3;
  t.subTest(testLocalFunction(f23));
  static_assert(f23.order() == linear);

  auto mf = -f;
  t.subTest(testLocalFunction(mf));
  static_assert(f.order() == mf.order());

  if constexpr (domainDim == worldDim) {
    auto eps = linearStrains(f);
    t.subTest(testLocalFunction(eps));
    static_assert(eps.order() == linear);
  }

  auto k = -dot(f + f, 3.0 * (g / 5.0) * 5.0);
  t.subTest(testLocalFunction(k));
  static_assert(k.order() == quadratic);
  static_assert(std::tuple_size_v<decltype(collectNonArithmeticLeafNodes(k))> == 3);
  static_assert(countNonArithmeticLeafNodes(k) == 3);

  auto dotfg = dot(f, g);
  t.subTest(testLocalFunction(dotfg));
  static_assert(countNonArithmeticLeafNodes(dotfg) == 2);
  static_assert(dotfg.order() == quadratic);
  static_assert(
      std::is_same_v<typename decltype(dotfg)::Ids, std::tuple<Dune::index_constant<0>, Dune::index_constant<0>>>);

  auto normSq = normSquared(f);
  t.subTest(testLocalFunction(normSq));
  static_assert(normSq.order() == quadratic);

  auto logg = log(dotfg);
  t.subTest(testLocalFunction(logg));

  auto powf = pow<3>(dotfg);
  t.subTest(testLocalFunction(powf));

  auto powfgsqrtdotfg = sqrt(powf);
  t.subTest(testLocalFunction(powfgsqrtdotfg));

  if constexpr (size > 1)  // Projection-Based only makes sense in 2d+
  {
    auto gP = Ikarus::ProjectionBasedLocalFunction(localBasis, vBlockedLocal3);
    static_assert(gP.order() == nonLinear);
    t.subTest(testLocalFunction(gP));
  }

  //  {
  //    auto localBasisNotBound = Ikarus::LocalBasis(fe.localBasis());
  //    auto fNotBound          = Ikarus::StandardLocalFunction(localBasisNotBound, vBlockedLocal);
  //    auto h1                 = f + fNotBound;
  //    //    (h1.viewOverIntegrationPoints(), "The basis of the leaf nodes are not in the same state.");
  //
  //    const auto &ruleHigher                 = Dune::QuadratureRules<double, domainDim>::rule(fe.type(), 7);
  //    auto localBasisBoundButToDifferentRule = Ikarus::LocalBasis(fe.localBasis());
  //    localBasisBoundButToDifferentRule.bind(ruleHigher, bindDerivatives(0, 1));
  //    auto fBoundButHigher = Ikarus::StandardLocalFunction(localBasisBoundButToDifferentRule, vBlockedLocal);
  ////    auto h2              = f + fBoundButHigher;
  //    //    (h2.viewOverIntegrationPoints(), "The basis of the leaf nodes are not in the same state.");
  //  }

  using namespace Ikarus::DerivativeDirections;

  const double tol = 1e-13;

  auto f2 = Ikarus::StandardLocalFunction(localBasis, vBlockedLocal, _0);
  auto g2 = Ikarus::StandardLocalFunction(localBasis, vBlockedLocal2, _1);
  static_assert(countNonArithmeticLeafNodes(f2) == 1);
  static_assert(countNonArithmeticLeafNodes(g2) == 1);

  auto k2 = dot(f2 + g2, g2);
  static_assert(countNonArithmeticLeafNodes(k2) == 3);
  static_assert(std::is_same_v<typename decltype(k2)::Ids,
                               std::tuple<Dune::index_constant<0>, Dune::index_constant<1>, Dune::index_constant<1>>>);

  auto b2 = collectNonArithmeticLeafNodes(k2);
  static_assert(std::tuple_size_v<decltype(b2)> == 3);

  for (int gpIndex = 0; [[maybe_unused]] auto &gp : rule) {
    const auto &N  = localBasis.evaluateFunction(gpIndex);
    const auto &dN = localBasis.evaluateJacobian(gpIndex);
    t.check(Dune::FloatCmp::eq(
        (f2.evaluateFunction(gpIndex) + g2.evaluateFunction(gpIndex)).dot(g2.evaluateFunction(gpIndex)),
        k2.evaluateFunction(gpIndex)[0]));
    auto resSingleSpatial
        = ((f2.evaluateDerivative(gpIndex, wrt(spatial(0))) + g2.evaluateDerivative(gpIndex, wrt(spatial(0))))
                   .transpose()
               * g2.evaluateFunction(gpIndex)
           + (f2.evaluateFunction(gpIndex) + g2.evaluateFunction(gpIndex)).transpose()
                 * g2.evaluateDerivative(gpIndex, wrt(spatial(0))))
              .eval();
    t.check(isApproxSame(resSingleSpatial, k2.evaluateDerivative(gpIndex, wrt(spatial(0))), tol));
    auto resSpatialAll
        = (((f2.evaluateDerivative(gpIndex, wrt(spatialAll)) + g2.evaluateDerivative(gpIndex, wrt(spatialAll)))
                .transpose()
            * g2.evaluateFunction(gpIndex))
               .transpose()
           + (f2.evaluateFunction(gpIndex) + g2.evaluateFunction(gpIndex)).transpose()
                 * g2.evaluateDerivative(gpIndex, wrt(spatialAll)))
              .eval();
    static_assert(resSpatialAll.cols() == domainDim);
    static_assert(resSpatialAll.rows() == 1);

    t.check(isApproxSame(resSpatialAll, k2.evaluateDerivative(gpIndex, wrt(spatialAll)), tol));

    for (size_t iC = 0; iC < fe.size(); ++iC) {
      const VectorType dfdi = g2.evaluateFunction(gpIndex) * N[iC];

      const VectorType dkdi = k2.evaluateDerivative(gpIndex, wrt(coeff(_0, iC)));

      t.check(isApproxSame(dfdi, dkdi, tol));

      for (size_t jC = 0; jC < fe.size(); ++jC) {
        const MatrixType dkdij         = k2.evaluateDerivative(gpIndex, wrt(coeff(_0, iC, _1, jC)));
        const MatrixType dkdijExpected = N[jC] * N[iC] * MatrixType::Identity();
        t.check(isApproxSame(dkdijExpected, dkdij, tol));

        const MatrixType dkdij2         = k2.evaluateDerivative(gpIndex, wrt(coeff(_0, iC, _0, jC)));
        const MatrixType dkdijExpected2 = MatrixType::Zero();
        t.check(isApproxSame(dkdijExpected2, dkdij2, tol));
        const MatrixType dkdij3         = k2.evaluateDerivative(gpIndex, wrt(coeff(_1, iC, _1, jC)));
        const MatrixType dkdijExpected3 = 2 * N[iC] * N[jC] * MatrixType::Identity();
        t.check(isApproxSame(dkdijExpected3, dkdij3, tol));

        const MatrixType dkdSij         = k2.evaluateDerivative(gpIndex, wrt(spatial(0), coeff(_0, iC, _1, jC)));
        const MatrixType dkdSijR        = k2.evaluateDerivative(gpIndex, wrt(coeff(_0, iC, _1, jC), spatial(0)));
        const MatrixType dkdSijExpected = (dN(jC, 0) * N[iC] + N[jC] * dN(iC, 0)) * MatrixType::Identity();
        t.check(isApproxSame(dkdSijR, dkdSij, tol));
        t.check(isApproxSame(dkdSijExpected, dkdSij, tol));
      }
    }
    ++gpIndex;
  }
  //  }
  return t;
}

// Most of the following tests are commented out due to very long compile times and long runtimes in debug mode we hope
//  to still capture most the bugs
using namespace Dune::GeometryTypes;
auto testExpressionsOnLine() {
  TestSuite t("testExpressionsOnLine");
  //    std::cout << "line with linear ansatz functions and 1d local function" << std::endl;
  //    localFunctionTestConstructor<1, 1, 1>(line);
  //  localFunctionTestConstructor<1, 2, 1>(line);  // line with linear ansatz functions and 2d lf
  //  std::cout << "line with linear ansatz functions and 3d local function" << std::endl;
  //  localFunctionTestConstructor<1, 3, 1>(line);
  //  std::cout << "line with quadratic ansatz functions and 1d local function" << std::endl;
  //  localFunctionTestConstructor<1, 1, 2>(line);
  //  localFunctionTestConstructor<1, 2, 2>(line);  // line with quadratic ansatz functions and 2d lf
  std::cout << "line with quadratic ansatz functions and 3d local function" << std::endl;
  t.subTest(localFunctionTestConstructor<1, 3, 2>(line));
  return t;
}

auto testExpressionsOnTriangle() {
  TestSuite t("testExpressionsOnTriangle");

  //  std::cout << "triangle with linear ansatz functions and 1d local function" << std::endl;
  //  localFunctionTestConstructor<2, 1, 1>(triangle);
  std::cout << "triangle with linear ansatz functions and 2d local function" << std::endl;
  t.subTest(localFunctionTestConstructor<2, 2, 1>(triangle));
  //  std::cout << "triangle with linear ansatz functions and 3d local function" << std::endl;
  //  localFunctionTestConstructor<2, 3, 1>(triangle);
  //  std::cout << "triangle with quadratic ansatz functions and 1d local function" << std::endl;
  //  localFunctionTestConstructor<2, 1, 2>(triangle);
  //  localFunctionTestConstructor<2, 2, 2>(triangle);  // triangle with quadratic ansatz functions and 2d lf
  std::cout << "triangle with quadratic ansatz functions and 3d local function" << std::endl;
  t.subTest(localFunctionTestConstructor<2, 3, 2>(triangle));
  return t;
}

auto testExpressionsOnQuadrilateral() {
  TestSuite t("testExpressionsOnQuadrilateral");
  //  std::cout << "quadrilateral with linear ansatz functions and 1d local function" << std::endl;
  //  localFunctionTestConstructor<2, 1, 1>(quadrilateral);
  std::cout << "quadrilateral with linear ansatz functions and 2d local function" << std::endl;
  t.subTest(localFunctionTestConstructor<2, 2, 1>(quadrilateral));
  //  std::cout << "quadrilateral with linear ansatz functions and 3d local function" << std::endl;
  //  localFunctionTestConstructor<2, 3, 1>(quadrilateral);
  //  std::cout << "quadrilateral with quadratic ansatz functions and 1d local function" << std::endl;
  //  localFunctionTestConstructor<2, 1, 2>(quadrilateral);
  //  localFunctionTestConstructor<2, 2, 2>(quadrilateral);  // quadrilateral with quadratic ansatz functions and 2d lf
  std::cout << "quadrilateral with quadratic ansatz functions and 3d local function" << std::endl;
  t.subTest(localFunctionTestConstructor<2, 3, 2>(quadrilateral));
  return t;
}

auto testExpressionsOnHexahedron() {
  TestSuite t("testExpressionsOnHexahedron");
  //  std::cout << "hexahedron with linear ansatz functions and 1d local function" << std::endl;
  //  localFunctionTestConstructor<3, 1, 1>(hexahedron);  // hexahedron with linear ansatz functions and 1d lf
  //  localFunctionTestConstructor<3, 2, 1>(hexahedron);  // hexahedron with linear ansatz functions and 2d lf
  std::cout << "hexahedron with linear ansatz functions and 3d local function" << std::endl;
  t.subTest(localFunctionTestConstructor<3, 3, 1>(hexahedron));
  //  std::cout << "hexahedron with quadratic ansatz functions and 1d local function" << std::endl;
  //  localFunctionTestConstructor<3, 1, 2>(hexahedron);
  //  localFunctionTestConstructor<3, 2, 2>(hexahedron);  // hexahedron with quadratic ansatz functions and 2d lf
  //  std::cout << "hexahedron with quadratic ansatz functions and 3d local function" << std::endl;
  //  localFunctionTestConstructor<3, 3, 2>(hexahedron);  // hexahedron with quadratic ansatz functions and 3d lf
  return t;
}

int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t;

  t.subTest(testExpressionsOnLine());
  t.subTest(testExpressionsOnTriangle());
  t.subTest(testExpressionsOnQuadrilateral());
  t.subTest(testExpressionsOnHexahedron());

  return t.exit();
}