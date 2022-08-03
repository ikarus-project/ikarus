
#include <config.h>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "common.hh"
#include "factories.hh"
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
#include <ikarus/utils/multiIndex.hh>

using namespace Dune::Functions::BasisFactory;

template <typename LB, bool isCopy = false>
void testLocalBasis(const LB& localBasis, const Dune::GeometryType& type) {
  const double tol = 1e-12;
  using namespace autodiff;
  using namespace Ikarus;

  constexpr int gridDim = LB::gridDim;
  const auto& rule      = Dune::QuadratureRules<double, gridDim>::rule(type, 3);

  for (const auto& gp : rule) {
    /// Check spatial derivatives
    {
      /// Check if spatial derivatives are really derivatives
      /// Perturb in a random direction in the elements parameter space and check spatial derivative
      auto func = [&](auto& gpOffset_) {
        Eigen::VectorXd N;
        localBasis.evaluateFunction(gp.position() + toFieldVector(gpOffset_), N);

        return N;
      };
      auto jacobianLambda = [&](auto& gpOffset_) {
        Eigen::Matrix<double, Eigen::Dynamic, gridDim> dN;
        localBasis.evaluateJacobian(gp.position() + toFieldVector(gpOffset_), dN);
        return dN.eval();
      };

      Eigen::Vector<double, gridDim> ipOffset = (Eigen::Vector<double, gridDim>::Random()).normalized() / 8;

      auto nonLinOpSpatialAll
          = Ikarus::NonLinearOperator(linearAlgebraFunctions(func, jacobianLambda), parameter(ipOffset));
      EXPECT_TRUE((checkJacobian<decltype(nonLinOpSpatialAll), Eigen::Vector<double, gridDim>>(
          nonLinOpSpatialAll, {.draw = false, .writeSlopeStatementIfFailed = true, .tolerance = 1e-2})));
      if constexpr (gridDim > 1) {
        std::cout << "Test Second Derivatives" << std::endl;

        constexpr auto iterC = voigtNotationContainer<gridDim>;
        auto iter            = iterC.begin();
        for (int i = 0; const auto [firstDirection, secondDirection] : voigtNotationContainer<gridDim>) {
          std::cout << "Test Mixed Directions: " << firstDirection << " " << secondDirection << std::endl;
          auto jacobianLambda1D = [&](const auto& gpOffset_) {
            Eigen::Matrix<double, Eigen::Dynamic, gridDim> dN;
            Dune::FieldVector<double, gridDim> gpOffset2D;
            std::ranges::fill(gpOffset2D, 0);
            gpOffset2D[firstDirection] = gpOffset_[0];
            localBasis.evaluateJacobian(gp.position() + gpOffset2D, dN);
            return dN.col(secondDirection).eval();
          };
          constexpr int secondDerivatives = gridDim * (gridDim + 1) / 2;
          auto hessianLambda              = [&](const auto& gpOffset_) {
            Eigen::Matrix<double, Eigen::Dynamic, secondDerivatives> ddN;
            Dune::FieldVector<double, gridDim> gpOffset2D;
            std::ranges::fill(gpOffset2D, 0);
            gpOffset2D[firstDirection] = gpOffset_[0];
            localBasis.evaluateSecondDerivatives(gp.position() + gpOffset2D, ddN);
            return ddN.col(i).eval();
          };

          Eigen::Vector<double, 1> ipOffset1D(1);

          auto nonLinOpHg = Ikarus::NonLinearOperator(linearAlgebraFunctions(jacobianLambda1D, hessianLambda),
                                                      parameter(ipOffset1D));

          EXPECT_TRUE((checkJacobian<decltype(nonLinOpHg), Eigen::Vector<double, 1>>(
              nonLinOpHg, {.draw = false, .writeSlopeStatementIfFailed = true, .tolerance = 1e-2})));
          ++i;
        }
      }
    }
  }
}

template <int domainDim, int order>
void localBasisTestConstructor(const Dune::GeometryType& geometryType, size_t nNodalTestPointsI = 6) {
  using namespace Ikarus;
  using namespace Dune::Indices;

  using FECache = Dune::PQkLocalFiniteElementCache<double, double, domainDim, order>;
  FECache feCache;
  const auto& fe  = feCache.get(geometryType);
  auto localBasis = Ikarus::LocalBasis(fe.localBasis());

  const auto& rule = Dune::QuadratureRules<double, domainDim>::rule(fe.type(), 3);
  testLocalBasis(localBasis, geometryType);
}

TEST(LocalBasisTests, TestLocalBasis) {
  using namespace Dune::GeometryTypes;
  std::cout << "Test line with linear ansatz functions" << std::endl;
  localBasisTestConstructor<1, 1>(line);
  std::cout << "Test line with quadratic ansatz functions" << std::endl;
  localBasisTestConstructor<1, 2>(line);
  std::cout << "Test triangle with linear ansatz functions" << std::endl;
  localBasisTestConstructor<2, 1>(triangle);
  std::cout << "Test triangle with quadratic ansatz functions" << std::endl;
  localBasisTestConstructor<2, 2>(triangle);
  std::cout << "Test quadrilateral with linear ansatz functions" << std::endl;
  localBasisTestConstructor<2, 1>(quadrilateral);
  std::cout << "Test quadrilateral with quadratic ansatz functions" << std::endl;
  localBasisTestConstructor<2, 2>(quadrilateral);
  std::cout << "Test hexahedron with linear ansatz functions" << std::endl;
  localBasisTestConstructor<3, 1>(hexahedron);
  std::cout << "Test hexahedron with quadratic ansatz functions" << std::endl;
  localBasisTestConstructor<3, 2>(hexahedron);
}
