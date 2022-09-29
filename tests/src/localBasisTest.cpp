//
#include <config.h>
#include <dune/common/test/testsuite.hh>
#include <dune/common/parallel/mpihelper.hh>
using Dune::TestSuite;

#include "common.hh"
#include <complex>

#include <dune/common/classname.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <Eigen/Core>


#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/localFunctions/expressions.hh>
#include <ikarus/manifolds/unitVector.hh>
#include <ikarus/utils/functionSanityChecks.hh>
#include <ikarus/utils/multiIndex.hh>

using namespace Dune::Functions::BasisFactory;

template <typename LB, bool isCopy = false>
auto testLocalBasis(LB& localBasis, const Dune::GeometryType& type) {
  TestSuite t("testLocalBasis");
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
      t.check((checkJacobian<decltype(nonLinOpSpatialAll), Eigen::Vector<double, gridDim>>(
          nonLinOpSpatialAll, {.draw = false, .writeSlopeStatementIfFailed = true, .tolerance = 1e-2})));
      if constexpr (gridDim > 1) {
        std::cout << "Test Second Derivatives" << std::endl;

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

            t.check((checkJacobian<decltype(nonLinOpHg), Eigen::Vector<double, 1>>(
              nonLinOpHg, {.draw = false, .writeSlopeStatementIfFailed = true, .tolerance = 1e-2})));
          ++i;
        }
      }
    }
  }
  // Unbound basis checks
  t.check(not localBasis.isBound());
    try {
        localBasis.evaluateFunction(0);
        t.check(false,"The prior function call should have thrown! You should not end up here.");
    }catch(const std::logic_error& )
    { }
    try {
        localBasis.evaluateJacobian(0);
        t.check(false,"The prior function call should have thrown! You should not end up here.");
    }catch(const std::logic_error& )
    { }
  if constexpr (gridDim > 1) {
      try {
          localBasis.evaluateSecondDerivatives(0);
          t.check(false,"The prior function call should have thrown! You should not end up here.");
      }catch(const std::logic_error& )
      { }
    localBasis.bind(rule, bindDerivatives(0, 1, 2));
  } else {
      localBasis.bind(rule, bindDerivatives(0, 1)); }
    t.check(localBasis.isBound());
    return t;
}

template <int domainDim, int order>
auto localBasisTestConstructor(const Dune::GeometryType& geometryType,[[maybe_unused]] size_t nNodalTestPointsI = 6) {
  using namespace Ikarus;
  using namespace Dune::Indices;

  using FECache = Dune::PQkLocalFiniteElementCache<double, double, domainDim, order>;
  FECache feCache;
  const auto& fe  = feCache.get(geometryType);
  auto localBasis = Ikarus::LocalBasis(fe.localBasis());

  return testLocalBasis(localBasis, geometryType);
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    TestSuite t;
    using namespace Dune::GeometryTypes;
    std::cout << "Test line with linear ansatz functions" << std::endl;
    t.subTest(localBasisTestConstructor<1, 1>(line));
    std::cout << "Test line with quadratic ansatz functions" << std::endl;
    t.subTest(localBasisTestConstructor<1, 2>(line));
    std::cout << "Test triangle with linear ansatz functions" << std::endl;
    t.subTest(localBasisTestConstructor<2, 1>(triangle));
    std::cout << "Test triangle with quadratic ansatz functions" << std::endl;
    t.subTest(localBasisTestConstructor<2, 2>(triangle));
    std::cout << "Test quadrilateral with linear ansatz functions" << std::endl;
    t.subTest(localBasisTestConstructor<2, 1>(quadrilateral));
    std::cout << "Test quadrilateral with quadratic ansatz functions" << std::endl;
    t.subTest(localBasisTestConstructor<2, 2>(quadrilateral));
    std::cout << "Test hexahedron with linear ansatz functions" << std::endl;
    t.subTest(localBasisTestConstructor<3, 1>(hexahedron));
    std::cout << "Test hexahedron with quadratic ansatz functions" << std::endl;
    t.subTest(localBasisTestConstructor<3, 2>(hexahedron));

    return t.exit();
}