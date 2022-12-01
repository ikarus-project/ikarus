#include <config.h>

#include "testHelpers.hh"

#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/fufem/dunepython.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
using Dune::TestSuite;
#include <dune/common/float_cmp.hh>
#include <dune/functions/functionspacebases/interpolate.hh>

#include <Eigen/Core>

#include <ikarus/utils/duneUtilities.hh>

auto dirichletBCTest() {
  TestSuite t("SimpleAssemblersTest");
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox       = {4, 2};
  std::array<int, 2> elementsPerDirection = {2, 1};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);

  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeConstSharedBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));

  Ikarus::DirichletValues dirichletValues1(basis);
  dirichletValues1.fixDOFs([](auto& basis_, auto& dirichFlags) {
    Dune::Functions::forEachBoundaryDOF(basis_, [&](auto&& indexGlobal) { dirichFlags[indexGlobal] = true; });
  });

  Ikarus::DirichletValues dirichletValues2(basis);
  dirichletValues2.fixBoundaryDOFs([](auto& dirichFlags, auto&& indexGlobal) { dirichFlags[indexGlobal] = true; });

  for (std::size_t i = 0; i < basis->size(); ++i)
    t.check(dirichletValues1.isConstrained(i) == dirichletValues2.isConstrained(i))
        << "Different dirichlet value creations didn't provide the same result. Index: i=" << i;

  auto inhomogeneousDisplacement = []<typename T>(const auto& globalCoord, const T& lambda) {
    Eigen::Vector<T, 2> localInhomogeneous;
    if (globalCoord[0] > 4 - 1e-8) {
      localInhomogeniuous[0] = 4 * lambda;
      localInhomogeniuous[1] = 4 * lambda;
    } else
      localInhomogeniuous.setZero();
    return localInhomogeniuous;
  };

  auto inhomogeneousDisplacement2 = []<typename T>(const auto& globalCoord, const T& lambda) {
    Eigen::Vector<T, 2> localInhomogeneous;
    if (globalCoord[0] < 1e-8) {
      localInhomogeniuous[0] = 7 * lambda;
      localInhomogeniuous[1] = 7 * lambda;
    } else
      localInhomogeniuous.setZero();
    return localInhomogeniuous;
  };

  dirichletValues1.storeInhomogeneousBoundaryCondition(inhomogeneousDisplacement);
  dirichletValues1.storeInhomogeneousBoundaryCondition(inhomogeneousDisplacement2);
  Eigen::VectorXd disps, dispDerivs;
  dirichletValues1.evaluateInhomogeneousBoundaryCondition(disps, 2);
  dirichletValues1.evaluateInhomogeneousBoundaryConditionDerivative(dispDerivs, 2);

  auto lambdaCheck = [&](auto&& localIndex, auto&& localView, auto&& intersection) {
    auto globalIndex = localView.index(localIndex);
    if (intersection.geometry().center()[0] > 4 - 1e-8) {
      t.check(Dune::FloatCmp::eq(disps[globalIndex[0]], 8.0)) << "Values differ disps[i]: " << disps[globalIndex[0]];
      t.check(Dune::FloatCmp::eq(dispDerivs[globalIndex[0]], 4.0))
          << "Values differ dispDerivs[i]: " << dispDerivs[globalIndex[0]];
    } else if (intersection.geometry().center()[0] < 1e-8) {
      t.check(Dune::FloatCmp::eq(disps[globalIndex[0]], 14.0)) << "Values differ disps[i]: " << disps[globalIndex[0]];
      t.check(Dune::FloatCmp::eq(dispDerivs[globalIndex[0]], 7.0))
          << "Values differ dispDerivs[i]: " << dispDerivs[globalIndex[0]];
    }
  };
  Dune::Functions::forEachBoundaryDOF(*basis, lambdaCheck);

  // Check that we can store lambda from python
  std::string inhomogeneousDisplacementFunction
      = std::string("lambda globalCoord,lam: ( numpy.array([1*lam*globalCoord[0], 2*lam*globalCoord[1]]) )");

  Python::start();
  Python::Reference main = Python::import("__main__");
  Python::run("import math");
  Python::run("import numpy");
  auto pythonFuncdouble
      = Python::make_function<Dune::FieldVector<double, 2>>(Python::evaluate(inhomogeneousDisplacementFunction));
  auto pythonFuncdual = Python::make_function<Dune::FieldVector<autodiff::real, 2>>(
      Python::evaluate(inhomogeneousDisplacementFunction));
  double lambda             = 7.5;
  auto resVal               = Ikarus::toEigen(pythonFuncdouble(Dune::FieldVector<double, 2>({1.0, 2.0}), lambda));
  autodiff::real lambdadual = lambda;
  lambdadual[1] = 1;

  auto res                      = pythonFuncdual(Dune::FieldVector<double, 2>({1.0, 2.0}), lambdadual);
  auto valueResult              = grad<0>(Ikarus::toEigen(res));
  auto derivativeResult         = grad(Ikarus::toEigen(res));
  auto derivativeResultExpected = Eigen::Vector<double, 2>({1.0, 4.0});
  t.check(valueResult.isApprox(resVal)) << "valueResult's value is not " << derivativeResultExpected << " but "
                                        << derivativeResult;
  t.check(derivativeResult.isApprox(derivativeResultExpected))
      << "derivativeResult's value is not " << derivativeResultExpected << " but " << derivativeResult;

  auto inhomogeneousDisplacementFromPython = [&]<typename T>(const auto& globalCoord, const T& lambda) {
    auto pythonFunc
        = Python::make_function<Dune::FieldVector<T, 2>>(Python::evaluate(inhomogeneousDisplacementFunction));

    return Ikarus::toEigen(pythonFunc(globalCoord, lambda));
  };
  //
  Ikarus::DirichletValues dirichletValues3(basis);
  dirichletValues3.fixBoundaryDOFs([](auto& dirichFlags, auto&& indexGlobal) { dirichFlags[indexGlobal] = true; });

  Eigen::VectorXd disps2, dispDerivs2;
  dirichletValues3.storeInhomogeneousBoundaryCondition(inhomogeneousDisplacementFromPython);
  double lambda3 = 6;
  dirichletValues3.evaluateInhomogeneousBoundaryCondition(disps2, lambda3);
  dirichletValues3.evaluateInhomogeneousBoundaryConditionDerivative(dispDerivs2, lambda3);
  t.check(disps2.isApprox(dispDerivs2 * lambda3));

  return t;
}

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t;

  t.subTest(dirichletBCTest());
  return t.exit();
}