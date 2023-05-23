// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#include <config.h>

#include "testHelpers.hh"

#include <vector>

#include <dune/common/float_cmp.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/fufem/dunepython.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/localfefunctions/eigenDuneTransformations.hh>

#include <Eigen/Core>

#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/duneUtilities.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

auto dirichletBCTest() {
  TestSuite t("SimpleAssemblersTest");
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox       = {4, 2};
  std::array<int, 2> elementsPerDirection = {2, 1};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);

  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>()));

  auto basisP = std::make_shared<const decltype(basis)>(basis);

  Ikarus::DirichletValues dirichletValues1(basisP->flat());
  dirichletValues1.fixDOFs([](auto& basis_, auto& dirichFlags) {
    Dune::Functions::forEachBoundaryDOF(basis_, [&](auto&& indexGlobal) { dirichFlags[indexGlobal] = true; });
  });

  Ikarus::DirichletValues dirichletValues2(basisP->flat());
  dirichletValues2.fixBoundaryDOFs([](auto& dirichFlags, auto&& indexGlobal) { dirichFlags[indexGlobal] = true; });

  for (std::size_t i = 0; i < basisP->flat().size(); ++i)
    t.check(dirichletValues1.isConstrained(i) == dirichletValues2.isConstrained(i))
        << "Different dirichlet value creations didn't provide the same result. Index: i=" << i;

  auto inhomogeneousDisplacement = []<typename T>(const auto& globalCoord, const T& lambda) {
    Eigen::Vector<T, 2> localInhomogeneous;
    if (globalCoord[0] > 4 - 1e-8) {
      localInhomogeneous[0] = 4 * lambda;
      localInhomogeneous[1] = 4 * lambda;
    } else
      localInhomogeneous.setZero();
    return localInhomogeneous;
  };

  auto inhomogeneousDisplacement2 = []<typename T>(const auto& globalCoord, const T& lambda) {
    Eigen::Vector<T, 2> localInhomogeneous;
    if (globalCoord[0] < 1e-8) {
      localInhomogeneous[0] = 7 * lambda;
      localInhomogeneous[1] = 7 * lambda;
    } else
      localInhomogeneous.setZero();
    return localInhomogeneous;
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
  Dune::Functions::forEachBoundaryDOF(basisP->flat(), lambdaCheck);

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
  auto resVal               = Dune::toEigen(pythonFuncdouble(Dune::FieldVector<double, 2>({1.0, 2.0}), lambda));
  autodiff::real lambdadual = lambda;
  lambdadual[1]             = 1;

  auto res                      = Dune::toEigen(pythonFuncdual(Dune::FieldVector<double, 2>({1.0, 2.0}), lambdadual));
  auto valueResult              = grad<0>(res);
  auto derivativeResult         = grad(res);
  auto derivativeResultExpected = Eigen::Vector<double, 2>({1.0, 4.0});
  t.check(valueResult.isApprox(resVal)) << "valueResult's value is not " << derivativeResultExpected << " but "
                                        << derivativeResult;
  t.check(derivativeResult.isApprox(derivativeResultExpected))
      << "derivativeResult's value is not " << derivativeResultExpected << " but " << derivativeResult;

  auto inhomogeneousDisplacementFromPython = [&]<typename T>(const auto& globalCoord, const T& lambda_) {
    auto pythonFunc
        = Python::make_function<Dune::FieldVector<T, 2>>(Python::evaluate(inhomogeneousDisplacementFunction));

    return Dune::toEigen(pythonFunc(globalCoord, lambda_));
  };
  //
  Ikarus::DirichletValues dirichletValues3(basisP->flat());
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
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(dirichletBCTest());
  return t.exit();
}
