//
// Created by Alex on 21.07.2021.
//
#include <../../config.h>
#include <numbers>

#include <dune/common/float_cmp.hh>
#include <dune/common/indices.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/onedgrid.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/LocalBasis/localBasis.h"
#include "ikarus/basis/basishelper.h"
#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/Solver/LinearSolver/LinearSolver.h>
using namespace Dune::Functions::BasisFactory;
using namespace Dune::Functions;
using namespace Dune::Functions::BasisBuilder;
using namespace Dune::Indices;

void exampleTrussElement() {
  const double L        = 1;
  const double EA       = 1.0;
  const int numElements = 10;
  const int numGP       = 2;
  Dune::OneDGrid grid(numElements, 0, L);
  auto gridView = grid.leafGridView();
  // draw(gridView);
  using namespace Dune::Functions::BasisFactory;
  auto basis     = makeBasis(gridView, lagrange<1>());
  auto localView = basis.localView();

  auto numDofs           = basis.size();
  Eigen::MatrixXd K_Glob = Eigen::MatrixXd::Zero(numDofs, numDofs);

  for (auto& ele : elements(gridView)) {
    localView.bind(ele);
    auto& fe = localView.tree().finiteElement();
    Ikarus::LocalBasis localBasis(fe.localBasis());
    const auto& rule = Dune::QuadratureRules<double, 1>::rule(ele.type(), numGP, Dune::QuadratureType::GaussLegendre);

    Eigen::VectorXd dNdxi = Eigen::VectorXd::Zero(2);
    Eigen::MatrixXd K     = Eigen::MatrixXd::Zero(2, 2);
    Eigen::VectorXd B     = Eigen::VectorXd::Zero(2);

    auto detJ = ele.geometry().integrationElement(0.0);
    for (auto& gp : rule) {
      localBasis.evaluateJacobian(gp.position(), dNdxi);
      B = dNdxi / detJ;
      K += EA * B * B.transpose() * detJ * gp.weight();
    }

    // Adding local stiffness the global stiffness
    for (auto i = 0U; i < localView.size(); ++i)
      for (auto j = 0U; j < localView.size(); ++j)
        K_Glob(localView.index(i)[0], localView.index(j)[0]) += K(i, j);
  }

  // external force (1 one the right end)
  Eigen::VectorXd F_ExtGlobRed = Eigen::VectorXd::Zero(numDofs - 1);
  F_ExtGlobRed(numDofs - 2)    = 1.0;

  // reduce stiffness matrix (fix left end)
  const Eigen::MatrixXd K_GlobRed = K_Glob(Eigen::seq(1, Eigen::last), Eigen::seq(1, Eigen::last));

  // solve the linear system
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT);
  linSolver.factorize(K_GlobRed);
  const Eigen::VectorXd D_GlobRed = linSolver.solve(F_ExtGlobRed);

  std::cout << "Displacement: " << D_GlobRed(Eigen::last) << "\n";
}

Eigen::MatrixXd TimoshenkoBeamStiffness(auto localView, auto gridElement, auto quadratureRule,
                                        const Eigen::Matrix2d& C) {
  using namespace Dune::Indices;

  Ikarus::LocalBasis basisW(localView.tree().child(_0).finiteElement().localBasis());
  Ikarus::LocalBasis basisPhi(localView.tree().child(_1).finiteElement().localBasis());

  // Determinant of Jacobian, obtained from gridElement
  auto detJ = gridElement.geometry().volume();

  // get number of DOFs for w and phi
  auto numDofsW      = basisW.size();
  auto numDofsPhi    = basisPhi.size();
  auto numDofsPerEle = numDofsW + numDofsPhi;

  // initialize quantities
  Eigen::MatrixXd K        = Eigen::MatrixXd::Zero(numDofsPerEle, numDofsPerEle);
  Eigen::VectorXd dNwDxi   = Eigen::VectorXd::Zero(numDofsW);
  Eigen::VectorXd NphiDxi  = Eigen::VectorXd::Zero(numDofsPhi);
  Eigen::VectorXd dNphiDxi = Eigen::VectorXd::Zero(numDofsPhi);

  // integration point loop
  for (auto& gp : quadratureRule) {
    // evaluate ansatz functions and their derivatives
    basisW.evaluateJacobian(gp.position(), dNwDxi);
    basisPhi.evaluateFunction(gp.position(), NphiDxi);
    basisPhi.evaluateJacobian(gp.position(), dNphiDxi);

    // setup B-operator
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2, numDofsPerEle);

    // fill columns of B-Operator related to w-DOFs
    for (unsigned int i = 0; i < localView.tree().child(_0).size(); ++i) {
      B(1, localView.tree().child(_0).localIndex(i)) = dNwDxi[i] / detJ;
    }

    // fill columns of B-Operator related to phi-DOFs
    for (unsigned int i = 0; i < localView.tree().child(_1).size(); ++i) {
      B(0, localView.tree().child(_1).localIndex(i)) = dNphiDxi[i] / detJ;
      B(1, localView.tree().child(_1).localIndex(i)) = NphiDxi[i];
    }

    // integration of stiffness matrix
    K += B.transpose() * C * B * detJ * gp.weight();
  }

  return K;
}

enum class TimoschenkoBeam { w, phi };

unsigned int getGlobalDofIdImpl(const auto& basis, const double position) {
  auto localView       = basis.localView();
  auto seDOFs          = subEntityDOFs(basis);
  const auto& gridView = basis.gridView();
  for (auto& element : elements(gridView)) {
    localView.bind(element);
    for (unsigned int i = 0; i < element.subEntities(1); ++i) {
      if (Dune::FloatCmp::eq(element.template subEntity<1>(i).geometry().center()[0], position, 1e-8)) {
        auto& localIndex = seDOFs.bind(localView, i, 1);
        assert(localIndex.size() == 1 && "It is expected that only one w-DOF is associated with a vertex");
        return localView.index(localIndex[0])[0];
      }
    }
  }
  throw std::runtime_error(
      "There is no displacement dof at the requested position. Currently, only DOFs at vertices are supported.");
}

unsigned int getGlobalDofId(TimoschenkoBeam requestedQuantity, const auto& basis, const double position) {
  using namespace Dune::Indices;
  if (requestedQuantity == TimoschenkoBeam::w)
    return getGlobalDofIdImpl(subspaceBasis(basis, _0), position);
  else if (requestedQuantity == TimoschenkoBeam::phi)
    return getGlobalDofIdImpl(subspaceBasis(basis, _1), position);
  else
    throw std::runtime_error("The requested quantity is not supported");
}

void plotDeformedTimoschenkoBeam(auto& gridView, auto& basis, auto& d_glob, double EI, double GA, double L, double F) {
  using namespace Dune::Indices;
  auto wGlobal   = makeDiscreteGlobalBasisFunction<double>(subspaceBasis(basis, _0), d_glob);
  auto phiGlobal = makeDiscreteGlobalBasisFunction<double>(subspaceBasis(basis, _1), d_glob);

  auto wSol = [&](auto x) {
    return -F * Dune::power(x[0], 3) / (6.0 * EI) + L * F * Dune::power(x[0], 2) / (2.0 * EI) + F * x[0] / GA;
  };
  auto phiSol            = [&](auto x) { return F * Dune::power(x[0], 2) / (2.0 * EI) - L * F * x[0] / EI; };
  auto wGlobalAnalytic   = makeAnalyticGridViewFunction(wSol, gridView);
  auto phiGlobalAnalytic = makeAnalyticGridViewFunction(phiSol, gridView);

  using namespace matplot;
  auto f = figure(true);
  tiledlayout(1, 2);
  auto ax1 = nexttile();
  auto ax2 = nexttile();
  hold(ax1, true);
  hold(ax2, true);

  auto localView        = basis.localView();
  auto wLocal           = localFunction(wGlobal);
  auto phiLocal         = localFunction(phiGlobal);
  auto wLocalAnalytic   = localFunction(wGlobalAnalytic);
  auto phiLocalAnalytic = localFunction(phiGlobalAnalytic);
  std::vector<double> x = linspace(0, 1, 10);
  std::vector<double> x_L;
  std::vector<double> yw, yphi, ywAna, yphiAna;
  for (auto& edge : elements(gridView)) {
    wLocal.bind(edge);
    wLocalAnalytic.bind(edge);
    phiLocal.bind(edge);
    phiLocalAnalytic.bind(edge);
    localView.bind(edge);
    x_L     = transform(x, [&](auto x) { return edge.geometry().global({x}); });
    yw      = transform(x, [&](auto x) { return wLocal({x}); });
    ywAna   = transform(x, [&](auto x) { return wLocalAnalytic({x}); });
    yphi    = transform(x, [&](auto x) { return phiLocal({x}); });
    yphiAna = transform(x, [&](auto x) { return phiLocalAnalytic({x}); });

    auto l0 = ax1->plot(x_L, yw);
    l0->line_width(2);
    l0->color("blue");

    auto l0_ana = ax1->plot(x_L, ywAna);
    l0_ana->line_width(2);
    l0_ana->color("red");

    auto l1 = ax2->plot(x_L, yphi);
    l1->line_width(2);
    l1->color("blue");

    auto l1_ana = ax2->plot(x_L, yphiAna);
    l1_ana->line_width(2);
    l1_ana->color("red");
  }

  f->show();
}

void exampleTimoshenkoBeam(const int polynomialOrderW, const int polynomialOrderPhi,const int numElements) {
  const double b  = 1;
  const double L  = 1e3;
  const double E  = 1e8;
  const double G  = E/2;
  const double t  = 1e-3;
  const double EI = E * b * t * t * t / 12.0;
  const double GA = G * b * t;
  const double F  = 0.1 * t * t * t;
  Eigen::Matrix2d C;
  C << EI, 0, 0, GA;
  const int maxOrderIntegration    = std::max(2*(polynomialOrderW-1),2*polynomialOrderPhi);
  Dune::OneDGrid grid(numElements, 0, L);
  auto gridView = grid.leafGridView();
  //draw(gridView);

  // Basis with different orders for w (first) and phi (second)
  auto basis     = makeBasis(gridView,composite(lagrange(polynomialOrderW), lagrange(polynomialOrderPhi), FlatLexicographic()));
  auto localView = basis.localView();

  // global stiffness matrix and force vector
  auto numDofs              = basis.size();
  Eigen::VectorXd F_ExtGlob = Eigen::VectorXd::Zero(numDofs);
  Eigen::MatrixXd K_Glob    = Eigen::MatrixXd::Zero(numDofs, numDofs);

  for (auto& ele : elements(gridView)) {
    localView.bind(ele);

    // Define the integration rule
    const auto& rule = Dune::QuadratureRules<double, 1>::
            rule(ele.type(), maxOrderIntegration, Dune::QuadratureType::GaussLegendre);

    // get local stiffness matrix
    auto K_local = TimoshenkoBeamStiffness(localView, ele, rule, C);

    // Adding local stiffness the global stiffness
    for (auto i = 0U; i < localView.size(); ++i)
      for (auto j = 0U; j < localView.size(); ++j)
        K_Glob(localView.index(i)[0], localView.index(j)[0]) += K_local(i, j);
  }

  // apply load on the right-hand side
  F_ExtGlob(getGlobalDofId(TimoschenkoBeam::w, basis, L)) = F;

  // clamp left-hand side
  std::vector<unsigned int> fixedDofs{getGlobalDofId(TimoschenkoBeam::w, basis, 0.0),
                                      getGlobalDofId(TimoschenkoBeam::phi, basis, 0.0)};
  for (auto dof : fixedDofs) {
    K_Glob.col(dof).setZero();
    K_Glob.row(dof).setZero();
    K_Glob(dof, dof) = 1.0;
  }

  // solve the linear system
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT);
  linSolver.factorize(K_Glob);
  const Eigen::VectorXd D_Glob = linSolver.solve(F_ExtGlob);
  // analytical solution
  //std::cout << "Bernoulli solution for displacement at L: " << F * L * L * L / (3.0 * EI) << "\n";

  // plot the result

  plotDeformedTimoschenkoBeam(gridView, basis, D_Glob, EI, GA, L, F);
}

int main() {
  //  exampleTrussElement();
  exampleTimoshenkoBeam(1,1,1);
  exampleTimoshenkoBeam(2,1,1);
  exampleTimoshenkoBeam(2,2,1);
  exampleTimoshenkoBeam(3,2,1);
}