//
// Created by Alex on 21.07.2021.
//
#include <../../config.h>
#include <numbers>

#include <dune/common/indices.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/onedgrid.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/LocalBasis/localBasis.h"
#include "ikarus/basis/basishelper.h"
#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/Solver/LinearSolver/LinearSolver.h>

void exampleTrussElement(){
  constexpr int griddim = 1;
  const double L = 1;
  const double EA = 1.0;
  const int numElements = 10;
  const int numGP = 2;
  constexpr int polynomialOrder = 2;
  Dune::OneDGrid grid(numElements, 0, L);
  //grid.globalRefine(2);
  auto gridView = grid.leafGridView();
  std::cout << gridView.dimensionworld << "\n";
  draw(gridView);
  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, lagrange<polynomialOrder>());
  auto basisView = basis.localView();


  for (auto& ele : elements(gridView)) {
    basisView.bind(ele);
    auto& fe = basisView.tree().finiteElement();
    Ikarus::LocalBasis localBasis(fe.localBasis());
    const auto& rule = Dune::QuadratureRules<double, 1>::rule(ele.type(),numGP,Dune::QuadratureType::GaussLegendre);

    Eigen::VectorXd dNdxi = Eigen::VectorXd::Zero(polynomialOrder+1);
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(polynomialOrder+1,polynomialOrder+1);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(polynomialOrder+1);

    auto detJ = ele.geometry().integrationElement(0.0);
    for (auto& gp : rule) {
      localBasis.evaluateJacobian(gp.position(),dNdxi);
      B = dNdxi/detJ;
      K += EA*B*B.transpose()*detJ*gp.weight();
    }
    std::cout << "The stiffness matrix is: " << "\n";
    std::cout << K << "\n";
  }
}

Eigen::MatrixXd TimoshenkoBeamStiffness (auto basisView, auto gridElement, auto quadratureRule, const Eigen::Matrix2d& C){
  using namespace Dune::Indices;

  Ikarus::LocalBasis basisW(basisView.tree().child(_0).finiteElement().localBasis());
  Ikarus::LocalBasis basisPhi(basisView.tree().child(_1).finiteElement().localBasis());

  // Determinant of Jacobian, obtained from gridElement
  auto detJ = gridElement.geometry().integrationElement(0.0);

  // get number of DOFs for w and phi
  auto numDofsW = basisW.size();
  auto numDofsPhi = basisPhi.size();
  auto numDofsPerEle = numDofsW + numDofsPhi;

  // initialize quantities
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(numDofsPerEle,numDofsPerEle);
  Eigen::VectorXd dNwDxi = Eigen::VectorXd::Zero(numDofsW);
  Eigen::VectorXd NphiDxi = Eigen::VectorXd::Zero(numDofsPhi);
  Eigen::VectorXd dNphiDxi = Eigen::VectorXd::Zero(numDofsPhi);

  // integration point loop
  for (auto& gp : quadratureRule) {
    // evaluate ansatz functions and their derivatives
    basisW.evaluateJacobian(gp.position(),dNwDxi);
    basisPhi.evaluateFunction(gp.position(),NphiDxi);
    basisPhi.evaluateJacobian(gp.position(),dNphiDxi);

    // setup B-operator
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2,numDofsPerEle);

    // fill columns of B-Operator related to w-DOFs
    for (unsigned int i = 0; i < basisView.tree().child(_0).size(); ++i) {
      B(1,basisView.tree().child(_0).localIndex(i)) = dNwDxi[i]/detJ;
    }

    // fill columns of B-Operator related to phi-DOFs
    for (unsigned int i = 0; i < basisView.tree().child(_1).size(); ++i) {
      B(0,basisView.tree().child(_1).localIndex(i)) = dNphiDxi[i]/detJ;
      B(1,basisView.tree().child(_1).localIndex(i)) = NphiDxi[i];
    }

    // integration of stiffness matrix
    K += B.transpose()*C*B*detJ*gp.weight();
  }

  return K;
}

enum class TimoschenkoBeam {w,phi};

unsigned int getGlobalDofIdImpl(const auto& basis, const double position){
  auto localView = basis.localView();
  auto seDOFs = subEntityDOFs(basis);
  const auto& gridView = basis.gridView();
  for(auto& element : elements(gridView)){
    localView.bind(element);
    for (unsigned int i = 0; i < element.subEntities(1); ++i) {
      if (element.template subEntity<1>(i).geometry().center() == position)
      {
        auto& localIndex = seDOFs.bind(localView,i,1);
        assert(localIndex.size() == 1 && "It is expected that only one w-DOF is associated with a vertex");
        return localView.index(localIndex[0])[0];
      }
    }
  }
  throw std::runtime_error("There is no displacement dof at the requested position. Currently, only DOFs at vertices are supported.");
}

unsigned int getGlobalDofId(TimoschenkoBeam requestedQuantity, const auto& basis, const double position){
  using namespace Dune::Indices;
  if (requestedQuantity == TimoschenkoBeam::w)
    return getGlobalDofIdImpl(subspaceBasis(basis,_0), position);
  else if (requestedQuantity == TimoschenkoBeam::phi)
    return getGlobalDofIdImpl(subspaceBasis(basis,_1), position);
  else
    throw std::runtime_error("The requested quantity is not supported");
}

void plotDeformedTimoschenkoBeam(auto& basis, auto& gridView, auto& dofVector){
  using namespace matplot;
  using namespace Dune::Indices;
  auto f  = figure(true);
  auto ax = gca();
  hold(ax, true);
  auto subBasisView = subspaceBasis(basis,_0).localView();
  auto seDOFs = subEntityDOFs(basis);
  for (auto& edge : elements(gridView)) {
    subBasisView.bind(edge);
    std::array<double, 2> xEdge{}, yEdge{};
    for (unsigned int i = 0; i < edge.subEntities(1); ++i) {
      auto localIndex = seDOFs.bind(subBasisView,i,1);
      assert(localIndex.size() == 1 && "It is expected that only one w-DOF is associated with a vertex");
      xEdge[i] = edge.template subEntity<1>(i).geometry().center();
      yEdge[i] = -dofVector[subBasisView.index(localIndex[0])[0]];
    }

    auto l = ax->plot(xEdge, yEdge, "-o");
    l->line_width(2);
    l->color("black");
    l->marker_size(10);
    l->marker_face_color("red");
  }

  f->show();
}


void exampleTimoshenkoBeam() {

  using namespace Dune::Functions::BasisFactory;
  using namespace Dune::Functions;
  using namespace Dune::Indices;
  const double b                   = 1;
  const double L                   = 5;
  const double E                   = 1;
  const double G                   = 1;
  const double t                   = 1;
  const double EI                  = E*b*t*t*t/12.0;
  const double GA                  = G*b*t;
  const double F                   = 0.1;
  Eigen::Matrix2d C;
  C << EI, 0,
       0, GA;
  const int numElements            = 10;
  const int numGP                  = 2;
  constexpr int polynomialOrderW   = 1;
  constexpr int polynomialOrderPhi = 1;
  Dune::OneDGrid grid(numElements, 0, L);
  auto gridView = grid.leafGridView();
  draw(gridView);

  // Basis with different orders for w (first) and phi (second)
  auto basis     = makeBasis(gridView, composite(lagrange<polynomialOrderW>(), lagrange<polynomialOrderPhi>(), FlatLexicographic()));
  auto basisView = basis.localView();

  // global stiffness matrix and force vector
  auto numDofs              = basis.size();
  Eigen::VectorXd F_ExtGlob = Eigen::VectorXd::Zero(numDofs);
  Eigen::MatrixXd K_Glob    = Eigen::MatrixXd::Zero(numDofs, numDofs);

  for (auto& ele : elements(gridView)) {
    basisView.bind(ele);

    // Define the integration rule
    const auto& rule = Dune::QuadratureRules<double, 1>::rule(ele.type(), numGP, Dune::QuadratureType::GaussLegendre);

    // get local stiffness matrix
    auto K_local = TimoshenkoBeamStiffness(basisView,ele,rule, C);

    // Adding local stiffness the global stiffness
    for (auto i = 0U; i < basisView.size(); ++i)
      for (auto j = 0U; j < basisView.size(); ++j)
        K_Glob(basisView.index(i)[0], basisView.index(j)[0]) += K_local(i,j);
  }

  // apply load on the right-hand side
  F_ExtGlob(getGlobalDofId(TimoschenkoBeam::w,basis,5)) = F;

  // clamp left-hand side
  std::vector<unsigned int> fixedDofs{getGlobalDofId(TimoschenkoBeam::w,basis,0.0),
                                     getGlobalDofId(TimoschenkoBeam::phi,basis,0.0)};
  for (auto dof : fixedDofs){
    for (unsigned int i = 0; i < numDofs; ++i) {
      K_Glob(i,dof) = 0.0;
      K_Glob(dof,i) = 0.0;
    }
    K_Glob(dof,dof) = 1.0;
  }

  // solve the linear system
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT);
  linSolver.factorize(K_Glob);
  const Eigen::VectorXd D_Glob = linSolver.solve(F_ExtGlob);

  // analytical solution
  std::cout << "Bernoulli solution for displacement at L: " << F*L*L*L/(3.0*EI) << "\n";

  // plot the result
  plotDeformedTimoschenkoBeam(basis,gridView,D_Glob);

}

  int main(){
  // exampleTrussElement();
  exampleTimoshenkoBeam();
}