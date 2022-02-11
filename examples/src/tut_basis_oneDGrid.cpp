//
// Created by Alex on 21.07.2021.
//
#include <../../config.h>
#include <numbers>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/onedgrid.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/LocalBasis/localBasis.h"
#include <ikarus/Grids/GridHelper/griddrawer.h>

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


void exampleTimoschenkoBeam(){
  constexpr int griddim = 1;
  const double L = 1;
  const double EA = 1.0;
  const int numElements = 10;
  const int numGP = 2;
  constexpr int polynomialOrderW = 1;
  constexpr int polynomialOrderPhi = 1;
  Dune::OneDGrid grid(numElements, 0, L);
  //grid.globalRefine(2);
  auto gridView = grid.leafGridView();
  std::cout << gridView.dimensionworld << "\n";
  draw(gridView);
  using namespace Dune::Functions::BasisFactory;
  // ToDo Basis mit unterschiedlichen Ordnungen
  auto basis = makeBasis(gridView, power<1>(lagrange<1>(), FlatInterleaved()));
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



int main(){
  // exampleTrussElement();
  exampleTimoschenkoBeam();
}