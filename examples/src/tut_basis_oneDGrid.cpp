//
// Created by Alex on 21.07.2021.
//
#include <../../config.h>
#include <numbers>


#include <dune/geometry/quadraturerules.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/grid/onedgrid.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/LocalBasis/localBasis.h"
#include <ikarus/Grids/GridHelper/griddrawer.h>




int main(){
  constexpr int griddim = 1;
  const double L = 1;
  Dune::OneDGrid grid(10, 0, L);
  //grid.globalRefine(2);
  auto gridView = grid.leafGridView();
  std::cout << gridView.dimensionworld << "\n";
  draw(gridView);
  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, lagrange<1>());
  std::vector<double> w(basis.size());
  auto basisView = basis.localView();


  for (auto& ele : elements(gridView)) {
    basisView.bind(ele);
    auto& fe = basisView.tree().finiteElement();
    Ikarus::LocalBasis localBasis(fe.localBasis());
    const auto& rule = Dune::QuadratureRules<double, 1>::rule(ele.type(),2,Dune::QuadratureType::GaussLegendre);
//    for (auto& gp : rule) {
//      std::cout << gp.position() << "\n";
//    }
//    for (auto& gp : rule) {
//      localBasis.evaluateFunctionAndJacobian()
//
//          energy
//          += localBasis.
//    }
  }
}