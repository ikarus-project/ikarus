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
  constexpr int griddim                                    = 2;
  const double L = 1;
  Dune::OneDGrid grid(10, 0, L);
  grid.globalRefine(2);
  auto gridView = grid.leafGridView();
  draw(gridView);
  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, lagrange<1>());
  std::vector<double> w(basis.size());
  auto localView = basis.localView();

  auto dirichletPredicate = [](auto p) { return std::sin(p[0]); };
  std::vector<double> uhat(basis.size());
  
//
//  for (auto& ele : elements(gridView)) {
//    localView.bind(ele);
//    auto& fe = localView.tree().finiteElement();
//    Ikarus::LocalBasis localBasis(fe.localBasis());
//    const auto& rule = Dune::QuadratureRules<double, 2>::rule(ele.type(), 2 * fe.order());
//    for (auto& gp : rule) {
//      localBasis.evaluateFunctionAndJacobian()
//
//          energy
//          += localBasis.
//    }
//  }
}