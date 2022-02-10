//
// Created by Alex on 21.07.2021.
//
#include <../../config.h>
#include <numbers>

#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/igaalgorithms.hh>
#include <dune/iga/nurbsgrid.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/LocalBasis/localBasis.h"
#include <ikarus/Grids/GridHelper/griddrawer.h>

int main() {
  constexpr int griddim                                    = 2;
  constexpr int dimworld                                   = 2;
  const std::array<std::vector<double>, griddim> knotSpans = {{{0, 0, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<griddim, dimworld>::ControlPointType;

  const double Lx = 1;
  const double Ly = 1;
  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0}, .w = 1}, {.p = {Lx, 0}, .w = 1}}, {{.p = {0, Ly}, .w = 1}, {.p = {Lx, Ly}, .w = 1}}};

  std::array<int, griddim> dimsize = {2, 2};

  auto controlNet = Dune::IGA::NURBSPatchData<griddim, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<griddim, dimworld>;

  Dune::IGA::NURBSPatchData<griddim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = {1, 1};
  patchData.controlPoints = controlNet;
  patchData               = Dune::IGA::degreeElevate(patchData, 1, 0);
  patchData               = Dune::IGA::degreeElevate(patchData, 1, 1);
  Grid grid(patchData);
  grid.globalRefine(2);
  auto gridView = grid.leafGridView();
  draw(gridView);
  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, gridView.getPreBasis());
  std::vector<double> w(basis.size());
  auto localView = basis.localView();

  auto dirichletPredicate = [](auto p) { return std::sin(p[0]); };
  std::vector<double> uhat(basis.size());

  Dune::Functions::interpolate(basis, uhat, dirichletPredicate);
  for (auto& ele : elements(gridView)) {
    localView.bind(ele);
    auto& fe = localView.tree().finiteElement();
    Ikarus::LocalBasis localBasis(fe.localBasis());
    const auto& rule = Dune::QuadratureRules<double, 2>::rule(ele.type(), 2 * fe.order());
    for (auto& gp : rule) {
      localBasis.evaluateFunctionAndJacobian()

          energy
          += localBasis.
    }
  }
}