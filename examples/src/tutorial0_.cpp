//
// Created by Alex on 21.07.2021.
//

#include <../../config.h>
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <numbers>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/foamgrid/foamgrid/foamgridfactory.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/igaalgorithms.hh>
#include <dune/iga/nurbsgrid.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/LocalBasis/localBasis.h"
#include <ikarus/Grids/GridHelper/griddrawer.h>

struct Truss {
  constexpr static const double EA = 1000;

  template <typename LocalView>
  static double calculateEnergy(const LocalView& localView, const std::vector<double>& d) {
    return calculateEnergyImpl<double>(localView, d);
  }

  template <typename LocalView>
  static auto stiffnessMatrix(const LocalView& localView, const std::vector<double>& d) {
    Eigen::Vector4dual2nd dx(localView.size());
    dx.setZero();
    auto f = [&](auto& x) { return calculateEnergyImpl(localView, d, x); };
    return hessian(f, wrt(dx), at(dx));
  }

  template <typename LocalView>
  static auto internalForces(const LocalView& localView, const std::vector<double>& d) {
    Eigen::Vector4dual dx(localView.size());
    dx.setZero();
    auto f = [&](auto& x) { return calculateEnergyImpl(localView, d, x); };
    return gradient(f, wrt(dx), at(dx));
  }

private:
  template <typename LocalView, class Scalar>
  static Scalar calculateEnergyImpl(const LocalView& localView, const std::vector<double>& d,
                                    const Eigen::Vector4<Scalar>& dx) {
    Scalar energy = 0.0;
    auto& ele     = localView.element();
    const auto X1 = Ikarus::toEigenVector(ele.geometry().corner(0));
    const auto X2 = Ikarus::toEigenVector(ele.geometry().corner(1));

    Eigen::Matrix<Scalar, 2, 2> u;
    u.setZero();
    for (int i = 0; i < 2; ++i)
      for (int k2 = 0; k2 < 2; ++k2) {
        //        std::cout<<"i,k2: "<<i<<" "<<k2<<" "<<localView.tree().child(k2).localIndex(i)<<std::endl;
        u.col(i)(k2) = dx[localView.tree().child(k2).localIndex(i)]
                       + d[localView.index(localView.tree().child(k2).localIndex(i))[0]];
      }

    const Eigen::Vector2<Scalar> x1 = X1 + u.col(0);
    const Eigen::Vector2<Scalar> x2 = X2 + u.col(1);

    const double LRefsquared = (X1 - X2).squaredNorm();
    const Scalar lsquared    = (x1 - x2).squaredNorm();

    const Scalar Egl = 1.0 / 2.0 * (lsquared - LRefsquared) / LRefsquared;

    energy = 1.0 / 2.0 * EA * Egl * Egl * sqrt(LRefsquared);
    return energy;
  }
};

int main() {
  Dune::GridFactory<Dune::FoamGrid<1, 2, double>> gridFactory;
  const double h = 1.0;
  const double L = 1.0;
  gridFactory.insertVertex({0, 0});
  gridFactory.insertVertex({L, h});
  gridFactory.insertVertex({2 * L, 0});
  gridFactory.insertElement(Dune::GeometryTypes::line, {0, 1});
  gridFactory.insertElement(Dune::GeometryTypes::line, {1, 2});
  auto grid     = gridFactory.createGrid();
  auto gridView = grid->leafGridView();
  draw(gridView);

  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));
  std::vector<double> u(basis.size());
  std::fill(u.begin(), u.end(), 0.0);
  std::cout << "Size: " << u.size() << std::endl;
  Eigen::Matrix<double, 6, 6> K;
  Eigen::Matrix<double, 6, 1> R;
  K.setZero();
  R.setZero();
  auto localView = basis.localView();

  for (auto& ele : elements(gridView)) {
    localView.bind(ele);
    const auto FintEle = Truss::internalForces(localView, u);
    const auto Kele    = Truss::stiffnessMatrix(localView, u);
    std::vector<double> indices;
    for (auto i = 0U; i < localView.size(); ++i) {
      R(localView.index(i)[0]) += FintEle(i);
      for (auto j = 0U; j < localView.size(); ++j) {
        K(localView.index(i)[0], localView.index(j)[0]) += Kele(i, j);
      }
    }
  }

  std::cout << K << std::endl;
  std::cout << K.fullPivLu().rank() << std::endl;
  std::cout << R << std::endl;

  Eigen::Matrix<double, 2, 2> Kred = K({2, 3}, {2, 3});
  std::cout << Kred << std::endl;
  std::cout << Kred.fullPivLu().rank() << std::endl;
}