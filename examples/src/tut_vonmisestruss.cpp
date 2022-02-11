//
// Created by Alex on 21.07.2021.
//

#include <../../config.h>
#include <autodiff/forward/dual/dual.hpp>

#include <matplot/matplot.h>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/FiniteElements/AutodiffFE.h>
#include <ikarus/Grids/GridHelper/griddrawer.h>

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

struct Truss {
  constexpr static const double EA = 100;
  template <typename LocalView, class Scalar>
  static Scalar calculateScalarImpl(const LocalView& localView, const Eigen::VectorXd& d,
                                    const Eigen::Vector4<Scalar>& dx) {
    Scalar energy = 0.0;
    auto& ele     = localView.element();
    const auto X1 = Ikarus::toEigenVector(ele.geometry().corner(0));
    const auto X2 = Ikarus::toEigenVector(ele.geometry().corner(1));

    Eigen::Matrix<Scalar, 2, 2> u;
    u.setZero();
    for (int i = 0; i < 2; ++i)
      for (int k2 = 0; k2 < 2; ++k2)
        u.col(i)(k2) = dx[localView.tree().child(k2).localIndex(i)]
                       + d[localView.index(localView.tree().child(k2).localIndex(i))[0]];

    const Eigen::Vector2<Scalar> x1 = X1 + u.col(0);
    const Eigen::Vector2<Scalar> x2 = X2 + u.col(1);

    const double LRefsquared = (X1 - X2).squaredNorm();
    const Scalar lsquared    = (x1 - x2).squaredNorm();

    const Scalar Egl = 1.0 / 2.0 * (lsquared - LRefsquared) / LRefsquared;

    energy = 1.0 / 2.0 * EA/sqrt(LRefsquared) * Egl * Egl ;
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
  Eigen::VectorXd u(basis.size());
  u.setZero();
  std::cout << "Size: " << u.size() << std::endl;
  Eigen::Matrix<double, 6, 6> K;
  Eigen::Matrix<double, 6, 1> R;

  auto localView          = basis.localView();
  auto stiffnessAndForces = [&](auto& uL, const double lambda) {
    K.setZero();
    R.setZero();
    for (auto& ele : elements(gridView)) {
      localView.bind(ele);
      const auto FintEle = Ikarus::AutoDiffFE<Truss>::calculateVector(localView, uL);
      const auto Kele    = Ikarus::AutoDiffFE<Truss>::calculateMatrix(localView, uL);
      for (auto i = 0U; i < localView.size(); ++i) {
        R(localView.index(i)[0]) += FintEle(i);
        for (auto j = 0U; j < localView.size(); ++j)
          K(localView.index(i)[0], localView.index(j)[0]) += Kele(i, j);
      }
    }
    Eigen::Matrix<double, 2, 2> Kred = K({2, 3}, {2, 3});
    Eigen::Vector<double, 2> Rred    = R({2, 3});
    Rred[1] -= -lambda;
    return std::make_tuple(Kred, Rred);
  };

  double lambla = 0;

  const auto [Kred_, Rred_] = stiffnessAndForces(u, lambla);

  std::cout << Kred_ << std::endl;
  std::cout << Rred_ << std::endl;
  std::cout << Kred_.fullPivLu().rank() << std::endl;

  const int loadSteps       = 60;
  const double lambdaFactor = 0.5;
  Eigen::Matrix3Xd lambdaAndDisp;
  lambdaAndDisp.setZero(Eigen::NoChange, loadSteps);
  for (int ls = 0; ls < loadSteps; ++ls) {
    std::cout << "====" << ls << "====" << std::endl;
    for (int iter = 0; iter < 100; ++iter) {
      const auto [Kred, Rred] = stiffnessAndForces(u, lambla);
      u.segment<2>(2) -= Kred.inverse() * Rred;
      std::cout << " Rnorm: " << Rred.norm() << std::endl;
      if (Rred.norm() < 1e-8) break;
    }
    lambdaAndDisp(0, ls) = lambla;
    lambdaAndDisp(1, ls) = u[2];
    lambdaAndDisp(2, ls) = u[3];
    lambla += lambdaFactor;

    auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,2>>(basis,u);
    Dune::VTKWriter vtkWriter(gridView);
    vtkWriter.addVertexData(disp, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
    vtkWriter.write("TestTruss" + std::to_string(ls));
  }

  using namespace matplot;
  Eigen::VectorXd lambdaVec = lambdaAndDisp.row(0);
  Eigen::VectorXd dVec      = -lambdaAndDisp.row(2);
  auto f                    = figure(true);

  title("Load-DisplacementCurve");
  xlabel("y-Displacement");
  ylabel("LoadFactor");

  auto analyticalLoadDisplacementCurve = [&](auto& w) {
    const double Ltruss = std::sqrt(h*h + L*L);
    return Truss::EA * Dune::power(h, 3) / Dune::power(Ltruss, 3)
           * (w / h - 1.5 * Dune::power(w / h, 2) + 0.5 * Dune::power(w / h, 3));
  };

  std::vector<double> x = linspace(0.0,dVec.maxCoeff());
  std::vector<double> y1 = transform(x, [&](auto x) { return analyticalLoadDisplacementCurve(x); });
  auto p = plot(x, y1,dVec, lambdaVec);
  p[0]->line_width(2);
  p[1]->line_width(2);
  p[1]->marker(line_spec::marker_style::asterisk);
  show();
}