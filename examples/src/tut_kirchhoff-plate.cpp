/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#include <config.h>

#include <matplot/matplot.h>
#include <numbers>

#include <dune/common/indices.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/igaalgorithms.hh>
#include <dune/iga/nurbsgrid.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlRoutines/loadControl.hh>
#include <ikarus/finiteElements/feBases/autodiffFE.hh>
#include <ikarus/finiteElements/feBases/scalarFE.hh>
#include <ikarus/finiteElements/feTraits.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/utils/algorithms.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/eigenDuneTransformations.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <ikarus/utils/observer/loadControlObserver.hh>
#include <ikarus/utils/observer/nonLinearSolverLogger.hh>

template <typename Basis>
struct KirchhoffPlate : Ikarus::ScalarFieldFE<Basis>, Ikarus::AutoDiffFE<KirchhoffPlate<Basis>, Basis> {
  using BaseDisp = Ikarus::ScalarFieldFE<Basis>;
  using BaseAD   = Ikarus::AutoDiffFE<KirchhoffPlate<Basis>, Basis>;
  using BaseAD::size;
  using LocalView         = typename Basis::LocalView;
  using FERequirementType = typename BaseAD::FERequirementType;

  KirchhoffPlate(const Basis& basis, const typename LocalView::Element& element, double p_Emodul, double p_nu,
                 double p_thickness)
      : BaseDisp(basis, element),
        BaseAD(basis, element),
        localView_{basis.localView()},
        Emodul{p_Emodul},
        nu{p_nu},
        thickness{p_thickness} {
    localView_.bind(element);
    geometry_ = localView_.element().geometry();
  }

  static Eigen::Matrix<double, 3, 3> constitutiveMatrix(double Emod, double p_nu, double p_thickness) {
    const double factor = Emod * Dune::power(p_thickness, 3) / (12.0 * (1.0 - p_nu * p_nu));
    Eigen::Matrix<double, 3, 3> D;
    D.setZero();
    D(0, 0) = 1;
    D(0, 1) = D(1, 0) = p_nu;
    D(1, 1)           = 1;
    D(2, 2)           = (1 - p_nu) / 2.0;
    D *= factor;
    return D;
  }

  template <class Scalar>
  [[nodiscard]] Scalar calculateScalarImpl(const FERequirementType& par, const Eigen::VectorX<Scalar>& dx) const {
    const auto& wGlobal = par.getSolution(Ikarus::FESolutions::displacement);
    const auto& lambda  = par.getParameter(Ikarus::FEParameter::loadfactor);
    const auto D        = constitutiveMatrix(Emodul, nu, thickness);
    Scalar energy       = 0.0;
    auto& ele           = localView_.element();
    auto& fe            = localView_.tree().finiteElement();
    Eigen::VectorX<Scalar> wNodal;
    wNodal.setZero(fe.size());
    for (auto i = 0U; i < fe.size(); ++i)
      wNodal(i) = dx[i] + wGlobal[localView_.index(localView_.tree().localIndex(i))[0]];

    const auto& localBasis = fe.localBasis();

    const auto& rule = Dune::QuadratureRules<double, 2>::rule(ele.type(), 2 * localBasis.order());
    /// Calculate Kirchhoff plate energy
    for (auto& gp : rule) {
      std::vector<Dune::FieldVector<double, 1>> dN_xixi;
      std::vector<Dune::FieldVector<double, 1>> dN_xieta;
      std::vector<Dune::FieldVector<double, 1>> dN_etaeta;
      std::vector<Dune::FieldVector<double, 1>> N_dune;
      Eigen::VectorXd N(fe.size());

      localBasis.evaluateFunction(gp.position(), N_dune);
      std::ranges::copy(N_dune, N.begin());
      localBasis.partial({2, 0}, gp.position(), dN_xixi);
      localBasis.partial({1, 1}, gp.position(), dN_xieta);
      localBasis.partial({0, 2}, gp.position(), dN_etaeta);

      const auto Jinv = Ikarus::toEigenMatrix(geometry_.jacobianInverseTransposed(gp.position())).transpose().eval();

      Eigen::VectorXd dN_xx(fe.size());
      Eigen::VectorXd dN_yy(fe.size());
      Eigen::VectorXd dN_xy(fe.size());
      using Dune::power;
      for (auto i = 0U; i < fe.size(); ++i) {
        dN_xx[i] = dN_xixi[i] * power(Jinv(0, 0), 2);
        dN_yy[i] = dN_etaeta[i] * power(Jinv(1, 1), 2);
        dN_xy[i] = dN_xieta[i] * Jinv(0, 0) * Jinv(1, 1);
      }
      Eigen::Vector<Scalar, 3> kappa;
      kappa(0) = dN_xx.dot(wNodal);
      kappa(1) = dN_yy.dot(wNodal);
      kappa(2) = 2 * dN_xy.dot(wNodal);
      Scalar w = N.dot(wNodal);

      energy += (0.5 * kappa.dot(D * kappa) - w * lambda) * geometry_.integrationElement(gp.position()) * gp.weight();
    }

    /// Clamp boundary using penalty method
    const double penaltyFactor = 1e8;
    if (ele.hasBoundaryIntersections())
      for (auto& intersection : intersections(localView_.globalBasis().gridView(), ele))
        if (intersection.boundary()) {
          const auto& rule1 = Dune::QuadratureRules<double, 1>::rule(intersection.type(), 2 * localBasis.order());
          for (auto& gp : rule1) {
            const auto& gpInElement = intersection.geometryInInside().global(gp.position());
            std::vector<Dune::FieldMatrix<double, 1, 2>> dN_xi_eta;
            localBasis.evaluateJacobian(gpInElement, dN_xi_eta);
            Eigen::VectorXd dN_x(fe.size());
            Eigen::VectorXd dN_y(fe.size());
            const auto Jinv
                = Ikarus::toEigenMatrix(geometry_.jacobianInverseTransposed(gpInElement)).transpose().eval();
            for (auto i = 0U; i < fe.size(); ++i) {
              dN_x[i] = dN_xi_eta[i][0][0] * Jinv(0, 0);
              dN_y[i] = dN_xi_eta[i][0][1] * Jinv(1, 1);
            }
            const Scalar w_x = dN_x.dot(wNodal);
            const Scalar w_y = dN_y.dot(wNodal);

            energy += 0.0 * 0.5 * penaltyFactor * (w_x * w_x + w_y * w_y);
          }
        }

    return energy;
  }

private:
  LocalView localView_;
  typename LocalView::Element::Geometry geometry_;
  double Emodul;
  double nu;
  double thickness;
};

int main() {
  /// Create 2D nurbs grid
  using namespace Ikarus;
  constexpr int griddim                                    = 2;
  constexpr int dimworld                                   = 2;
  const std::array<std::vector<double>, griddim> knotSpans = {{{0, 0, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<griddim, dimworld>::ControlPointType;

  const double Lx = 1;
  const double Ly = 1;
  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0}, .w = 1}, {.p = {0, Ly}, .w = 1}}, {{.p = {Lx, 0}, .w = 1}, {.p = {Lx, Ly}, .w = 1}}};

  std::array<int, griddim> dimsize = {2, 2};

  std::vector<double> dofsVec;
  std::vector<double> l2Evcector;
  auto controlNet = Dune::IGA::NURBSPatchData<griddim, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<griddim, dimworld>;

  Dune::IGA::NURBSPatchData<griddim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = {1, 1};
  patchData.controlPoints = controlNet;
  /// Increate polynomial degree in each direction
  patchData = Dune::IGA::degreeElevate(patchData, 0, 1);
  patchData = Dune::IGA::degreeElevate(patchData, 1, 1);
  Grid grid(patchData);

  for (int ref = 0; ref < 5; ++ref) {
    auto gridView = grid.leafGridView();
    //    draw(gridView);
    using namespace Dune::Functions::BasisFactory;
    /// Create nurbs basis with extracted preBase from grid
    auto basis = makeBasis(gridView, gridView.getPreBasis());
    /// Fix complete boundary (simply supported plate)
    std::vector<bool> dirichletFlags(basis.size(), false);
    Dune::Functions::forEachBoundaryDOF(basis, [&](auto&& index) { dirichletFlags[index] = true; });

    /// Create finite elements
    auto localView         = basis.localView();
    const double Emod      = 2.1e8;
    const double nu        = 0.3;
    const double thickness = 0.1;
    std::vector<KirchhoffPlate<decltype(basis)>> fes;
    for (auto& ele : elements(gridView))
      fes.emplace_back(basis, ele, Emod, nu, thickness);

    /// Create assembler
    auto denseAssembler = DenseFlatAssembler(basis, fes, dirichletFlags);

    /// Create non-linear operator with potential energy
    Eigen::VectorXd w;
    w.setZero(basis.size());

    const double totalLoad = 2000;

    auto kFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
      Ikarus::FErequirements req = FErequirementsBuilder()
                                       .insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
                                       .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                       .addAffordance(Ikarus::MatrixAffordances::stiffness)
                                       .build();
      return denseAssembler.getMatrix(req);
    };

    auto rFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
      Ikarus::FErequirements req = FErequirementsBuilder()
                                       .insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
                                       .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                       .addAffordance(Ikarus::VectorAffordances::forces)
                                       .build();
      return denseAssembler.getVector(req);
    };

    const auto& K = kFunction(w, totalLoad);
    const auto& R = rFunction(w, totalLoad);
    Eigen::LDLT<Eigen::MatrixXd> solver;
    solver.compute(K);
    w -= solver.solve(R);

    // Output solution to vtk
    auto wGlobalFunc = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(basis, w);
    Dune::SubsamplingVTKWriter vtkWriter(gridView, Dune::refinementLevels(2));
    vtkWriter.addVertexData(wGlobalFunc, Dune::VTK::FieldInfo("w", Dune::VTK::FieldInfo::Type::scalar, 1));
    vtkWriter.write("Test_KPlate");

    /// Create analytical solution function for the simply supported case
    const double D = Emod * Dune::power(thickness, 3) / (12 * (1 - Dune::power(nu, 2)));
    // https://en.wikipedia.org/wiki/Bending_of_plates#Simply-supported_plate_with_uniformly-distributed_load
    auto wAna = [&](auto x) {
      double w                = 0.0;
      const int seriesFactors = 40;
      const double pi         = std::numbers::pi;
      auto oddFactors
          = std::ranges::iota_view(1, seriesFactors) | std::views::filter([](auto i) { return i % 2 != 0; });
      for (auto m : oddFactors)
        for (auto n : oddFactors)
          w += sin(m * pi * x[0] / Lx) * sin(n * pi * x[1] / Ly)
               / (m * n * Dune::power(m * m / (Lx * Lx) + n * n / (Ly * Ly), 2));

      return 16 * totalLoad / (Dune::power(pi, 6) * D) * w;
    };
    //    std::cout << wxy(Lx / 2.0, Ly / 2.0) << std::endl;

    /// Displacement at center of clamped square plate
    // clamped sol http://faculty.ce.berkeley.edu/rlt/reports/clamp.pdf
    const double wCenterClamped = 1.265319087 / (D / (totalLoad * Dune::power(Lx, 4)) * 1000.0);
    //    std::cout << wCenterClamped << std::endl;
    auto wGlobalFunction = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 1>>(basis, w);
    auto wGlobalAnalyticFunction = Dune::Functions::makeAnalyticGridViewFunction(wAna, gridView);
    auto localw                  = localFunction(wGlobalFunction);
    auto localwAna               = localFunction(wGlobalAnalyticFunction);

    /// Calculate L_2 error for simply supported case
    double l2_error = 0.0;
    for (auto& ele : elements(gridView)) {
      localView.bind(ele);
      localw.bind(ele);
      localwAna.bind(ele);
      const auto geo   = localView.element().geometry();
      const auto& rule = Dune::QuadratureRules<double, 2>::rule(
          ele.type(), 2 * localView.tree().finiteElement().localBasis().order());
      for (auto gp : rule) {
        const auto gpGlobalPos = geo.global(gp.position());

        const auto w_ex = localwAna(gp.position());
        const auto w_fe = localw(gp.position());
        l2_error += Dune::power(w_ex - w_fe, 2) * ele.geometry().integrationElement(gp.position()) * gp.weight();
      }
    }

    l2_error = std::sqrt(l2_error);
    std::cout << "l2_error: " << l2_error << " Dofs:: " << basis.size() << std::endl;
    dofsVec.push_back(basis.size());
    l2Evcector.push_back(l2_error);
    grid.globalRefine(1);
  }
  /// Draw L_2 error over dofs count
  using namespace matplot;
  auto f  = figure(true);
  auto ax = gca();
  ax->y_axis().label("L2_error");

  ax->x_axis().label("#Dofs");
  auto p = ax->loglog(dofsVec, l2Evcector);
  p->line_width(2);
  p->marker(line_spec::marker_style::asterisk);
  show();
}