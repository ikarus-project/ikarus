//
// Created by Alex on 21.07.2021.
//
#include <../../config.h>
#include <matplot/matplot.h>
#include <numbers>

#include <dune/common/indices.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/igaalgorithms.hh>
#include <dune/iga/nurbsgrid.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/Controlroutines/LoadControl.h"
#include "ikarus/LocalBasis/localBasis.h"
#include "ikarus/Solver/NonLinearSolver/NewtonRaphson.hpp"
#include "ikarus/utils/Observer/controlVTKWriter.h"
#include "ikarus/utils/Observer/nonLinearSolverLogger.h"
#include "ikarus/utils/utils/algorithms.h"
#include <ikarus/FiniteElements/AutodiffFE.h>
#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/LinearAlgebra/NonLinearOperator.h>
#include <ikarus/utils/concepts.h>

struct KirchhoffPlate {
  static constexpr double Emodul    = 2.1e8;
  static constexpr double nu        = 0.3;
  static constexpr double thickness = 0.1;

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

  template <typename LocalView, class Scalar>
  static Scalar calculateScalarImpl(const LocalView& localView, const Eigen::VectorXd& x,
                                    const Eigen::VectorX<Scalar>& dx, const double lambda) {
    const auto D  = constitutiveMatrix(Emodul, nu, thickness);
    Scalar energy = 0.0;
    auto& ele     = localView.element();
    auto& fe      = localView.tree().finiteElement();
    Eigen::VectorX<Scalar> wNodal;
    wNodal.setZero(fe.size());
    for (auto i = 0U; i < fe.size(); ++i)
      wNodal(i) = dx[localView.tree().localIndex(i)] + x[localView.index(localView.tree().localIndex(i))[0]];

    const auto& localBasis = fe.localBasis();

    const auto& rule = Dune::QuadratureRules<double, 2>::rule(ele.type(), 2 * localBasis.order());
    /// Calculate Kirchhoff plate energy
    for (auto& gp : rule) {
      //      std::vector<Dune::FieldMatrix<double, 1, 2>> dN_xi_eta;
      std::vector<Dune::FieldVector<double, 1>> dN_xixi;
      std::vector<Dune::FieldVector<double, 1>> dN_xieta;
      std::vector<Dune::FieldVector<double, 1>> dN_etaeta;
      std::vector<Dune::FieldVector<double, 1>> N_dune;
      Eigen::VectorXd N(fe.size());

      //      localBasis.evaluateJacobian(gp.position(), dN_xi_eta);
      localBasis.evaluateFunction(gp.position(), N_dune);
      std::ranges::copy(N_dune, N.begin());
      localBasis.partial({2, 0}, gp.position(), dN_xixi);
      localBasis.partial({1, 1}, gp.position(), dN_xieta);
      localBasis.partial({0, 2}, gp.position(), dN_etaeta);

      const auto Jinv
          = Ikarus::toEigenMatrix(ele.geometry().jacobianInverseTransposed(gp.position())).transpose().eval();
      //      std::cout << Jinv << std::endl;
      Eigen::VectorXd dN_xx(fe.size());
      Eigen::VectorXd dN_yy(fe.size());
      Eigen::VectorXd dN_xy(fe.size());
      for (auto i = 0U; i < fe.size(); ++i) {
        dN_xx[i] = dN_xixi[i] * Dune::power(Jinv(0, 0), 2) + dN_etaeta[i] * Dune::power(Jinv(0, 1), 2);
        dN_yy[i] = dN_etaeta[i] * Dune::power(Jinv(1, 1), 2) + dN_xixi[i] * Dune::power(Jinv(1, 0), 2);
        dN_xy[i] = (Jinv(1, 0) * dN_etaeta[i] + dN_xieta[i] * Jinv(0, 0)) * Jinv(1, 1)
                   + Jinv(0, 1) * (Jinv(1, 0) * dN_xieta[i] + dN_xixi[i] * Jinv(0, 0));
      }
      Eigen::Vector<Scalar, 3> kappa;
      kappa(0) = dN_xx.dot(wNodal);
      kappa(1) = dN_yy.dot(wNodal);
      kappa(2) = 2 * dN_xy.dot(wNodal);
      Scalar w = N.dot(wNodal);

      energy
          += (0.5 * kappa.dot(D * kappa) - w * lambda) * ele.geometry().integrationElement(gp.position()) * gp.weight();
    }

    /// Clamp boundary using penalty method
    const double penaltyFactor = 1e8;
    if (ele.hasBoundaryIntersections())
      for (auto& intersection : intersections(localView.globalBasis().gridView(), ele))
        if (intersection.boundary()) {
          const auto& rule1 = Dune::QuadratureRules<double, 1>::rule(intersection.type(), 2 * localBasis.order());
          for (auto& gp : rule1) {
            const auto& gpInElement = intersection.geometryInInside().global(gp.position());
            std::vector<Dune::FieldMatrix<double, 1, 2>> dN_xi_eta;
            localBasis.evaluateJacobian(gpInElement, dN_xi_eta);
            Eigen::VectorXd dN_x(fe.size());
            Eigen::VectorXd dN_y(fe.size());
            const auto Jinv
                = Ikarus::toEigenMatrix(ele.geometry().jacobianInverseTransposed(gpInElement)).transpose().eval();
            for (auto i = 0U; i < fe.size(); ++i) {
              dN_x[i] = dN_xi_eta[i][0][0] * Jinv(0, 0) + dN_xi_eta[i][0][1] * Jinv(0, 1);
              dN_y[i] = dN_xi_eta[i][0][0] * Jinv(1, 0) + dN_xi_eta[i][0][1] * Jinv(1, 1);
            }
            Scalar w_x = dN_x.dot(wNodal);
            Scalar w_y = dN_y.dot(wNodal);

            energy += 0.0*0.5 * penaltyFactor * (w_x * w_x + w_y * w_y);
          }
        }

    return energy;
  }
};

template <typename Basis>
class DenseFlatAssembler {
public:
  explicit DenseFlatAssembler(const Basis& basis, const std::vector<bool>& dirichFlags)
      : basis_{&basis}, dirichletFlags{&dirichFlags} {}

  Eigen::MatrixXd& getMatrix(const Eigen::VectorXd& displacement, const double& lambda) {
    return getMatrixImpl(displacement, lambda);
  }

  Eigen::VectorXd& getVector(const Eigen::VectorXd& displacement, const double& lambda) {
    return getVectorImpl(displacement, lambda);
  }

  double getScalar(const Eigen::VectorXd& displacement, const double& lambda) {
    return getScalarImpl(displacement, lambda);
  }

private:
  Eigen::MatrixXd& getMatrixImpl(const Eigen::VectorXd& displacement, const double& lambda) {
    mat.setZero(basis_->size(), basis_->size());
    auto localView = basis_->localView();
    for (auto& ge : elements(basis_->gridView())) {
      localView.bind(ge);
      auto matLoc = Ikarus::AutoDiffFE<KirchhoffPlate>::calculateMatrix(localView, displacement, lambda);
      for (auto i = 0U; i < localView.size(); ++i)
        for (auto j = 0U; j < localView.size(); ++j) {
          mat(localView.index(i)[0], localView.index(j)[0]) += matLoc(i, j);
        }
    }
    for (auto i = 0U; i < basis_->size(); ++i)
      if (dirichletFlags->at(i)) mat.col(i).setZero();
    for (auto i = 0U; i < basis_->size(); ++i)
      if (dirichletFlags->at(i)) mat.row(i).setZero();
    for (auto i = 0U; i < basis_->size(); ++i)
      if (dirichletFlags->at(i)) mat(i, i) = 1;
    return mat;
  }

  Eigen::VectorXd& getVectorImpl(const Eigen::VectorXd& displacement, const double& lambda) {
    vec.setZero(basis_->size());
    auto localView = basis_->localView();

    for (auto& ge : elements(basis_->gridView())) {
      localView.bind(ge);
      auto vecLocal = Ikarus::AutoDiffFE<KirchhoffPlate>::calculateVector(localView, displacement, lambda);
      for (auto i = 0U; i < localView.size(); ++i)
        vec(localView.index(i)[0]) += vecLocal(i);
    }
    for (auto i = 0U; i < basis_->size(); ++i)
      if (dirichletFlags->at(i)) vec[i] = 0;

    return vec;
  }

  double getScalarImpl(const Eigen::VectorXd& displacement, const double& lambda) {
    double scalar  = 0.0;
    auto localView = basis_->localView();

    for (auto& ge : elements(basis_->gridView())) {
      localView.bind(ge);
      for (auto i = 0U; i < localView.size(); ++i)
        scalar += Ikarus::AutoDiffFE<KirchhoffPlate>::calculateScalar(localView, displacement, lambda);
    }

    return scalar;
  }
  Basis const* basis_;
  std::vector<bool> const* dirichletFlags;
  Eigen::MatrixXd mat{};
  Eigen::VectorXd vec{};
};

int main() {
  ///Create IGA Grid
  using namespace Ikarus;
  constexpr int griddim                                    = 2;
  constexpr int dimworld                                   = 2;
  const std::array<std::vector<double>, griddim> knotSpans = {{{0, 0, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<griddim, dimworld>::ControlPointType;

  const double Lx = 1;
  const double Ly = 1;
  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0}, .w = 1}, {.p = {Lx, 0}, .w = 1}}, {{.p = {0, Ly}, .w = 1}, {.p = {Lx, Ly}, .w = 1}}};

  std::array<int, griddim> dimsize = {2, 2};

  std::vector<double> dofsVec;
  std::vector<double> l2Evcector;

    auto controlNet = Dune::IGA::NURBSPatchData<griddim, dimworld>::ControlPointNetType(dimsize, controlPoints);
    using Grid      = Dune::IGA::NURBSGrid<griddim, dimworld>;

    Dune::IGA::NURBSPatchData<griddim, dimworld> patchData;
    patchData.knotSpans     = knotSpans;
    patchData.degree        = {1, 1};
    patchData.controlPoints = controlNet;
    patchData               = Dune::IGA::degreeElevate(patchData, 0, 1);
    patchData               = Dune::IGA::degreeElevate(patchData, 1, 1);
    Grid grid(patchData);
    grid.globalRefine(2);
    auto gridView = grid.leafGridView();
    draw(gridView);
    using namespace Dune::Functions::BasisFactory;
    auto basis = makeBasis(gridView, gridView.getPreBasis());
    std::vector<double> w(basis.size());

    auto dirichletPredicate = [](auto p) { return std::sin(p[0]); };
    std::vector<double> uhat(basis.size());

    std::vector<bool> dirichletFlags(basis.size());
    std::fill(dirichletFlags.begin(), dirichletFlags.end(), false);

    Dune::Functions::forEachBoundaryDOF(basis, [&](auto&& index) { dirichletFlags[index] = true; });
    auto denseAssembler = DenseFlatAssembler(basis, dirichletFlags);

    Eigen::VectorXd d;
    d.setZero(basis.size());
    double lambda = 0.0;

    auto energyFunction = [&](auto&& lambda, auto&& disp) -> auto { return denseAssembler.getScalar(disp, lambda); };
    auto fintFunction   = [&](auto&& lambda, auto&& disp) -> auto& { return denseAssembler.getVector(disp, lambda); };
    auto KFunction      = [&](auto&& lambda, auto&& disp) -> auto& { return denseAssembler.getMatrix(disp, lambda); };

    auto nonLinOp  = Ikarus::NonLinearOperator(linearAlgebraFunctions(energyFunction, fintFunction, KFunction),
                                               parameter(lambda, d));
    auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT);

    auto nr                      = Ikarus::NewtonRaphson(nonLinOp.subOperator<1, 2>(), std::move(linSolver));
    auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();

    auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<decltype(basis)>>(basis, d, 2);
    vtkWriter->setFileNamePrefix("TestKLplate");
    vtkWriter->setVertexSolutionName("displacement");
    nr.subscribeAll(nonLinearSolverObserver);

    const double totalLoad = 2000;
    auto lc                = Ikarus::LoadControl(std::move(nr), 20, {0, totalLoad});

    lc.subscribe(ControlMessages::SOLUTION_CHANGED,vtkWriter);
    std::cout << "Energy before: " << nonLinOp.value() << std::endl;
    lc.run();
    nonLinOp.update<0>();
    std::cout << "Energy after: " << nonLinOp.value() << std::endl;

    const double D = KirchhoffPlate::Emodul * Dune::power(KirchhoffPlate::thickness, 3)
                     / (12 * (1 - Dune::power(KirchhoffPlate::nu, 2)));
    auto wxy =
        [&](auto x,
            auto
                y) {  // https://en.wikipedia.org/wiki/Bending_of_plates#Simply-supported_plate_with_uniformly-distributed_load
          double w                = 0.0;
          const int seriesFactors = 40;
          const double pi         = std::numbers::pi;
          auto isOdd              = std::views::filter([](auto i) { return i % 2 != 0; });
          auto oddFactors         = std::ranges::iota_view(1, seriesFactors) | isOdd;
          for (auto m : oddFactors) {
            for (auto n : oddFactors) {
              w += std::sin(m * pi * x / Lx) * std::sin(n * pi * y / Ly)
                   / (m * n * Dune::power(m * m / (Lx * Lx) + n * n / (Ly * Ly), 2));
            }
          }
          return 16 * totalLoad / (Dune::power(pi, 6) * D) * w;
        };
//    std::cout << wxy(Lx / 2.0, Ly / 2.0) << std::endl;

    const double wCenterClamped = 1.265319087
                                  / (D / (totalLoad * Dune::power(Lx, 4))
                                     * 1000.0);  // clamped sol http://faculty.ce.berkeley.edu/rlt/reports/clamp.pdf
//    std::cout << wCenterClamped << std::endl;
    auto disp      = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 1>>(basis, d);
    auto localDisp = localFunction(disp);
    auto localView = basis.localView();

    double l2_error = 0.0;
    for (auto& ele : elements(gridView)) {
      localDisp.bind(ele);
      localView.bind(ele);
      const auto geo   = localView.element().geometry();
      const auto& rule = Dune::QuadratureRules<double, 2>::rule(
          ele.type(), 2 * localView.tree().finiteElement().localBasis().order());
      for (auto gp : rule) {
        const auto gpGlobalPos = geo.global(gp.position());
        const auto w_ex        = wxy(gpGlobalPos[0], gpGlobalPos[1]);
        const auto w_fe        = localDisp(gp.position());
        l2_error += Dune::power(w_ex - w_fe, 2) * ele.geometry().integrationElement(gp.position()) * gp.weight();
      }
    }
    l2_error = std::sqrt(l2_error);
    std::cout << "l2_error: " << l2_error << "Dofs:: " << basis.size() << std::endl;
    dofsVec.push_back(basis.size());
    l2Evcector.push_back(l2_error);

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