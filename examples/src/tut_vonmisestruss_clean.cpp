//
// Created by Alex on 21.07.2021.
//

#include <../../config.h>
#include <autodiff/forward/dual/dual.hpp>
#include <matplot/matplot.h>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/Controlroutines/LoadControl.h"
#include "ikarus/LinearAlgebra/NonLinearOperator.h"
#include "ikarus/Solver/LinearSolver/LinearSolver.h"
#include "ikarus/Solver/NonLinearSolver/NewtonRaphson.hpp"
#include "ikarus/utils/Observer/controlVTKWriter.h"
#include "ikarus/utils/Observer/genericControlObserver.h"
#include "ikarus/utils/Observer/nonLinearSolverLogger.h"
#include <ikarus/Assembler/SimpleAssemblers.h>
#include <ikarus/FiniteElements/AutodiffFE.h>
#include <ikarus/FiniteElements/FEPolicies.h>
#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/Variables/ParameterFactory.h>

using namespace Ikarus;
template <typename Basis>
struct Truss : Ikarus::FiniteElements::FEDisplacement<Basis>, Ikarus::AutoDiffFEClean<Truss<Basis>, Basis> {
  using BaseDisp = Ikarus::FiniteElements::FEDisplacement<Basis>;
  using BaseAD   = Ikarus::AutoDiffFEClean<Truss<Basis>, Basis>;
  using BaseAD::size;
  friend BaseAD;
  using LocalView         = typename Basis::LocalView;
  using FERequirementType = typename BaseAD::FERequirementType;
  using Traits            = TraitsFromLocalView<LocalView>;
  Truss(const Basis& basis, const typename LocalView::Element& element, double p_EA)
      : BaseDisp(basis, element), BaseAD(basis, element), localView_{basis.localView()}, EA{p_EA} {
    localView_.bind(element);
  }

private:
  template <class Scalar>
  Scalar calculateScalarImpl(const FERequirementType& par, const Eigen::VectorX<Scalar>& dx) const {
    const auto& d      = par.sols[0].get();
    const auto& lambda = par.parameter.at(FEParameter::loadfactor);

    auto& ele     = localView_.element();
    const auto X1 = Ikarus::toEigenVector(ele.geometry().corner(0));
    const auto X2 = Ikarus::toEigenVector(ele.geometry().corner(1));

    Eigen::Matrix<Scalar, Traits::worlddim, 2> u;
    u.setZero();
    for (int i = 0; i < 2; ++i)
      for (int k2 = 0; k2 < Traits::worlddim; ++k2)
        u.col(i)(k2)
            = dx[Traits::worlddim * i + k2] + d[localView_.index(localView_.tree().child(k2).localIndex(i))[0]];

    const Eigen::Vector2<Scalar> x1 = X1 + u.col(0);
    const Eigen::Vector2<Scalar> x2 = X2 + u.col(1);

    const double LRefsquared = (X1 - X2).squaredNorm();
    const Scalar lsquared    = (x1 - x2).squaredNorm();

    const Scalar Egl = 0.5 * (lsquared - LRefsquared) / LRefsquared;

    return 0.5 * EA / sqrt(LRefsquared) * Egl * Egl;
  }

private:
  LocalView localView_;
  double EA;
};

int main() {
  /// Construct grid
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
  /// Construct basis
  auto basis = makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));

  /// Create finite elements
  const double EA = 100;
  std::vector<Truss<decltype(basis)>> fes;
  for (auto& ele : elements(gridView))
    fes.emplace_back(basis, ele, EA);

  /// Collect dirichlet nodes
  std::vector<bool> dirichletFlags(basis.size(), false);
  Dune::Functions::forEachBoundaryDOF(basis, [&](auto&& index) { dirichletFlags[index] = true; });

  /// Create assembler
  auto denseFlatAssembler = DenseFlatSimpleAssembler(basis, fes, dirichletFlags);

  /// Create non-linear operator
  double lambda = 0;
  Eigen::VectorXd d;
  d.setZero(basis.size());

  auto RFunction = [&](auto&& u, auto&& lambdalocal) -> auto& {
    auto& R = denseFlatAssembler.getVector(forces, u, lambdalocal);
    R[3] -= -lambdalocal;
    return R;
  };
  auto KFunction = [&](auto&& u, auto&& lambdalocal) -> auto& {
    return denseFlatAssembler.getMatrix(stiffness, u, lambdalocal);
  };

  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(RFunction, KFunction), parameter(d, lambda));

  /// Choose linear solver
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT);

  /// Create Nonlinear solver for controlroutine, i.e. a Newton-Rahpson object
  auto nr = Ikarus::NewtonRaphson(nonLinOp, std::move(linSolver));
  nr.setup({.tol = 1e-8, .maxIter = 100});

  /// Create Observer to write information of the non-linear solver on the console
  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();

  const int loadSteps = 10;
  Eigen::Matrix3Xd lambdaAndDisp;
  lambdaAndDisp.setZero(Eigen::NoChange, loadSteps + 1);
  /// Create Observer which executes when control routines messages SOLUTION_CHANGED
  auto lvkObserver = std::make_shared<Ikarus::GenericControlObserver>(ControlMessages::SOLUTION_CHANGED, [&](int step) {
    lambdaAndDisp(0, step) = lambda;
    lambdaAndDisp(1, step) = d[2];
    lambdaAndDisp(2, step) = d[3];
  });

  /// Create Observer which writes vtk files when control routines messages SOLUTION_CHANGED
  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<decltype(basis)>>(basis, d, 2);
  vtkWriter->setFieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2);
  vtkWriter->setFileNamePrefix("TestTruss");
  nr.subscribeAll(nonLinearSolverObserver);

  /// Create loadcontrol
  auto lc = Ikarus::LoadControl(std::move(nr), loadSteps, {0, 30});
  lc.subscribeAll({vtkWriter, lvkObserver});

  /// Execute!
  lc.run();

  /// Postprocess
  using namespace matplot;
  Eigen::VectorXd lambdaVec = lambdaAndDisp.row(0);
  Eigen::VectorXd dVec      = -lambdaAndDisp.row(2);
  auto f                    = figure(true);

  title("Load-Displacement Curve");
  xlabel("y-Displacement");
  ylabel("LoadFactor");

  auto analyticalLoadDisplacementCurve = [&](auto& w) {
    const double Ltruss = std::sqrt(h * h + L * L);
    return EA * Dune::power(h, 3) / Dune::power(Ltruss, 3)
           * (w / h - 1.5 * Dune::power(w / h, 2) + 0.5 * Dune::power(w / h, 3));
  };

  std::vector<double> x  = linspace(0.0, dVec.maxCoeff());
  std::vector<double> y1 = transform(x, [&](auto x) { return analyticalLoadDisplacementCurve(x); });
  auto p                 = plot(x, y1, dVec, lambdaVec);
  p[0]->line_width(2);
  p[1]->line_width(2);
  p[1]->marker(line_spec::marker_style::asterisk);
  show();

}