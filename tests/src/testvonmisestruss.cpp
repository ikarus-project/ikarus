// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <matplot/matplot.h>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/localfefunctions/eigenDuneTransformations.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <autodiff/forward/dual/dual.hpp>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlroutines/adaptivestepsizing.hh>
#include <ikarus/controlroutines/pathfollowingtechnique.hh>
#include <ikarus/finiteelements/feBases/autodifffe.hh>
#include <ikarus/finiteelements/feBases/powerbasisfe.hh>
#include <ikarus/linearalgebra/nonlinearoperator.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/solver/nonLinearsolver/newtonraphson.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/duneutilities.hh>
#include <ikarus/utils/eigendunetransformations.hh>
#include <ikarus/utils/fancyactivationfunctions.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controllogger.hh>
#include <ikarus/utils/observer/controlvtkwriter.hh>
#include <ikarus/utils/observer/genericobserver.hh>
#include <ikarus/utils/observer/nonlinearsolverlogger.hh>

using namespace Ikarus;
template <typename Basis_, typename FERequirements_ = FErequirements<>, bool useEigenRef = false>
class Truss : public PowerBasisFE<typename Basis_::FlatBasis> {
public:
  using Basis             = Basis_;
  using FlatBasis         = typename Basis::FlatBasis;
  using BaseDisp          = Ikarus::PowerBasisFE<FlatBasis>;
  using LocalView         = typename FlatBasis::LocalView;
  using Element           = typename LocalView::Element;
  using Geometry          = typename Element::Geometry;
  using FERequirementType = FERequirements_;
  using Traits            = TraitsFromLocalView<LocalView, useEigenRef>;
  Truss(const Basis &basis, const typename LocalView::Element &element, double p_EA)
      : BaseDisp(basis.flat(), element), EA{p_EA} {
    this->localView().bind(element);
  }

  inline double calculateScalar(const FERequirementType &par) const { return calculateScalarImpl<double>(par); }

protected:
  template <typename ScalarType>
  auto calculateScalarImpl(const FERequirementType &par, const std::optional<const Eigen::VectorX<ScalarType>> &dx
                                                         = std::nullopt) const -> ScalarType {
    const auto &d      = par.getGlobalSolution(Ikarus::FESolutions::displacement);
    const auto &lambda = par.getParameter(FEParameter::loadfactor);

    auto &ele     = this->localView().element();
    const auto X1 = Dune::toEigen(ele.geometry().corner(0));
    const auto X2 = Dune::toEigen(ele.geometry().corner(1));

    Eigen::Matrix<ScalarType, Traits::worlddim, 2> u;
    u.setZero();
    if (dx) {
      for (int i = 0; i < 2; ++i)
        for (int k2 = 0; k2 < Traits::worlddim; ++k2)
          u.col(i)(k2) = dx.value()[Traits::worlddim * i + k2]
                         + d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
    } else {
      for (int i = 0; i < 2; ++i)
        for (int k2 = 0; k2 < Traits::worlddim; ++k2)
          u.col(i)(k2) = d[this->localView().index(this->localView().tree().child(k2).localIndex(i))[0]];
    }

    const Eigen::Vector2<ScalarType> x1 = X1 + u.col(0);
    const Eigen::Vector2<ScalarType> x2 = X2 + u.col(1);

    const double LRefsquared  = (X1 - X2).squaredNorm();
    const ScalarType lsquared = (x1 - x2).squaredNorm();

    const ScalarType Egl = 0.5 * (lsquared - LRefsquared) / LRefsquared;

    return 0.5 * EA * Egl * Egl * sqrt(LRefsquared);
  }

private:
  double EA;
};

int main(int argc, char **argv) {
  Ikarus::init(argc, argv);
  /// Construct grid
  Dune::GridFactory<Dune::FoamGrid<1, 2, double>> gridFactory;
  const double h = 25.0;
  const double L = 10.0;
  gridFactory.insertVertex({0, 0});
  gridFactory.insertVertex({L, h});
  gridFactory.insertVertex({2 * L, 0});
  gridFactory.insertElement(Dune::GeometryTypes::line, {0, 1});
  gridFactory.insertElement(Dune::GeometryTypes::line, {1, 2});
  auto grid     = gridFactory.createGrid();
  auto gridView = grid->leafGridView();

  /// Construct basis
  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>()));

  /// Create finite elements
  const double EA = 3000;
  std::vector<AutoDiffFE<Truss<decltype(basis)>>> fes;
  for (auto &ele : elements(gridView))
    fes.emplace_back(basis, ele, EA);

  /// Collect dirichlet nodes
  auto basisP = std::make_shared<const decltype(basis)>(basis);
  Ikarus::DirichletValues dirichletValues(basisP->flat());
  dirichletValues.fixBoundaryDOFs(
      [&](auto &dirichletFlags, auto &&globalIndex) { dirichletFlags[globalIndex] = true; });

  /// Create assembler
  auto denseFlatAssembler = DenseFlatAssembler(fes, dirichletValues);

  /// Create non-linear operator
  double lambda = 0;
  Eigen::VectorXd d;
  d.setZero(basis.flat().size());

  auto req = FErequirements().addAffordance(Ikarus::AffordanceCollections::elastoStatics);

  auto RFunction = [&](auto &&u, auto &&lambdaLocal) -> auto & {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, u)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    auto &R = denseFlatAssembler.getVector(req);
    R[3] -= -lambdaLocal;
    return R;
  };
  auto KFunction = [&](auto &&u, auto &&lambdaLocal) -> auto & {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, u)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return denseFlatAssembler.getMatrix(req);
  };

  auto nonLinOp = Ikarus::NonLinearOperator(functions(RFunction, KFunction), parameter(d, lambda));

  /// Choose linear solver
  auto linSolver = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);

  /// Create Nonlinear solver for control routine, i.e. a Newton-Raphson object
  auto nr = Ikarus::makeNewtonRaphsonWithSubsidiaryFunction(nonLinOp, std::move(linSolver));
  nr->setup({.tol = 1e-6, .maxIter = 20});

  /// Create Observer to write information of the non-linear solver on the
  /// console
  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
  auto controlObserver         = std::make_shared<ControlLogger>();
  auto pft                     = Ikarus::StandardArcLength{};
  auto afArgs
      = Ikarus::AdaptiveStepSizing::ActivationFunctionArgs{.r_0 = 0.01, .r_i = 0.05, .r_1 = 0.2, .lambda_i = 0.1};
  auto af         = Ikarus::AdaptiveStepSizing::LinearPiecewiseFunction(afArgs);
  auto adaptiveSS = Ikarus::FancyAdaptiveStepSizing<>();
  int loadSteps   = 13;
  double stepSize = 0.1;
  auto alc        = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft, adaptiveSS, af);
  Eigen::Matrix3Xd lambdaAndDisp;
  lambdaAndDisp.setZero(Eigen::NoChange, loadSteps);
  /// Create Observer which executes when control routines messages
  /// SOLUTION_CHANGED
  auto lvkObserver
      = std::make_shared<Ikarus::GenericObserver<ControlMessages>>(ControlMessages::SOLUTION_CHANGED, [&](int step) {
          lambdaAndDisp(0, step) = lambda;
          lambdaAndDisp(1, step) = d[2];
          lambdaAndDisp(2, step) = d[3];
        });

  /// Create Observer which writes vtk files when control routines messages
  /// SOLUTION_CHANGED
  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>>>(
      basis.flat(), d, 2);
  vtkWriter->setFieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2);
  vtkWriter->setFileNamePrefix("vonMisesTruss");
  nr->subscribeAll(nonLinearSolverObserver);
  alc.subscribeAll({controlObserver, vtkWriter, lvkObserver});

  /// Execute!
  const auto controlInfo = alc.run();

  /// Postprocess
  using namespace matplot;
  Eigen::VectorXd lambdaVec = lambdaAndDisp.row(0);
  Eigen::VectorXd dVec      = -lambdaAndDisp.row(2);  // vertical displacement
  auto f                    = figure(true);
  auto ax                   = gca();
  title("Load-Displacement Curve");
  xlabel("y-Displacement");
  ylabel("LoadFactor");

  auto analyticalLoadDisplacementCurve = [&](auto &w) {
    const double Ltruss = std::sqrt(h * h + L * L);
    return (2.0 * EA * Dune::power(h, 3) / Dune::power(Ltruss, 3))
           * (w / h - 1.5 * Dune::power(w / h, 2) + 0.5 * Dune::power(w / h, 3));
  };

  std::cout << std::setprecision(16) << "Max Lambda:\t" << lambdaVec.maxCoeff() << std::endl;
  std::cout << std::setprecision(16) << "Max Disp  :\t" << dVec.maxCoeff() << std::endl;
  std::vector<std::string> legends;
  legends.push_back("Analytical");
  legends.push_back("FE");

  std::vector<double> x  = linspace(0.0, dVec.maxCoeff());
  std::vector<double> y1 = transform(x, [&](auto x) { return analyticalLoadDisplacementCurve(x); });
  auto p                 = plot(x, y1, dVec, lambdaVec);
  p[0]->line_width(3);
  p[1]->line_width(2);
  p[1]->marker(line_spec::marker_style::asterisk);
  ax->legend(legends);
  auto legend = ax->legend();
  legend->location(legend::general_alignment::topleft);
  save("vonMisesTruss.png");
  // f->draw();
  //  using namespace std::chrono_literals;
  //  std::this_thread::sleep_for(5s);
}