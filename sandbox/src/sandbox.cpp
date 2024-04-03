// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <iostream>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/localfefunctions/eigenDuneTransformations.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/finiteelements/autodifffe.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/eigendunetransformations.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/nonlinearoperator.hh>
#include <ikarus/utils/observer/controlvtkwriter.hh>
#include <ikarus/utils/observer/nonlinearsolverlogger.hh>
using namespace Ikarus;
template <typename PreFE, typename FE>
class Truss;

struct TrussPre
{
  double EA;

  template <typename PreFE, typename FE>
  using Skill = Truss<PreFE, FE>;
};

template <typename PreFE, typename FE>
class Truss
{
public:
  using Traits            = typename PreFE::Traits;
  using BasisHandler      = typename Traits::BasisHandler;
  using FlatBasis         = typename Traits::FlatBasis;
  using FERequirementType = typename Traits::FERequirementType;
  using LocalView         = typename Traits::LocalView;
  using Geometry          = typename Traits::Geometry;
  using Element           = typename Traits::Element;
  using Pre               = TrussPre;

  Truss(Pre pre)
      : EA{pre.EA} {}

protected:
  template <template <typename, int, int> class RT>
  requires Dune::AlwaysFalse<RT<double, 1, 1>>::value
  auto calculateAtImpl(const FERequirementType& req, const Dune::FieldVector<double, Traits::mydim>& local,
                       Dune::PriorityTag<0>) const {}

  template <typename ScalarType>
  auto calculateScalarImpl(const FERequirementType& par,
                           const std::optional<std::reference_wrapper<const Eigen::VectorX<ScalarType>>>& dx =
                               std::nullopt) const -> ScalarType {
    const auto& d         = par.getGlobalSolution(Ikarus::FESolutions::displacement);
    const auto& lambda    = par.getParameter(FEParameter::loadfactor);
    const auto& localView = underlying().localView();
    const auto& tree      = localView.tree();
    auto& ele             = localView.element();
    const auto X1         = Dune::toEigen(ele.geometry().corner(0));
    const auto X2         = Dune::toEigen(ele.geometry().corner(1));

    Eigen::Matrix<ScalarType, Traits::worlddim, 2> u;
    u.setZero();
    if (dx) {
      for (int i = 0; i < 2; ++i)
        for (int k2 = 0; k2 < Traits::worlddim; ++k2)
          u.col(i)(k2) =
              dx.value().get()[Traits::worlddim * i + k2] + d[localView.index(tree.child(k2).localIndex(i))[0]];
    } else {
      for (int i = 0; i < 2; ++i)
        for (int k2 = 0; k2 < Traits::worlddim; ++k2)
          u.col(i)(k2) = d[localView.index(tree.child(k2).localIndex(i))[0]];
    }

    const Eigen::Vector2<ScalarType> x1 = X1 + u.col(0);
    const Eigen::Vector2<ScalarType> x2 = X2 + u.col(1);

    const double LRefsquared  = (X1 - X2).squaredNorm();
    const ScalarType lsquared = (x1 - x2).squaredNorm();

    const ScalarType Egl = 0.5 * (lsquared - LRefsquared) / LRefsquared;

    return 0.5 * EA * sqrt(LRefsquared) * Egl * Egl;
  }

private:
  //> CRTP
  const auto& underlying() const { return static_cast<const FE&>(*this); }
  auto& underlying() { return static_cast<FE&>(*this); }
  double EA;
};

auto truss(double EA) { return TrussPre(EA); }

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  constexpr int worldDimension= 2;
  Dune::GridFactory<Dune::FoamGrid<1, worldDimension, double>> gridFactory;

  const double h = 1.0;
  const double L = 2.0;
  gridFactory.insertVertex({0, 0});
  gridFactory.insertVertex({L, h});
  gridFactory.insertVertex({2 * L, 0});
  gridFactory.insertElement(Dune::GeometryTypes::line, {0, 1});
  gridFactory.insertElement(Dune::GeometryTypes::line, {1, 2});
  auto grid     = gridFactory.createGrid();
  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>()));


  const double EA  = 100;
  auto sk          = Ikarus::skills(truss(EA));
  using AutoDiffFE = Ikarus::AutoDiffFE<decltype(makeFE(basis, sk))>;
  std::vector<AutoDiffFE> fes;
  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(AutoDiffFE(makeFE(basis, sk)));
    fes.back().bind(ge);
  }

  Ikarus::DirichletValues dirichletValues(basis.flat());
dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& globalIndex) { dirichletFlags[globalIndex] = true; });

  auto assembler = Ikarus::SparseFlatAssembler(fes, dirichletValues);

  double lambda = 0;
  Eigen::VectorXd d;
  d.setZero(basis.flat().size());

  auto req = Ikarus::FERequirements().addAffordance(Ikarus::AffordanceCollections::elastoStatics);

  auto RFunction = [&](auto&& u, auto&& lambdaLocal) -> auto {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, u)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    auto R = assembler.getVector(req);
    R[3] -= -lambdaLocal;
    return R;
  };

  auto KFunction = [&](auto&& u, auto&& lambdaLocal) -> auto& {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, u)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return assembler.getMatrix(req);
  };


  auto nonLinOp = Ikarus::NonLinearOperator(functions(RFunction, KFunction), parameter(d, lambda));

  auto linSolver = Ikarus::LinearSolver(Ikarus::SolverTypeTag::sd_SimplicialLDLT);

  /// Create Nonlinear solver for controlroutine, i.e. a Newton-Rahpson object
  auto nr = Ikarus::makeNewtonRaphson(nonLinOp, std::move(linSolver));
  nr->setup({.tol = 1e-8, .maxIter = 100});

  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
  nr->subscribeAll(nonLinearSolverObserver);

  // const int loadSteps = 10;
  // auto lc = Ikarus::LoadControl(std::move(nr), loadSteps, {0, 30});

  lambda+=7;

  nonLinOp.derivative();
  auto& K = KFunction(d,lambda); //nonLinOp.value();
  auto R = RFunction(d,lambda); //nonLinOp.derivative();
  Eigen::SparseLU<std::remove_cvref_t<decltype(K)>> ld;
  ld.compute(K);
  if (ld.info() != Eigen::Success)
    DUNE_THROW(Dune::MathError, "Failed Compute");

  d -= ld.solve(R);
  if (ld.info() != Eigen::Success)
    DUNE_THROW(Dune::MathError, "Failed Solve");

  auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, worldDimension>>(basis.flat(), d);


  Dune::VTKWriter vtkWriter(gridView, Dune::VTK::nonconforming);
  vtkWriter.addVertexData(disp, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, worldDimension));
  vtkWriter.write("iks007_vonMisesTruss_GeoLin");
  // auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>>>(
  // basis.flat(), d, 1);
  // vtkWriter->setFieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, worldDimension);
  // vtkWriter->setFileNamePrefix("iks007_vonMisesTruss");
  // lc.subscribeAll(vtkWriter);
  //
  // lc.run();

}
