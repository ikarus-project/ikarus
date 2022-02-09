//
// Created by Alex on 21.07.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"

#include <dune/alugrid/grid.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/iga/nurbsgrid.hh>

#include "spdlog/spdlog.h"

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/Controlroutines/LoadControl.h"
#include "ikarus/FiniteElements/NonLinearElasticityFEwithBasis.h"
#include "ikarus/Solver/NonLinearSolver/NewtonRaphson.hpp"
#include "ikarus/basis/basishelper.h"
#include "ikarus/utils/Observer/controlVTKWriter.h"
#include "ikarus/utils/Observer/nonLinearSolverLogger.h"
#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/LinearAlgebra/NonLinearOperator.h>
#include <ikarus/utils/concepts.h>
#include <ikarus/utils/utils/algorithms.h>

template <Ikarus::Concepts::FlatIndexBasis Basis, typename FEContainer>
class DenseFlatAssembler {
public:
  using RequirementType = typename FEContainer::value_type::FERequirementType;
  explicit DenseFlatAssembler(const Basis& basis, const FEContainer& fes, const std::vector<bool>& dirichFlags)
      : basis_{&basis}, feContainer{fes}, dirichletFlags{&dirichFlags} {}

  Eigen::MatrixXd& getMatrix(const Ikarus::MatrixAffordances& p_matrixAffordances,const Eigen::VectorXd& displacement, const double& lambda) {
    return getMatrixImpl(p_matrixAffordances,displacement, lambda);
  }

  Eigen::VectorXd& getVector(const Ikarus::VectorAffordances& p_vectorAffordances,const Eigen::VectorXd& displacement, const double& lambda) {
    return getVectorImpl(p_vectorAffordances,displacement, lambda);
  }

  double getScalar(const Ikarus::ScalarAffordances & p_scalarAffordances,const Eigen::VectorXd& displacement, const double& lambda) {
    return getScalarImpl(p_scalarAffordances,displacement, lambda);
  }

private:
  Eigen::MatrixXd& getMatrixImpl(const Ikarus::MatrixAffordances& p_matrixAffordances,const Eigen::VectorXd& displacement, const double& lambda) {
    mat.setZero(basis_->size(), basis_->size());
    RequirementType requirements;
    requirements.matrixAffordances = p_matrixAffordances;
    requirements.sols.emplace_back(displacement);
    requirements.parameter.insert({Ikarus::FEParameter::loadfactor,lambda});
    for (auto& fe : feContainer) {
      auto matLoc      = fe.calculateMatrix(requirements);
      auto globalIndices = fe.globalIndices();
      for (auto i = 0; auto idi : fe.globalIndices()) {
        for (auto j = 0; auto idj : fe.globalIndices()) {
          mat(idi[0], idj[0]) += matLoc(i, j);
          ++j;
        }
        ++i;
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

  Eigen::VectorXd& getVectorImpl(const Ikarus::VectorAffordances& p_vectorAffordances, const Eigen::VectorXd& displacement, const double& lambda) {
    vec.setZero(basis_->size());
    auto localView = basis_->localView();
    RequirementType requirements;
    requirements.vectorAffordances = p_vectorAffordances;
    requirements.sols.emplace_back(displacement);
    requirements.parameter.insert({Ikarus::FEParameter::loadfactor,lambda});
    for (auto& fe : feContainer) {
//      Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis<decltype(localView)> fe(localView, 1000, 0.3);

      auto vecLocal = fe.calculateVector(requirements);
      for (int i = 0;auto id : fe.globalIndices()) {
        vec(id[0]) += vecLocal(i);
        ++i;
      }
    }
    for (auto i = 0U; i < basis_->size(); ++i) {
      if (dirichletFlags->at(i)) vec[i] = 0;
    }

    return vec;
  }

  double getScalarImpl(const Ikarus::ScalarAffordances& p_scalarAffordances,const Eigen::VectorXd& displacement, const double& lambda) {
    double scalar = 0.0;
    vec.setZero(basis_->size());
    RequirementType requirements;
    requirements.scalarAffordances = p_scalarAffordances;
    requirements.sols.emplace_back(displacement);
    requirements.parameter.insert({Ikarus::FEParameter::loadfactor,lambda});
    auto localView = basis_->localView();

    for (auto& ge : elements(basis_->gridView())) {
      localView.bind(ge);
      Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis<decltype(localView)> fe(localView, 1000, 0.3);
      for (auto i = 0U; i < localView.size(); ++i)
        scalar += fe.calculateScalar(requirements);
    }

    return scalar;
  }
  Basis const* basis_;
  FEContainer const& feContainer;
  std::vector<bool> const* dirichletFlags;
  Eigen::MatrixXd mat{};
  Eigen::VectorXd vec{};
};

GTEST_TEST(LoadControlTestWithUGGrid, GridLoadControlTestWithUGGrid) {
  using namespace Ikarus;
  constexpr int gridDim = 2;
  /// ALUGrid Example
//        using Grid            = Dune::ALUGrid<gridDim, 2, Dune::simplex, Dune::conforming>;
//        auto grid             =
//        Dune::GmshReader<Grid>::read("../../tests/src/testFiles/unstructuredTrianglesfine.msh", false);

  /// IGA Grid Example
//  constexpr auto dimworld              = 2;
//  const std::array<int, gridDim> order = {2, 2};
//
//  const std::array<std::vector<double>, gridDim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};
//
//  using ControlPoint = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointType;
//
//  const std::vector<std::vector<ControlPoint>> controlPoints
//      = {{{.p = {0, 0}, .w = 5}, {.p = {0.5, 0}, .w = 1}, {.p = {1, 0}, .w = 1}},
//         {{.p = {0, 0.5}, .w = 1}, {.p = {0.5, 0.5}, .w = 10}, {.p = {1, 0.5}, .w = 1}},
//         {{.p = {0, 1}, .w = 1}, {.p = {0.5, 1}, .w = 1}, {.p = {1, 1}, .w = 1}}};
//
//  std::array<int, gridDim> dimsize
//      = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};
//
//  auto controlNet = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);
//  using Grid      = Dune::IGA::NURBSGrid<gridDim, dimworld>;
//
//  Dune::IGA::NURBSPatchData<gridDim, dimworld> patchData;
//  patchData.knotSpans     = knotSpans;
//  patchData.degree        = order;
//  patchData.controlPoints = controlNet;
//  auto grid               = std::make_shared<Grid>(patchData);
//  grid->globalRefine(5);

  /// YaspGrid Example
  using Grid        = Dune::YaspGrid<gridDim>;
  const double L    = 1;
  const double h    = 1;
  const size_t elex = 1;
  const size_t eley = 1;

  Dune::FieldVector<double, 2> bbox = {L, h};
  std::array<int, 2> eles           = {elex, eley};
  auto grid                         = std::make_shared<Grid>(bbox, eles);

  using GridView    = typename Grid::LeafGridView;
  GridView gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  //  auto basis = makeBasis(gridView, power<gridDim>(gridView.getPreBasis(), FlatInterleaved()));
  auto basis = makeBasis(gridView, power<gridDim>(lagrange<1>(), FlatInterleaved()));
  std::cout << "This gridview cotains: " << std::endl;
  std::cout << gridView.size(2) << " vertices" << std::endl;
  std::cout << gridView.size(1) << " edges" << std::endl;
  std::cout << gridView.size(0) << " elements" << std::endl;
  std::cout << basis.size() << " Dofs" << std::endl;

  draw(gridView);


  std::vector<Ikarus::FiniteElements::IFiniteElement<decltype(basis)::LocalView,Eigen::VectorXd>> fes;
  auto localView = basis.localView();
  for (auto& element : elements(basis.gridView())) {
    localView.bind(element);
//    decltype(localView)::MultiIndex
    fes.emplace_back(Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis<decltype(localView)>(localView, 1000, 0.3));
  }

  std::vector<bool> dirichletFlags(basis.size());

  Ikarus::markDirichletBoundaryDofs(basis, dirichletFlags,
                                    [](auto&& centerCoord) { return (std::abs(centerCoord[1]) < 1e-8); });
  Ikarus::utils::printContent(std::cout, dirichletFlags);
  auto denseAssembler = DenseFlatAssembler(basis, fes, dirichletFlags);

  Eigen::VectorXd d;
  d.setZero(basis.size());
  double lambda = 0.0;

  auto fintFunction   = [&](auto&& lambda, auto&& disp) -> auto& { return denseAssembler.getVector(forces,disp, lambda); };
  auto KFunction      = [&](auto&& lambda, auto&& disp) -> auto& { return denseAssembler.getMatrix(stiffness,disp, lambda); };
  auto energyFunction = [&](auto&& lambda, auto&& disp) -> auto { return denseAssembler.getScalar(potentialEnergy,disp, lambda); };
  std::cout<<KFunction(lambda,d)<<std::endl;
  auto nonLinOp  = Ikarus::NonLinearOperator(linearAlgebraFunctions(energyFunction, fintFunction, KFunction),
                                             parameter(lambda, d));
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT);

  auto nr                      = Ikarus::NewtonRaphson(nonLinOp.subOperator<1, 2>(), std::move(linSolver));
  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();

  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<decltype(basis)>>(basis, d, 2);
  vtkWriter->setFileNamePrefix("TestIGA");
  vtkWriter->setVertexSolutionName("displacement");
  nr.subscribeAll(nonLinearSolverObserver);

  auto lc = Ikarus::LoadControl(std::move(nr), 20, {0, 2000});

  lc.subscribeAll(vtkWriter);
  std::cout << "Energy before: " << nonLinOp.value() << std::endl;
  lc.run();
  nonLinOp.update<0>();
  std::cout << "Energy after: " << nonLinOp.value() << std::endl;
}