//
// Created by Alex on 21.07.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"

#include <dune/alugrid/grid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/iga/nurbsgrid.hh>

#include "spdlog/spdlog.h"

#include <Eigen/Core>
#include <Eigen/Dense>

//#include "ikarus/Assembler/SimpleAssemblers.h"
//#include "ikarus/Controlroutines/LoadControl.h"
//#include "ikarus/FEManager/DefaultFEManager.h"
//#include "ikarus/FiniteElements/ElasticityFE.h"
//#include "ikarus/FiniteElements/NonLinearElasticityFE.h"
//#include "ikarus/LinearAlgebra/DirichletConditionManager.h"
//#include "ikarus/utils/Observer/controlLogger.h"
//#include "ikarus/utils/Observer/gridDrawerObserver.h"
//#include "ikarus/utils/Observer/nonLinearSolverLogger.h"
//#include <ikarus/FiniteElements/ForceLoad.h>
#include <ikarus/Grids/GridHelper/griddrawer.h>
//#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>
#include <ikarus/LinearAlgebra/NonLinearOperator.h>

// TEST(LoadControlTestWithYaspGrid, GridLoadControlTestWithYaspGrid) {
//   using namespace Ikarus::Grid;
//   using Grid = Dune::YaspGrid<2>;
//   //   gridFactory;
//   using vertexType = Eigen::Vector2d;
//   std::vector<vertexType> verticesVec;
//   const double L    = 10;
//   const double h    = 1;
//   const size_t elex = 2;
//   const size_t eley = 2;
//
//   Dune::FieldVector<double, 2> bbox = {L, h};
//   std::array<int, 2> eles           = {elex, eley};
//   Grid grid(bbox, eles);
//   auto gridView = grid.leafGridView();
//
//   draw(gridView);
//   std::cout << "1" << std::endl;
//   std::vector<Ikarus::FiniteElements::IFiniteElement> feContainer;
//   std::cout << "2" << std::endl;
//   //  for (auto&& ge : rootEntities(gridView))
//   //    feContainer.emplace_back(Ikarus::FiniteElements::ElasticityFE(ge, gridView.indexSet(), 1000, 0.0));
//   std::cout << "==============================" << std::endl;
//   const auto& indexSet = gridView.indexSet();
//   for (auto&& ge : elements(gridView)) {
//     for (unsigned int id = 0; id < ge.subEntities(2); ++id) {
//       auto vertex = ge.template subEntity<2>(id);
//       std::cout << "id: " << id << " "
//                 << "Coords: " << vertex.geometry().corner(0) << " globID: " << indexSet.index(vertex) << " "
//                 << indexSet.subIndex(ge, id, 2) << std::endl;
//     }
//     feContainer.emplace_back(Ikarus::FiniteElements::NonLinearElasticityFE(ge, gridView.indexSet(), 1000, 0.3));
//   }
//   std::cout << "==============================" << std::endl;
//   std::cout << "3" << std::endl;
//   auto spaceFunction = [](const Eigen::Vector2d&) -> Eigen::Vector2d {
//     Eigen::Vector2d f{};
//     f[1] = 1;
//     return f;
//   };
//   std::cout << "4" << std::endl;
//   for (auto& edge : edges(gridView)) {
//     if (std::abs(edge.geometry().center()[1]) > h - 0.01) {
//       std::cout << "edge.geometry().center() " << edge.geometry().center() << std::endl;
//       feContainer.emplace_back(Ikarus::FiniteElements::ForceLoad(edge, gridView.indexSet(), spaceFunction));
//     }
//   }
//   std::cout << "5" << std::endl;
//   auto feManager = Ikarus::FEManager::DefaultFEManager(feContainer, gridView);
//   std::cout << "6" << std::endl;
//   Ikarus::DirichletConditionManager dirichletConditionManager(feManager);
//   for (auto& vertex : vertices(gridView)) {
//     if (std::abs(vertex.geometry().corner(0)[1]) < 1e-8) {
//       //      std::cout << "AddConstraint at " << vertex.geometry().corner(0) << std::endl;
//       dirichletConditionManager.addConstraint(vertex, 0);
//       dirichletConditionManager.addConstraint(vertex, 1);
//     }
//   }
//   std::cout << "====asdasdasdasd=====" << std::endl;
//   dirichletConditionManager.finalize();
//   std::cout << "7" << std::endl;
//   auto vectorAssembler = Ikarus::Assembler::VectorAssembler(feManager, dirichletConditionManager);
//
//   auto denseMatrixAssembler  = Ikarus::Assembler::DenseMatrixAssembler(feManager, dirichletConditionManager);
//   auto sparseMatrixAssembler = Ikarus::Assembler::SparseMatrixAssembler(feManager, dirichletConditionManager);
//   std::cout << "8" << std::endl;
//   [[maybe_unused]] auto fintFunction = [&](auto&& lambda) -> auto& {
//     return vectorAssembler.getReducedVector(Ikarus::FiniteElements::forces, lambda);
//   };
//   [[maybe_unused]] auto KFunction = [&](auto&& lambda) -> auto& {
//     return denseMatrixAssembler.getReducedMatrix(Ikarus::FiniteElements::stiffness, lambda);
//   };
//   [[maybe_unused]] auto KFunctionSparse = [&](auto&& lambda) -> auto& {
//     return sparseMatrixAssembler.getReducedMatrix(Ikarus::FiniteElements::stiffness, lambda);
//   };
//
//   auto controlObserver         = std::make_shared<ControlLogger>();
//   auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
//   auto gridDrawerObserver
//       = std::make_shared<GridDrawerObserver<decltype(gridView), decltype(feManager)>>(gridView, feManager);
//   auto time          = Ikarus::FEParameterFactory::createParameter(Ikarus::FEParameter::time, 1);
//   time.value[0]      = 0.0;
//   auto& Fint         = fintFunction(time);
//   const int dofsFree = ((elex + 1) * (eley + 1) - (elex + 1)) * 2;
//   EXPECT_EQ(Fint.size(), dofsFree);
//   auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::SparseLU);
//   auto nonLinOp  = Ikarus::NonLinearOperator(linearAlgebraFunctions(fintFunction, KFunctionSparse), parameter(time));
//   auto nr        = Ikarus::NewtonRaphson(
//              nonLinOp, std::move(linSolver),
//              [&dirichletConditionManager](decltype(feManager.getVariables())& x, const Eigen::VectorX<double>& D) {
//         x += dirichletConditionManager.viewAsFullVector(D);
//              });
//   nr.subscribeAll(nonLinearSolverObserver);
//
//   auto lc = Ikarus::LoadControl(feManager, dirichletConditionManager, std::move(nr), 10, {0, 5000});
//   //  auto lc = makeLoadControl<Ikarus::NewtonRaphson>(
//   //      feManager, dirichletConditionManager, linearAlgebraFunctions(fintFunction, KFunctionSparse),
//   //      Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::SparseLU), 10, {0, 1});
//   lc.subscribeAll(controlObserver);
//   //  lc.subscribeToNonLinearSolver(nonLinearSolverObserver);
//   lc.subscribe(ControlMessages::SOLUTION_CHANGED, gridDrawerObserver);
//   lc.run();
// }

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include "ikarus/FiniteElements/NonLinearElasticityFEwithBasis.h"
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/utils/concepts.h>

template <Ikarus::Concepts::FlatIndexBasis Basis, typename FEContainer>
class DenseFlatAssembler {
public:
  explicit DenseFlatAssembler(const Basis& basis, const FEContainer& fes, const std::vector<bool>& dirichFlags)
      : basis_{&basis}, feContainer{&fes}, dirichletFlags{&dirichFlags} {}

  Eigen::MatrixXd& getMatrix(const Eigen::VectorXd& displacement, const double& lambda) {
    return getMatrixImpl(displacement, lambda);
  }

  Eigen::VectorXd& getVector(const Eigen::VectorXd& displacement, const double& lambda) {
    return getVectorImpl(displacement, lambda);
  }

private:
  Eigen::MatrixXd& getMatrixImpl(const Eigen::VectorXd& displacement, const double& lambda) {
    mat.setZero(basis_->size(), basis_->size());
    auto localView = basis_->localView();
    for (auto& ge : elements(basis_->gridView())) {
      localView.bind(ge);
      Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis fe(localView, 1000, 0.0);
      auto matLoc      = fe.calculateMatrix(displacement, lambda);
      auto first_child = localView.tree().child(0);
      for (auto i = 0U; i < localView.size(); ++i)
        for (auto j = 0U; j < localView.size(); ++j) {
          mat(localView.index(i)[0], localView.index(j)[0]) += matLoc(i, j);
        }
    }
    localView.unbind();
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
      Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis<decltype(localView)> fe(localView, 1000, 0.0);
      auto vecLocal = fe.calculateVector(displacement, lambda);
      for (auto i = 0U; i < localView.size(); ++i)
        vec(localView.index(i)[0]) += vecLocal(i);
    }
    for (auto i = 0U; i < basis_->size(); ++i) {
      if (dirichletFlags->at(i)) vec[i] = 0;
    }

    return vec;
  }
  Basis const* basis_;
  FEContainer const* feContainer;
  std::vector<bool> const* dirichletFlags;
  Eigen::MatrixXd mat{};
  Eigen::VectorXd vec{};
};

GTEST_TEST(LoadControlTestWithUGGrid, GridLoadControlTestWithUGGrid) {
  constexpr int gridDim = 2;
  //    using Grid            = Dune::ALUGrid<gridDim, 2, Dune::simplex, Dune::conforming>;
  //    auto grid             = Dune::GmshReader<Grid>::read("../../tests/src/testFiles/unstructuredTrianglesfine.msh",
  //    false);
  //
  //
  constexpr auto dimworld              = 2;
  const std::array<int, gridDim> order = {2, 2};

  const std::array<std::vector<double>, gridDim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};
  //  const std::vector<std::vector<FieldVector<double, dimworld> > > controlPointsold
  //      = {{{0, 0, 1}, {1, 0, 1}, {2, 0, 2}}, {{0, 1, 0}, {1, 1, 0}, {2, 1, 0}}, {{0, 2, 1}, {1, 2, 2}, {2, 2, 2}}};
  using ControlPoint = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0}, .w = 2}, {.p = {1, 0}, .w = 2}, {.p = {2, 0}, .w = 1}},
         {{.p = {0, 1}, .w = 1}, {.p = {1, 1}, .w = 4}, {.p = {2, 1}, .w = 1}},
         {{.p = {0, 2}, .w = 1}, {.p = {1, 2}, .w = 2}, {.p = {2, 2}, .w = 4}}};

  std::array<int, gridDim> dimsize
      = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<gridDim, dimworld>;

  Dune::IGA::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  auto grid               = std::make_shared<Grid>(patchData);
  grid->globalRefine(5);

  //  using Grid = Dune::YaspGrid<gridDim>;
  //  const double L    = 1;
  //  const double h    = 1;
  //  const size_t elex = 20;
  //  const size_t eley = 20;
  //
  //  Dune::FieldVector<double, 2> bbox = {L, h};
  //  std::array<int, 2> eles           = {elex, eley};
  //  auto grid                         = std::make_shared<Grid>(bbox, eles);

  //  auto nurbs = grid.getPreBasis();

  using GridView    = typename Grid::LeafGridView;
  GridView gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, power<gridDim>(gridView.getPreBasis(), FlatInterleaved()));
  //  auto basis = makeBasis(gridView, power<gridDim>(lagrange<1>(), FlatInterleaved()));
  std::cout << "This gridview cotains: " << std::endl;
  std::cout << gridView.size(2) << " vertices" << std::endl;
  std::cout << gridView.size(1) << " edges" << std::endl;
  std::cout << gridView.size(0) << " elements" << std::endl;

  draw(gridView);
  std::vector<Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis<decltype(basis)::LocalView>> fes;
  auto localView = basis.localView();
  for (auto& ge : elements(basis.gridView())) {
    localView.bind(ge);
    fes.emplace_back(localView, 1000, 0.0);
  }

  std::vector<bool> dirichletFlags(basis.size());

  Dune::Functions::forEachBoundaryDOF(basis, [&](auto&& localIndex, auto&& localView) {
    const Dune::FieldVector<double, 2> localCoordinate = {0.0, 0.5};
    const auto inElementLocalCoords = localView.element().geometry().global(localCoordinate);
//    std::cout << "edgeCenter: " << inElementLocalCoords << std::endl;
    if (std::abs(inElementLocalCoords[1]) > 1e-8)  // Only Fix Lower Edge
      return;
    auto& fe = localView.tree().child(0).finiteElement();
    std::vector<Dune::FieldVector<double, 1>> N;
    fe.localBasis().evaluateFunction(localCoordinate, N);
    for (size_t j = 0; j < N.size(); ++j) {
//      std::cout << N[j][0] << " ";
      if (N[j][0] > 1e-8) {
        for (size_t k = 0; k < localView.tree().CHILDREN; ++k) {
          auto localFixedDof                                = localView.tree().child(k).localIndex(j);
          dirichletFlags[localView.index(localFixedDof)[0]] = true;
        }
      }
    }
//    std::cout << std::endl;
  });
//  for (auto&& flag : dirichletFlags)
//    std::cout << flag << std::endl;

  auto denseAssembler = DenseFlatAssembler(basis, fes, dirichletFlags);

  Eigen::VectorXd d(basis.size());
  d.setZero();
  double lambda = 0.0;

  Dune::SubsamplingVTKWriter<GridView> vtkWriter(gridView, Dune::refinementLevels(1));
  //  Dune::VTKWriter<GridView> vtkWriter(gridView);
  vtkWriter.addVertexData(d, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, gridDim));
  vtkWriter.write("TestDuneBasisUSIGA_0");

  double fac = 100;
  for (int ls = 1; ls < 20; ++ls) {
    Eigen::FullPivLU<Eigen::MatrixXd> lu;
    for (int i = 0; i < 20; ++i) {
      const auto& K = denseAssembler.getMatrix(d, lambda);
      lu.compute(K);
      const auto& r = denseAssembler.getVector(d, lambda);
      const auto dd = lu.solve(r);
      d -= lu.solve(r);
      std::cout << "Rnorm: " << r.norm() << "    "
                << "dnorm: " << dd.norm() << "    "  << std::endl;
      if (r.norm() < 1e-8) break;
    }
    std::cout << "==============" << std::endl;
    auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis, d);
    Dune::SubsamplingVTKWriter<GridView> vtkWriterI(gridView, Dune::refinementLevels(1));
    //    Dune::VTKWriter<GridView> vtkWriterI(gridView);
    vtkWriterI.addVertexData(disp, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, gridDim));
    vtkWriterI.write("TestDuneBasisUSIGA_" + std::to_string(ls));
    lambda += fac;
  }
}