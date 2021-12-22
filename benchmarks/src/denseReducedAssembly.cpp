#include <benchmark/benchmark.h>
//#define EIGEN_SPARSEMATRIX_PLUGIN "libAddons/eigen/eigenSparseAddon.h"
#include <fstream>
#include <vector>
#include <Eigen/Core>

#include "ikarus/LinearAlgebra/DirichletConditionManager.h"
#include <ikarus/Assembler/SimpleAssemblers.h>
#include <ikarus/FEManager/DefaultFEManager.h>
#include <ikarus/FiniteElements/ElasticityFE.h>
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Geometries/GeometryType.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>

static void byhand(benchmark::State& state) {
  using namespace Ikarus::Grid;
  using Grid = SimpleGrid<2, 2>;
  SimpleGridFactory<2, 2> gridFactory;
  using vertexType = Eigen::Vector2d;
  std::vector<vertexType> verticesVec;
  const double L    = 10;
  const double h    = 1;
  const size_t elex = state.range(0);
  const size_t eley = state.range(0);
  for (size_t j = 0; j < eley + 1; ++j) {
    for (size_t i = 0; i < elex + 1; ++i)
      verticesVec.emplace_back(vertexType{i * L / (elex), j * h / (eley)});
  }
  for (auto&& vert : verticesVec)
    gridFactory.insertVertex(vert);
  for (size_t i = 0; i < elex; ++i) {
    for (size_t j = 0; j < eley; ++j) {
      std::vector<size_t> elementIndices;
      elementIndices.resize(4);
      elementIndices
          = {i + j * (elex + 1), i + j * (elex + 1) + 1, i + (j + 1) * (elex + 1), i + (j + 1) * (elex + 1) + 1};
      gridFactory.insertElement(Ikarus::GeometryType::linearQuadrilateral, elementIndices);
    }
  }

  Grid grid = gridFactory.createGrid();

  auto gridView = grid.leafGridView();

  const auto indexSet = gridView.indexSet();

  std::vector<Ikarus::FiniteElements::IFiniteElement> fes;

  for (auto&& ge : surfaces(gridView))
    fes.emplace_back(Ikarus::FiniteElements::ElasticityFE(ge, indexSet, 1000, 0.3));

  auto feManager = Ikarus::FEManager::DefaultFEManager(fes, gridView);

  Ikarus::DirichletConditionManager dirichletConditionManager(feManager);

  dirichletConditionManager.addConstraint(vertices(gridView).front(), 0);
  dirichletConditionManager.addConstraint(vertices(gridView).back(), 1);
  dirichletConditionManager.addConstraint(vertices(gridView).at(3), 1);
  dirichletConditionManager.finalize();

  auto denseMatrixAssembler = Ikarus::Assembler::DenseMatrixAssembler(feManager, dirichletConditionManager);

  for (auto _ : state) {
    Eigen::MatrixXd K;
    benchmark::DoNotOptimize(K= denseMatrixAssembler.getReducedMatrix(Ikarus::FiniteElements::stiffness));
    benchmark::ClobberMemory();
  }
}
BENCHMARK(byhand)->RangeMultiplier(2)->Range(1<<1, 1<<10)->Complexity();

static void eigenIndexing(benchmark::State& state) {
  using namespace Ikarus::Grid;
  using Grid = SimpleGrid<2, 2>;
  SimpleGridFactory<2, 2> gridFactory;
  using vertexType = Eigen::Vector2d;
  std::vector<vertexType> verticesVec;
  const double L    = 10;
  const double h    = 1;
  const size_t elex = state.range(0);
  const size_t eley = state.range(0);
  for (size_t j = 0; j < eley + 1; ++j) {
    for (size_t i = 0; i < elex + 1; ++i)
      verticesVec.emplace_back(vertexType{i * L / (elex), j * h / (eley)});
  }
  for (auto&& vert : verticesVec)
    gridFactory.insertVertex(vert);
  for (size_t i = 0; i < elex; ++i) {
    for (size_t j = 0; j < eley; ++j) {
      std::vector<size_t> elementIndices;
      elementIndices.resize(4);
      elementIndices
          = {i + j * (elex + 1), i + j * (elex + 1) + 1, i + (j + 1) * (elex + 1), i + (j + 1) * (elex + 1) + 1};
      gridFactory.insertElement(Ikarus::GeometryType::linearQuadrilateral, elementIndices);
    }
  }

  Grid grid = gridFactory.createGrid();

  auto gridView = grid.leafGridView();

  const auto indexSet = gridView.indexSet();

  std::vector<Ikarus::FiniteElements::IFiniteElement> fes;

  for (auto&& ge : surfaces(gridView))
    fes.emplace_back(Ikarus::FiniteElements::ElasticityFE(ge, indexSet, 1000, 0.3));

  auto feManager = Ikarus::FEManager::DefaultFEManager(fes, gridView);

  Ikarus::DirichletConditionManager dirichletConditionManager(feManager);

  dirichletConditionManager.addConstraint(vertices(gridView).front(), 0);
  dirichletConditionManager.addConstraint(vertices(gridView).back(), 1);
  dirichletConditionManager.addConstraint(vertices(gridView).at(3), 1);

  auto denseMatrixAssembler = Ikarus::Assembler::DenseMatrixAssembler(feManager, dirichletConditionManager);
  const std::vector<size_t> keepIndices(dirichletConditionManager.freeIndices().begin(),
                                        dirichletConditionManager.freeIndices().end());
  for (auto _ : state) {

    Eigen::MatrixXd K;
    Eigen::MatrixXd Kred;

    benchmark::DoNotOptimize(K=denseMatrixAssembler.getMatrix(Ikarus::FiniteElements::stiffness));
    benchmark::DoNotOptimize(Kred = K(keepIndices, keepIndices));
//    benchmark::DoNotOptimize(K = K(Eigen::all, keepIndices).eval());
    benchmark::ClobberMemory();
  }
}
BENCHMARK(eigenIndexing)->RangeMultiplier(2)->Range(1<<1, 1<<10)->Complexity();
