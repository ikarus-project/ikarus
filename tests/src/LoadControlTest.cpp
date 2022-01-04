//
// Created by Alex on 21.07.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"

#include "spdlog/spdlog.h"

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/Assembler/SimpleAssemblers.h"
#include "ikarus/Controlroutines/LoadControl.h"
#include "ikarus/FEManager/DefaultFEManager.h"
#include "ikarus/FiniteElements/ElasticityFE.h"
#include "ikarus/LinearAlgebra/DirichletConditionManager.h"
#include "ikarus/utils/Observer/controlLogger.h"
#include "ikarus/utils/Observer/gridDrawerObserver.h"
#include <ikarus/FiniteElements/ForceLoad.h>
#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>
#include <ikarus/LinearAlgebra/NonLinearOperator.h>

TEST(LoadControlTest, GridLoadControlTest) {
  using namespace Ikarus::Grid;
  using Grid = SimpleGrid<2, 2>;
  SimpleGridFactory<2, 2> gridFactory;
  using vertexType = Eigen::Vector2d;
  std::vector<vertexType> verticesVec;
  const double L    = 10;
  const double h    = 1;
  const size_t elex = 2;
  const size_t eley = 2;
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

  Grid grid     = gridFactory.createGrid();
  auto gridView = grid.leafGridView();

  //  draw(gridView);

  std::vector<Ikarus::FiniteElements::IFiniteElement> feContainer;

  for (auto&& ge : rootEntities(gridView))
    feContainer.emplace_back(Ikarus::FiniteElements::ElasticityFE(ge, gridView.indexSet(), 1000, 0.0));

  auto spaceFunction = [](const Eigen::Vector2d& v)-> Eigen::Vector2d{Eigen::Vector2d f{}; f[1]=1; return f;};
  feContainer.emplace_back(Ikarus::FiniteElements::ForceLoad(edges(gridView).back(), gridView.indexSet(),spaceFunction));

  auto feManager = Ikarus::FEManager::DefaultFEManager(feContainer, gridView);

  Ikarus::DirichletConditionManager dirichletConditionManager(feManager);
  dirichletConditionManager.addConstraint(vertices(gridView).front(), 0);
  dirichletConditionManager.addConstraint(vertices(gridView).back(), 1);
  dirichletConditionManager.addConstraint(vertices(gridView).at(3), 1);
  dirichletConditionManager.finalize();

  auto vectorAssembler = Ikarus::Assembler::VectorAssembler(feManager, dirichletConditionManager);

  auto denseMatrixAssembler  = Ikarus::Assembler::DenseMatrixAssembler(feManager, dirichletConditionManager);
  auto sparseMatrixAssembler = Ikarus::Assembler::SparseMatrixAssembler(feManager,dirichletConditionManager);

  auto& x = feManager.getVariables();

  [[maybe_unused]] auto fintFunction
      = [&](auto&& lambda) { return vectorAssembler.getReducedVector(Ikarus::FiniteElements::forces, lambda); };
  [[maybe_unused]] auto KFunction
      = [&](auto&& lambda) { return denseMatrixAssembler.getReducedMatrix(Ikarus::FiniteElements::stiffness, lambda); };
  [[maybe_unused]] auto KFunctionSparse
      = [&](auto&& lambda) { return sparseMatrixAssembler.getReducedMatrix(Ikarus::FiniteElements::stiffness, lambda); };

  auto controlObserver = std::make_shared<ControlLogger>();
  //  auto gridDrawerObserver =
  //  std::make_shared<GridDrawerObserver<decltype(gridView),decltype(feManager)>>(gridView,feManager);
  Ikarus::LoadControl lc(feManager,dirichletConditionManager, linearAlgebraFunctions(fintFunction, KFunction), 10,{0,1});
  lc.subscribeAll(controlObserver);
  //  lc.subscribe(ControlMessages::SOLUTION_CHANGED,gridDrawerObserver);
  lc.run();
}