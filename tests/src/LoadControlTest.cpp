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
  const size_t elex = 5;
  const size_t eley = 5;
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

  feContainer.emplace_back(Ikarus::FiniteElements::ForceLoad(*(edges(gridView).end() - 1), gridView.indexSet()));

  auto feManager = Ikarus::FEManager::DefaultFEManager(feContainer, gridView);

  Ikarus::DirichletConditionManager dirichletConditionManager(feManager);

  auto vectorAssembler = Ikarus::Assembler::VectorAssembler(feManager, dirichletConditionManager);

  auto denseMatrixAssembler  = Ikarus::Assembler::DenseMatrixAssembler(feManager, dirichletConditionManager);
  auto sparseMatrixAssembler = Ikarus::Assembler::SparseMatrixAssembler(feManager,dirichletConditionManager);

  auto& x = feManager.getVariables();

  //  auto lambdaTest = Ikarus::FEParameterFactory::createParameter(Ikarus::FEParameter::loadfactor,1);

  [[maybe_unused]] auto fintFunction
      = [&](auto&& lambda) { return vectorAssembler.getVector(Ikarus::FiniteElements::forces, lambda); };
  [[maybe_unused]] auto KFunction
      = [&](auto&& lambda) { return denseMatrixAssembler.getMatrix(Ikarus::FiniteElements::stiffness, lambda); };
  [[maybe_unused]] auto KFunctionSparse
      = [&](auto&& lambda) { return sparseMatrixAssembler.getMatrix(Ikarus::FiniteElements::stiffness, lambda); };
  //  Ikarus::NonLinearOperator nonLinearOperator(fintFunction, derivatives(KFunction), parameter());
  //  Ikarus::NonLinearOperator nonLinearOperatorWithSparseMatrix(fintFunction, derivatives(KFunctionSparse),
  //  parameter());
  auto controlObserver = std::make_shared<ControlLogger>();
  //  auto gridDrawerObserver =
  //  std::make_shared<GridDrawerObserver<decltype(gridView),decltype(feManager)>>(gridView,feManager);
  Ikarus::LoadControl lc(feManager, linearAlgebraFunctions(fintFunction, KFunction), 10);
  lc.subscribeAll(controlObserver);
  //  lc.subscribe(ControlMessages::SOLUTION_CHANGED,gridDrawerObserver);
  lc.run();
}