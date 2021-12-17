//
// Created by Alex on 21.07.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/Assembler/SimpleAssemblers.h"
#include "ikarus/FEManager/DefaultFEManager.h"
#include "ikarus/FiniteElements/ElasticityFE.h"
#include <ikarus/FiniteElements/ForceLoad.h>
#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/LinearAlgebra/NonLinearOperator.h>

#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>

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

  feContainer.emplace_back(Ikarus::FiniteElements::ForceLoad(*(edges(gridView).end()), gridView.indexSet()));

  auto feManager = Ikarus::FEManager::DefaultFEManager(feContainer, gridView);

  auto vectorAssembler = Ikarus::Assembler::VectorAssembler(feManager);

  auto denseMatrixAssembler  = Ikarus::Assembler::DenseMatrixAssembler(feManager);
  auto sparseMatrixAssembler = Ikarus::Assembler::SparseMatrixAssembler(feManager);

  auto& x = feManager.getVariables();

  const auto fintFunction = [&]() { return vectorAssembler.getVector(Ikarus::FiniteElements::forces); };
  auto KFunction          = [&]() { return denseMatrixAssembler.getMatrix(Ikarus::FiniteElements::stiffness); };
  auto KFunctionSparse    = [&]() { return sparseMatrixAssembler.getMatrix(Ikarus::FiniteElements::stiffness); };
  Ikarus::NonLinearOperator nonLinearOperator(fintFunction, derivatives(KFunction), parameter());
  Ikarus::NonLinearOperator nonLinearOperatorWithSparseMatrix(fintFunction, derivatives(KFunctionSparse), parameter());

// Ikarus::LoadControl lc(nonLinearOperator)
}