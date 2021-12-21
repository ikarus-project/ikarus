//
// Created by Alex on 21.07.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"
#define EIGEN_SPARSEMATRIX_PLUGIN "eigenSparseAddon.h"
#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/Assembler/SimpleAssemblers.h"
#include "ikarus/FEManager/DefaultFEManager.h"
#include "ikarus/FiniteElements/ElasticityFE.h"
#include "ikarus/Grids/GridHelper/griddrawer.h"
#include "ikarus/LinearAlgebra/DirichletConditionManager.h"
#include <ikarus/LinearAlgebra/NonLinearOperator.h>

auto f(double& x) { return 0.5 * x * x + x - 2; }
auto df(double& x) { return x + 1; }

TEST(NonLinearOperator, SimpleOperator) {
  double x = 13;

  auto fvLambda  = [&](auto&& x) { return f(x); };
  auto dfvLambda = [&](auto&& x) { return df(x); };
  Ikarus::NonLinearOperator nonLinOp(linearAlgebraFunctions(fvLambda,dfvLambda), parameter(x));

  auto& val      = nonLinOp.value();
  auto& gradient = nonLinOp.derivative();

  // Newton method test
  const double eps  = 1e-14;
  const int maxIter = 20;
  int iter          = 0;
  while (val > eps && iter < maxIter) {
    x -= val / gradient;
    nonLinOp.updateAll();
    ++iter;
  }
  EXPECT_DOUBLE_EQ(val, 0.0);
  double xExpected = std::sqrt(5.0) - 1.0;
  EXPECT_DOUBLE_EQ(gradient, df(xExpected));
  EXPECT_DOUBLE_EQ(x, xExpected);
}

Eigen::VectorXd fv(Eigen::VectorXd& x, Eigen::MatrixXd& A, Eigen::VectorXd& b) { return b + A * x; }
Eigen::MatrixXd dfv([[maybe_unused]] Eigen::VectorXd& x, Eigen::MatrixXd& A, [[maybe_unused]] Eigen::VectorXd& b) {
  return A;
}

TEST(NonLinearOperator, VectorValuedOperator) {
  Eigen::VectorXd x(3);

  x << 1, 2, 3;
  Eigen::VectorXd b(3);
  b << 5, 7, 8;
  Eigen::MatrixXd A(3, 3);
  A = Eigen::MatrixXd::Identity(3, 3) * 13;

  auto fvLambda  = [&](auto&& x) { return fv(x, A, b); };
  auto dfvLambda = [&](auto&& x) { return dfv(x, A, b); };
  auto nonLinOp  = Ikarus::NonLinearOperator(linearAlgebraFunctions(fvLambda,dfvLambda), parameter(x));

  auto& val      = nonLinOp.value();
  auto& jacobian = nonLinOp.derivative();

  // Newton method test
  const double eps  = 1e-14;
  const int maxIter = 20;
  int iter          = 0;
  while (val.norm() > eps && iter < maxIter) {
    x -= jacobian.inverse() * val;
    nonLinOp.updateAll();
    ++iter;
  }
  EXPECT_EQ(iter, 1);  // Linear System should be solved in one step
  EXPECT_THAT(b, EigenApproxEqual(-A * x, 1e-15));
}

double f2v(Eigen::VectorXd& x, Eigen::MatrixXd& A, Eigen::VectorXd& b) { return x.dot(b + A * x); }
Eigen::VectorXd df2v([[maybe_unused]] Eigen::VectorXd& x, Eigen::MatrixXd& A, [[maybe_unused]] Eigen::VectorXd& b) {
  return 2 * A * x + b;
}
Eigen::MatrixXd ddf2v([[maybe_unused]] Eigen::VectorXd& x, Eigen::MatrixXd& A, [[maybe_unused]] Eigen::VectorXd& b) {
  return 2 * A;
}

TEST(NonLinearOperator, SecondOrderVectorValuedOperator) {
  Eigen::VectorXd x(3);

  x << 1, 2, 3;
  Eigen::VectorXd b(3);
  b << 5, 7, 8;
  Eigen::MatrixXd A(3, 3);
  A = Eigen::MatrixXd::Identity(3, 3) * 13;

  auto fvLambda   = [&](auto&& x) { return f2v(x, A, b); };
  auto dfvLambda  = [&](auto&& x) { return df2v(x, A, b); };
  auto ddfvLambda = [&](auto&& x) { return ddf2v(x, A, b); };
  auto nonLinOp   = Ikarus::NonLinearOperator(linearAlgebraFunctions(fvLambda, dfvLambda, ddfvLambda), parameter(x));

  auto& val      = nonLinOp.value();
  auto& residual = nonLinOp.derivative();
  auto& hessian  = nonLinOp.secondDerivative();

  // Newton method test find root of first derivative
  const double eps  = 1e-14;
  const int maxIter = 20;
  int iter          = 0;
  while (residual.norm() > eps && iter < maxIter) {
    x -= hessian.inverse() * residual;
    nonLinOp.updateAll();
    ++iter;
  }

  EXPECT_EQ(iter, 1);                   // Linear System should be solved in one step
  EXPECT_EQ(val, -2.6538461538461533);  // Linear System should be solved in one step
  EXPECT_THAT(b, EigenApproxEqual(-2 * A * x, 1e-15));
}

#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>

TEST(NonLinearOperator, GridLoadControlTest) {
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

  auto feManager = Ikarus::FEManager::DefaultFEManager(feContainer, gridView);
  Ikarus::DirichletConditionManager dirichletConditionManager(feManager);
  auto vectorAssembler = Ikarus::Assembler::VectorAssembler(feManager,dirichletConditionManager);

  auto denseMatrixAssembler  = Ikarus::Assembler::DenseMatrixAssembler(feManager,dirichletConditionManager);
  auto sparseMatrixAssembler = Ikarus::Assembler::SparseMatrixAssembler(feManager);

  auto& x = feManager.getVariables();

  auto fintFunction = [&]() { return vectorAssembler.getVector(Ikarus::FiniteElements::forces); };
  auto KFunction          = [&]() { return denseMatrixAssembler.getMatrix(Ikarus::FiniteElements::stiffness); };
  auto KFunctionSparse    = [&]() { return sparseMatrixAssembler.getMatrix(Ikarus::FiniteElements::stiffness); };
  Ikarus::NonLinearOperator nonLinearOperator(linearAlgebraFunctions(fintFunction, KFunction), parameter());
  Ikarus::NonLinearOperator nonLinearOperatorWithSparseMatrix(linearAlgebraFunctions(fintFunction, KFunctionSparse), parameter());

  auto& K             = nonLinearOperator.derivative();
  auto& Ksparse = nonLinearOperatorWithSparseMatrix.derivative();
  EXPECT_THAT(Ksparse, EigenApproxEqual(K, 1e-15));

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> decomp(K);
  auto rank = decomp.rank();
  EXPECT_EQ(rank, K.cols() - 3);
  EXPECT_EQ(rank, K.rows() - 3);
  for (int i = 1; i < K.cols(); ++i) {
    K.coeffRef(i, 0) = 0.0;
    K.coeffRef(i, 1) = 0.0;
    K.coeffRef(i, 3) = 0.0;
    K.coeffRef(0, i) = 0.0;
    K.coeffRef(1, i) = 0.0;
    K.coeffRef(3, i) = 0.0;
  }
  K.coeffRef(0, 0) = 1;
  K.coeffRef(1, 1) = 1;
  K.coeffRef(3, 3) = 1;
  decomp.compute(K);
  rank = decomp.rank();
  EXPECT_EQ(rank, K.cols());
  EXPECT_EQ(rank, K.rows());

  auto f = nonLinearOperator.value();
  f.setZero();
  *(f.end() - 1) = 1;

  const Eigen::VectorXd D = K.ldlt().solve(f);
  x += D;

  EXPECT_EQ(x.size(), K.rows());

  const Eigen::VectorXd DExpected
      = (Eigen::VectorXd(D.size()) <<  // clang-format off
                     0,                     0,    0.0164851890835871,                     0,  0.001965773386162552,
 -0.001667201218069932,   0.01128335946342303,  0.001162955868264214,  0.003066026312401127, -0.002957921126533793,
  0.005928558345036764,  0.002053086270565305,  0.003749546029349699, -0.003843850843882211, 0.0005055017822511749,
  0.002661296000886878,  0.004438801077748433, -0.004337313972296621, -0.004889302054841124,  0.003001103912624741,
  0.005544771996353223, -0.004490409434796363,  -0.01015635854765239,  0.003107888488213445,   0.04468173607972595,
     0.116003062578087,   0.02790239873679473,    0.1157057299171808,   0.01136310906385944,     0.115495641419241,
 -0.005065813560613874,    0.1153566989681915,  -0.02150698196808424,    0.1152764645623072,  -0.03808209501596624,
     0.115246427550179,   0.06451829087213176,    0.3248432838240538,   0.03997298490470038,    0.3249186967747374,
   0.01546750172434345,    0.3249672853269041, -0.009033580351837775,    0.3250016827881556,   -0.0335738164890635,
    0.3250298818791605,  -0.05819940484074704,     0.325055450585741,   0.07657130165145397,     0.598521780858297,
   0.04714486190324706,     0.598495791002132,   0.01783386607777375,    0.5984646626256116,  -0.01144017663905125,
    0.5984192164801975,  -0.04074972562032137,    0.5983549304996514,  -0.07016388748708198,    0.5982716943882292,
   0.08049911555092175,    0.9036056530449433,   0.04954760592999906,    0.9036427616592595,   0.01866785645690688,
     0.903738477407708,  -0.01220627723080757,    0.9039079699348371,  -0.04314837701125365,    0.9041544598306273,
  -0.07423566623294452,    0.9044693949758897).finished();  // clang-format on

  EXPECT_THAT(D, EigenApproxEqual(DExpected, 1e-12));
  //  drawDeformed(gridView, feManager);
}