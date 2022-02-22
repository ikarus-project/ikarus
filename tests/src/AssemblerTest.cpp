//
// Created by Alex on 21.04.2021.
//
#include <gmock/gmock.h>

#include "../../config.h"
#include "testHelpers.h"

#include <fstream>
#include <vector>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <Eigen/Core>

#include <ikarus/Assembler/SimpleAssemblers.h>
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/FiniteElements/NonLinearElasticityFEwithBasis.h>

TEST(Assembler, SimpleAssemblersTest) {
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox = {4, 2};
  std::array<int, 2> eles           = {2, 1};
  auto grid                         = std::make_shared<Grid>(bbox, eles);

  for (int i = 0; i < 4; ++i) {
    auto gridView = grid->leafGridView();

    using namespace Dune::Functions::BasisFactory;
    auto basis = makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));

    const auto& indexSet = gridView.indexSet();

    std::vector<Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis<decltype(basis)>> fes;
    const double Emodul = 1000;

    for (auto&& ge : elements(gridView))
      fes.emplace_back(basis, ge, Emodul, 0.3);

    std::vector<bool> dirichFlags(basis.size(), false);

    Dune::Functions::forEachBoundaryDOF(basis, [&](auto&& indexGlobal) { dirichFlags[indexGlobal] = true; });

    Ikarus::SparseFlatAssembler sparseFlatAssembler(basis, fes, dirichFlags);
    Ikarus::DenseFlatAssembler denseFlatAssembler(basis, fes, dirichFlags);

    Ikarus::FErequirements req;
    Eigen::VectorXd d(basis.size());
    d.setZero();
    req.sols.emplace_back(d);
    req.parameter.insert({Ikarus::FEParameter::loadfactor, 0});
    req.matrixAffordances = Ikarus::MatrixAffordances::stiffness;
    auto& Kdense          = denseFlatAssembler.getMatrix(req);
    auto& K               = sparseFlatAssembler.getMatrix(req);
//    std::cout << K << std::endl;
//    std::cout << Kdense << std::endl;
    EXPECT_THAT(K, EigenApproxEqual(Kdense, 1e-15));
    EXPECT_THAT(K.rows(), 2*gridView.size(2));
    EXPECT_THAT(K.cols(), 2*gridView.size(2));
    const int boundaryNodes = (eles[0]*Dune::power(2,i)+1)*2+(eles[1]*Dune::power(2,i)+1)*2-4;
    EXPECT_EQ(std::ranges::count(dirichFlags,true), 2*boundaryNodes);
    grid->globalRefine(1);
  }
  //
  //  Ikarus::DirichletConditionManager dirichletConditionManager(feManager);
  //
  //  dirichletConditionManager.addConstraint(*vertices(gridView).begin(), 0);
  //  dirichletConditionManager.finalize();
  //
  //  auto vectorAssembler = Ikarus::Assembler::VectorAssembler(feManager, dirichletConditionManager);
  //  auto fint            = vectorAssembler.getVector(Ikarus::FiniteElements::forces);
  //  EXPECT_EQ(fint.size(), 12);
  //  const Eigen::VectorXd fintExpected = (Eigen::VectorXd(12) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0).finished();
  //  EXPECT_THAT(fint, EigenApproxEqual(fintExpected, 1e-15));
  //
  //  auto denseMatrixAssembler = Ikarus::Assembler::DenseMatrixAssembler(feManager, dirichletConditionManager);
  //  auto K                    = denseMatrixAssembler.getMatrix(Ikarus::FiniteElements::stiffness);
  //  EXPECT_EQ(K.rows(), 12);
  //  EXPECT_EQ(K.cols(), 12);
  //  const Eigen::MatrixXd KExpected
  //      = (Eigen::MatrixXd(12, 12) <<  // clang-format off
  //  494.5054945054946,  178.5714285714286, -302.1978021978022,
  //  -13.73626373626373,  54.94505494505495,  13.73626373626373, -247.2527472527473, -178.5714285714286, 0, 0, 0, 0,
  //  178.5714285714286,  494.5054945054946,  13.73626373626373,  54.94505494505495, -13.73626373626373,
  //  -302.1978021978022, -178.5714285714286, -247.2527472527473,                  0,                  0, 0, 0,
  // -302.1978021978023,  13.73626373626372,  989.0109890109891,                  0, -247.2527472527473,
  // 178.5714285714286,  109.8901098901099, -7.10542735760e-15, -302.1978021978022, -13.73626373626373,
  // -247.2527472527473, -178.5714285714286, -13.73626373626373,  54.94505494505495,                  0,
  // 989.0109890109891,  178.5714285714286, -247.2527472527473, -3.55271367880e-15,
  // -604.3956043956044,  13.73626373626373,  54.94505494505495, -178.5714285714286, -247.2527472527473,
  //  54.94505494505495, -13.73626373626373, -247.2527472527473,  178.5714285714286,  494.5054945054945,
  //  -178.5714285714286, -302.1978021978022,  13.73626373626374,                  0,                  0, 0, 0,
  //  13.73626373626373, -302.1978021978023,  178.5714285714286, -247.2527472527473, -178.5714285714286,
  //  494.5054945054945, -13.73626373626373,  54.94505494505495,                  0,                  0, 0, 0,
  // -247.2527472527473, -178.5714285714286,  109.8901098901099, -2.48689957516e-14, -302.1978021978022,
  // -13.73626373626373,  989.0109890109891,                  0, -247.2527472527473,  178.5714285714286,
  // -302.1978021978022,  13.73626373626374, -178.5714285714286, -247.2527472527473, -1.06581410364e-14,
  // -604.3956043956046,  13.73626373626373,  54.94505494505495,                  0,  989.0109890109891,
  // 178.5714285714286, -247.2527472527473, -13.73626373626372,  54.94505494505496,
  //                  0,                  0, -302.1978021978023,  13.73626373626373,                  0, 0,
  //                  -247.2527472527473,  178.5714285714286,  494.5054945054945,
  //                  -178.5714285714286,  54.94505494505495, -13.73626373626373, 0,                  0,
  //                  -13.73626373626373,  54.94505494505495,                  0,                  0, 178.5714285714286,
  //                  -247.2527472527473, -178.5714285714286,  494.5054945054945,  13.73626373626374,
  //                  -302.1978021978022, 0,                  0, -247.2527472527473, -178.5714285714286, 0, 0,
  //                  -302.1978021978022, -13.73626373626372,  54.94505494505495,  13.73626373626373, 494.5054945054945,
  //                  178.5714285714286, 0,                  0, -178.5714285714286, -247.2527472527473, 0,
  //                  0,  13.73626373626373,  54.94505494505496, -13.73626373626373, -302.1978021978022,
  //                  178.5714285714286,  494.5054945054945).finished();  // clang-format on
  //  EXPECT_THAT(K, EigenApproxEqual(KExpected, 1e-15));
  //
  //  auto sparseMatrixAssembler = Ikarus::Assembler::SparseMatrixAssembler(feManager, dirichletConditionManager);
  //  auto KSparse               = sparseMatrixAssembler.getMatrix(Ikarus::FiniteElements::stiffness);
  //
  //  EXPECT_THAT(KSparse, EigenApproxEqual(KExpected, 1e-15));
  //  EXPECT_THAT(KSparse, EigenApproxEqual(K, 1e-15));
  //
  //  auto scalarAssembler = Ikarus::Assembler::ScalarAssembler(feManager);
  //  auto w               = scalarAssembler.getScalar(Ikarus::FiniteElements::potentialEnergy);
  //  EXPECT_DOUBLE_EQ(w, 0.0);
  //
  //  // Reduced tests
  //  const auto& fintRed = vectorAssembler.getReducedVector(Ikarus::FiniteElements::forces);
  //
  //  EXPECT_EQ(fintRed.size(), fint.size() - 3);
  //
  //  Eigen::VectorXd testVector = Eigen::VectorXd::LinSpaced(9, 1, 9);
  //  const auto fullVector      = vectorAssembler.createFullVector(testVector);
  //
  //  EXPECT_EQ(fullVector.size(), fint.size());
  //  const Eigen::VectorXd fullVectorExpected = (Eigen::VectorXd(12) << 0, 1, 2, 3, 4, 5, 6, 0, 7, 8, 9, 0).finished();
  //  EXPECT_THAT(fullVector, EigenApproxEqual(fullVectorExpected, 1e-15));
  //
  //  const auto& KRed = denseMatrixAssembler.getReducedMatrix(Ikarus::FiniteElements::stiffness);
  //
  //  Eigen::MatrixXd KExpectedRed = KExpected;
  //  std::vector<size_t> keepIndices(dirichletConditionManager.freeIndices().begin(),
  //                                  dirichletConditionManager.freeIndices().end());
  //
  //  KExpectedRed = KExpectedRed(keepIndices, keepIndices).eval();
  //
  //  EXPECT_THAT(KRed, EigenApproxEqual(KExpectedRed, 1e-15));
  //
  //  auto& x = feManager.getVariables();
  //  Eigen::VectorXd D(feManager.numberOfDegreesOfFreedom());
  //  D.setOnes();
  //  x += D;
  //
  //  w = scalarAssembler.getScalar(Ikarus::FiniteElements::potentialEnergy);
  //  EXPECT_NEAR(w, 0.0, 1e-16 * Emodul);
  //
  //  const auto& KRedSparse = sparseMatrixAssembler.getReducedMatrix(Ikarus::FiniteElements::stiffness);
  //
  //  EXPECT_THAT(KRedSparse, EigenApproxEqual(KExpectedRed, 1e-15));
}
