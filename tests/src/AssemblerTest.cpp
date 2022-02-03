//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>

#include "testHelpers.h"

#include <fstream>
#include <vector>

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <Eigen/Core>

#include "ikarus/LinearAlgebra/DirichletConditionManager.h"
#include <ikarus/Assembler/SimpleAssemblers.h>
#include <ikarus/FEManager/DefaultFEManager.h>
#include <ikarus/FiniteElements/ElasticityFE.h>
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Geometries/GeometryType.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>






TEST(Assembler, SimpleAssemblersTest) {
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox = {4, 2};
  std::array<int, 2> eles           = {2, 1};
  auto grid                         = std::make_shared<Grid>(bbox, eles);

  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView,  power<2>(lagrange<1>(), FlatInterleaved()));

  const auto& indexSet = gridView.indexSet();
  auto localView = basis.localView();
  std::vector<Ikarus::FiniteElements::ElasticityFE<decltype(localView)>> fes;
  const double Emodul = 1000;

  for (auto&& ge : elements(gridView)) {
    localView.bind(ge);
    fes.emplace_back(Ikarus::FiniteElements::ElasticityFE(localView, Emodul, 0.3));
  }
  auto basis1Factory = lagrange<1>();
  auto basis2Factory = lagrange<1>();

  auto basis2 = makeBasis(gridView, composite(power<2>(lagrange<1>(),FlatInterleaved()),power<2>(lagrange<1>(),FlatInterleaved()),FlatLexicographic()));

  using namespace Dune::Indices;

  auto localView2 = basis2.localView();

  auto dispBasis = subspaceBasis(basis2,_0);
  auto pressureBasis = subspaceBasis(basis2,_1);
  auto localViewOfDisplacement = dispBasis.localView();
  auto localViewOfPressure = pressureBasis.localView();

  for (auto&& ge : elements(gridView)) {
    localView2.bind(ge);
    std::cout<<"===="<<std::endl;
    for (auto i = 0U; i < localView2.size(); ++i) {
      std::cout<<localView2.index(i)<<std::endl;
    }
    std::cout<<"==D=="<<std::endl;
    localViewOfDisplacement.bind(ge);
    for (size_t i = 0; i < localViewOfDisplacement.tree().child(0).finiteElement().size(); ++i) {
      for (size_t j = 0; j < localViewOfDisplacement.tree().CHILDREN; ++j) {
        std::cout<<localViewOfDisplacement.index(localViewOfDisplacement.tree().child(j).localIndex(i))<<" ";
      }
      std::cout<<std::endl;
    }
    std::cout<<"==P=="<<std::endl;
    localViewOfPressure.bind(ge);
    const auto& pressureFe=localViewOfPressure.tree().child(0).finiteElement();
    for (size_t i = 0; i < pressureFe.size(); ++i) {
      for (size_t j = 0; j < localViewOfPressure.tree().CHILDREN; ++j)
        std::cout<<localViewOfPressure.index(localViewOfPressure.tree().child(j).localIndex(i))  <<" ";
      std::cout<<std::endl;
    }
    std::cout<<"==End=="<<std::endl;

    const int nDofs0 = localView2.tree().child(_0,0).finiteElement().size();
    const int nDofs1 = localView2.tree().child(_1,0).finiteElement().size();
    for (int i=0; i<nDofs0+nDofs1; i++)
    {
      int localIndexI = 0;
      if (i < nDofs0) {
        auto& node = localView2.tree().child(_0,0);
        localIndexI = node.localIndex(i);
      } else {
        auto& node = localView2.tree().child(_1,0);
        localIndexI = node.localIndex(i);

      }
      auto multiIndex = localView2.index(localIndexI);
      std::cout<<multiIndex<<std::endl;
      //CompositeBasis number is contained in multiIndex[0], the Subspacebasis is contained in multiIndex[2]
      //multiIndex[1] contains the actual index
//      if (multiIndex[0] == 0)
//        localConfiguration0[i] = configuration0[multiIndex[1]];
//      else if (multiIndex[0] == 1)
//        localConfiguration1[i-nDofs0] = configuration1[multiIndex[1]];

    }

  }


//  Ikarus::DefaultFEManager feManager(fes,gridView)

//  std::cout<< decltype(localView2)<<std::endl;

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
//  494.5054945054946,  178.5714285714286, -302.1978021978022, -13.73626373626373,  54.94505494505495,  13.73626373626373, -247.2527472527473, -178.5714285714286,                  0,                  0,                  0,                  0,
//  178.5714285714286,  494.5054945054946,  13.73626373626373,  54.94505494505495, -13.73626373626373, -302.1978021978022, -178.5714285714286, -247.2527472527473,                  0,                  0,                  0,                  0,
// -302.1978021978023,  13.73626373626372,  989.0109890109891,                  0, -247.2527472527473,  178.5714285714286,  109.8901098901099, -7.10542735760e-15, -302.1978021978022, -13.73626373626373, -247.2527472527473, -178.5714285714286,
// -13.73626373626373,  54.94505494505495,                  0,  989.0109890109891,  178.5714285714286, -247.2527472527473, -3.55271367880e-15, -604.3956043956044,  13.73626373626373,  54.94505494505495, -178.5714285714286, -247.2527472527473,
//  54.94505494505495, -13.73626373626373, -247.2527472527473,  178.5714285714286,  494.5054945054945, -178.5714285714286, -302.1978021978022,  13.73626373626374,                  0,                  0,                  0,                  0,
//  13.73626373626373, -302.1978021978023,  178.5714285714286, -247.2527472527473, -178.5714285714286,  494.5054945054945, -13.73626373626373,  54.94505494505495,                  0,                  0,                  0,                  0,
// -247.2527472527473, -178.5714285714286,  109.8901098901099, -2.48689957516e-14, -302.1978021978022, -13.73626373626373,  989.0109890109891,                  0, -247.2527472527473,  178.5714285714286, -302.1978021978022,  13.73626373626374,
// -178.5714285714286, -247.2527472527473, -1.06581410364e-14, -604.3956043956046,  13.73626373626373,  54.94505494505495,                  0,  989.0109890109891,  178.5714285714286, -247.2527472527473, -13.73626373626372,  54.94505494505496,
//                  0,                  0, -302.1978021978023,  13.73626373626373,                  0,                  0, -247.2527472527473,  178.5714285714286,  494.5054945054945, -178.5714285714286,  54.94505494505495, -13.73626373626373,
//                  0,                  0, -13.73626373626373,  54.94505494505495,                  0,                  0,  178.5714285714286, -247.2527472527473, -178.5714285714286,  494.5054945054945,  13.73626373626374, -302.1978021978022,
//                  0,                  0, -247.2527472527473, -178.5714285714286,                  0,                  0, -302.1978021978022, -13.73626373626372,  54.94505494505495,  13.73626373626373,  494.5054945054945,  178.5714285714286,
//                  0,                  0, -178.5714285714286, -247.2527472527473,                  0,                  0,  13.73626373626373,  54.94505494505496, -13.73626373626373, -302.1978021978022,  178.5714285714286,  494.5054945054945).finished();  // clang-format on
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
