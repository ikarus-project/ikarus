//
// Created by Alex on 21.07.2021.
//
#include <gmock/gmock.h>

#include "testHelpers.h"
#define EIGEN_SPARSEMATRIX_PLUGIN "eigenSparseAddon.h"

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include "common.h"

#include "spdlog/spdlog.h"

#include <Eigen/Core>

#include "ikarus/Controlroutines/LoadControl.h"
#include "ikarus/FiniteElements/Mechanics/NonLinearElasticityFE.h"
#include "ikarus/Solver/NonLinearSolver/NewtonRaphson.hpp"
#include "ikarus/Solver/NonLinearSolver/TrustRegion.hpp"
#include "ikarus/utils/Observer/controlVTKWriter.h"
#include <ikarus/Assembler/SimpleAssemblers.h>
#include "ikarus/utils/drawing/griddrawer.h"
#include <ikarus/LinearAlgebra/NonLinearOperator.h>
#include <ikarus/utils/utils/algorithms.h>



template <typename T>
class NonLinearElasticityLoadControlNRandTR : public testing::Test {
public:
  using GridId = T;
  NonLinearElasticityLoadControlNRandTR() : value_{createGrid<T>()} {}
  decltype(createGrid<T>()) value_;
};
using GridTypes = ::testing::Types<Grids::Yasp, Grids::Alu, Grids::Iga>;
// using GridTypes = ::testing::Types<Grids::Yasp>;

TYPED_TEST_SUITE(NonLinearElasticityLoadControlNRandTR, GridTypes);

TYPED_TEST(NonLinearElasticityLoadControlNRandTR, ComputeMaxDisp) {
  auto gridView = this->value_->leafGridView();

  using namespace Ikarus;

  using namespace Dune::Functions::BasisFactory;
  //  auto basis = makeBasis(gridView, power<gridDim>(gridView.getPreBasis(), FlatInterleaved()));
  auto basis = makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));

  auto localView = basis.localView();
  std::vector<Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis<decltype(basis)>> fes;
  auto volumeLoad = [](auto& globalCoord, auto& lamb) {
    Eigen::Vector2d fext;
    fext.setZero();
    fext[1] = 2 * lamb;
    fext[0] = lamb;
    return fext;
  };
  for (auto& element : elements(gridView))
    fes.emplace_back(basis, element, 1000, 0.3, volumeLoad);

  std::vector<bool> dirichletFlags(basis.size(), false);

  Dune::Functions::forEachBoundaryDOF(basis, [&](auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[1]) < 1e-8) {
      dirichletFlags[localView.index(localIndex)[0]] = true;
    }
  });

  auto denseAssembler  = DenseFlatSimpleAssembler(basis, fes, dirichletFlags);
  auto sparseAssembler = SparseFlatAssembler(basis, fes, dirichletFlags);

  Eigen::VectorXd d;
  d.setZero(basis.size());
  double lambda = 0.0;

  auto residualFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    return denseAssembler.getVector(forces, disp, lambdaLocal);
  };

  auto KFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements<Eigen::VectorXd> req;
    req.sols.emplace_back(disp);
    req.parameter.insert({Ikarus::FEParameter::loadfactor, lambdaLocal});
    req.matrixAffordances = Ikarus::MatrixAffordances::stiffness;
    return sparseAssembler.getMatrix(req);
  };

  auto energyFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto {
    return denseAssembler.getScalar(potentialEnergy, disp_, lambdaLocal);
  };

  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(energyFunction, residualFunction, KFunction),
                                            parameter(d, lambda));

  const double gradTol = 1e-8;

  auto tr = Ikarus::makeTrustRegion(nonLinOp);
  tr->setup({.verbosity = 1,
             .maxiter   = 1000,
             .grad_tol  = gradTol,
             .corr_tol  = 1e-16,  // everything should converge to the gradient tolerance
             .useRand   = false,
             .rho_reg   = 1e8,
             .Delta0    = 1});

  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<decltype(basis)>>(basis, d, 2);
  vtkWriter->setFileNamePrefix("Test2Dsolid");
  vtkWriter->setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 2);

  auto lc = Ikarus::LoadControl(tr, 1, {0, 2000});
  lc.subscribeAll(vtkWriter);
  const auto controlInfo = lc.run();
  nonLinOp.template update<0>();
  const auto maxDisp = std::ranges::max(d);
  if constexpr (std::is_same_v<TypeParam, Grids::Yasp>) {
    EXPECT_DOUBLE_EQ(nonLinOp.value(), -1.4809559783564966e+03);
    EXPECT_NEAR(maxDisp, 0.786567027108460048,1e-12);
  } else if constexpr (std::is_same_v<TypeParam, Grids::Alu>) {
    EXPECT_NEAR(nonLinOp.value(), -1.4842107484533601e+03,1e-12);
    EXPECT_NEAR(maxDisp, 0.78426066482258983,1e-15);
  } else if constexpr (std::is_same_v<TypeParam, Grids::Iga>) {
    EXPECT_NEAR(nonLinOp.value(), -8.1142552237939071e+02,1e-12);
    EXPECT_NEAR(maxDisp, 0.615624125459537153,1e-15);
  }

  nonLinOp.template update<1>();
  EXPECT_TRUE(controlInfo.sucess);
  EXPECT_LE(nonLinOp.derivative().norm(), gradTol);
}