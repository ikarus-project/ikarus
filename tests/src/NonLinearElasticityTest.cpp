//
// Created by Alex on 21.07.2021.
//
#include <gmock/gmock.h>
#include "../../config.h"

#include "testHelpers.h"
#define EIGEN_SPARSEMATRIX_PLUGIN "eigenSparseAddon.h"
#include <dune/alugrid/grid.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/iga/nurbsgrid.hh>

#include "spdlog/spdlog.h"

#include <Eigen/Core>

#include "ikarus/Controlroutines/LoadControl.h"
#include "ikarus/FiniteElements/NonLinearElasticityFEwithBasis.h"
#include "ikarus/Solver/NonLinearSolver/NewtonRaphson.hpp"
#include "ikarus/Solver/NonLinearSolver/TrustRegion.hpp"
#include "ikarus/basis/basishelper.h"
#include "ikarus/utils/Observer/controlVTKWriter.h"
#include <ikarus/Assembler/SimpleAssemblers.h>
#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/LinearAlgebra/NonLinearOperator.h>
#include <ikarus/utils/utils/algorithms.h>

namespace Grids {
  struct Yasp {};
  struct Alu {};
  struct Iga {};
}  // namespace Grids

template <typename GridType>
auto createGrid() {
  //  //  /// ALUGrid Example
  if constexpr (std::is_same_v<GridType, Grids::Alu>) {
    using Grid = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
    auto grid  = Dune::GmshReader<Grid>::read("../../tests/src/testFiles/unstructuredTrianglesfine.msh", false);
    grid->globalRefine(0);
    return grid;
  } else if constexpr (std::is_same_v<GridType, Grids::Yasp>) {
    using Grid        = Dune::YaspGrid<2>;
    const double L    = 1;
    const double h    = 1;
    const size_t elex = 10;
    const size_t eley = 10;

    Dune::FieldVector<double, 2> bbox = {L, h};
    std::array<int, 2> eles           = {elex, eley};
    auto grid                         = std::make_shared<Grid>(bbox, eles);
    return grid;
  } else if constexpr (std::is_same_v<GridType, Grids::Iga>) {
    constexpr auto dimworld        = 2;
    const std::array<int, 2> order = {2, 2};

    const std::array<std::vector<double>, 2> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};

    using ControlPoint = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointType;

    const std::vector<std::vector<ControlPoint>> controlPoints
        = {{{.p = {0, 0}, .w = 5}, {.p = {0.5, 0}, .w = 1}, {.p = {1, 0}, .w = 1}},
           {{.p = {0, 0.5}, .w = 1}, {.p = {0.5, 0.5}, .w = 10}, {.p = {1, 0.5}, .w = 1}},
           {{.p = {0, 1}, .w = 1}, {.p = {0.5, 1}, .w = 1}, {.p = {1, 1}, .w = 1}}};

    std::array<int, 2> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

    auto controlNet = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointNetType(dimsize, controlPoints);
    using Grid      = Dune::IGA::NURBSGrid<2, dimworld>;

    Dune::IGA::NURBSPatchData<2, dimworld> patchData;
    patchData.knotSpans     = knotSpans;
    patchData.degree        = order;
    patchData.controlPoints = controlNet;
    auto grid               = std::make_shared<Grid>(patchData);
    grid->globalRefine(1);
    return grid;
  }
}

template <typename T>
class NonLinearElasticityLoadControlNRandTR : public testing::Test {
public:
  using GridId = T;
  NonLinearElasticityLoadControlNRandTR() : value_{createGrid<T>()} {}
  decltype(createGrid<T>()) value_;
};
using GridTypes = ::testing::Types<Grids::Yasp, Grids::Alu, Grids::Iga>;
//using GridTypes = ::testing::Types<Grids::Yasp>;

TYPED_TEST_SUITE(NonLinearElasticityLoadControlNRandTR, GridTypes);

TYPED_TEST(NonLinearElasticityLoadControlNRandTR, ComputeMaxDisp) {
  auto gridView = this->value_->leafGridView();

  using namespace Ikarus;

  using namespace Dune::Functions::BasisFactory;
  //  auto basis = makeBasis(gridView, power<gridDim>(gridView.getPreBasis(), FlatInterleaved()));
  auto basis = makeBasis(gridView, power<2>(lagrange<1>(), FlatInterleaved()));

  auto localView = basis.localView();
  std::vector<Ikarus::FiniteElements::NonLinearElasticityFEWithLocalBasis<decltype(basis)>> fes;
  auto volumeLoad = [](auto& globalCoord, auto& lamb)
  {Eigen::Vector2d fext;
    fext.setZero();
    fext[1] = 2 * lamb;
    fext[0] = lamb;
    return fext;
  };
  for (auto& element : elements(gridView))
    fes.emplace_back(basis, element, 1000, 0.3,volumeLoad);

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
             .corr_tol  = 1e-16, //everything should converge to the gradient tolerance
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
  if(std::is_same_v<TypeParam , Grids::Yasp>) {
    EXPECT_DOUBLE_EQ(nonLinOp.value(), -1.4809559783564966e+03);
    EXPECT_DOUBLE_EQ(maxDisp, 0.786567027108460048);
  } else if(std::is_same_v<TypeParam , Grids::Alu>) {
    EXPECT_DOUBLE_EQ(nonLinOp.value(), -1.4842107484533601e+03 );
    EXPECT_DOUBLE_EQ(maxDisp, 0.78426066482258983);
  } else if(std::is_same_v<TypeParam , Grids::Iga>) {
    EXPECT_DOUBLE_EQ(nonLinOp.value(), -8.1142552237939071e+02);
    EXPECT_DOUBLE_EQ(maxDisp, 0.615624125459537153);
  }

  nonLinOp.template update<1>();
  EXPECT_TRUE(controlInfo.sucess);
  EXPECT_LE(nonLinOp.derivative().norm(), gradTol);
}