//
// Created by Alex on 21.07.2021.
//

#include "../../config.h"

#include <dune/alugrid/grid.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/iga/nurbsgrid.hh>

#include "spdlog/spdlog.h"

#include <Eigen/Core>

#include "ikarus/Controlroutines/LoadControl.h"
#include "ikarus/FiniteElements/MicroMangeticsWithVectorPotential.h"
#include "ikarus/Solver/NonLinearSolver/NewtonRaphson.hpp"
#include "ikarus/Solver/NonLinearSolver/TrustRegion.hpp"
#include "ikarus/basis/basishelper.h"
#include "ikarus/utils/Observer/controlVTKWriter.h"
#include "ikarus/utils/Observer/nonLinearSolverLogger.h"
#include <ikarus/Assembler/SimpleAssemblers.h>
#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/LinearAlgebra/NonLinearOperator.h>
#include <ikarus/utils/utils/algorithms.h>
#include <ikarus/utils/functionSanityChecks.h>

constexpr int magnetizationOrder = 1;
constexpr int vectorPotOrder = 1;
constexpr int gridDim = 2;
constexpr int directorDim = 3;
constexpr int directorCorrectionDim = directorDim-1;

struct UnitCircleBoundary : Dune::BoundarySegment<2, 2, double> {
  UnitCircleBoundary(const Dune::FieldVector<double, 2>& a, const Dune::FieldVector<double, 2>& b) : corners{{a, b}} {}
  Dune::FieldVector<double, 2> operator()(const Dune::FieldVector<double, 1>& local) const override {
    Dune::FieldVector<double, 2> result = {0, 0};
    double omega                        = std::acos(corners[0] * corners[1]);
    return std::sin((1 - local[0]) * omega) / sin(omega) * corners[0] + sin(local[0] * omega) / sin(omega) * corners[1];
  }

  std::array<Dune::FieldVector<double, 2>, 2> corners;
};

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  using namespace Ikarus;

//  //  /// ALUGrid Example
//  using Grid = Dune::ALUGrid<gridDim, 2, Dune::simplex, Dune::conforming>;
//  auto grid             = Dune::GmshReader<Grid>::read("../../examples/src/testFiles/circleCoarse.msh", false);
//  grid->globalRefine(1);
  /// IGA Grid Example
  Dune::GridFactory<Dune::ALUGrid<gridDim, 2, Dune::cube, Dune::nonconforming>> gridFactory;
  //  std::array<FieldVector<double, 2>, 4> corners0 = {{{-sqrt(2) / 2, -sqrt(2) / 2}, {sqrt(2) / 2, -sqrt(2) / 2},
  //  {sqrt(2) / 2, sqrt(2) / 2}, {-sqrt(2) / 2, sqrt(2) / 2}}};
  Eigen::Vector2d v(1, 0);
  std::array<Dune::FieldVector<double, 2>, 6> corners0;
  Eigen::Rotation2D<double> R;
  R.angle() = 0.0;
  for (auto& corner : corners0) {
    Eigen::Vector2d a = R * v;
    corner[0]         = a[0];
    corner[1]         = a[1];
    R.angle() += 60.0 / 180.0 * std::numbers::pi;
  }

  gridFactory.insertVertex({0, 0});
  gridFactory.insertVertex({1, 0});
  gridFactory.insertVertex({0, 1});
  gridFactory.insertVertex({1, 1});
//  gridFactory.insertVertex(corners0[2]);
//  gridFactory.insertVertex(corners0[3]);
//  gridFactory.insertVertex(corners0[4]);
//  gridFactory.insertVertex(corners0[5]);

  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {0, 1, 2,3});
//  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 2, 3});
//  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 3, 4});
//  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 4, 5});
//  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 5, 6});
//  gridFactory.insertElement(Dune::GeometryTypes::triangle, {0, 6, 1});

  /// Create boundary segments which map the boundaries onto the unit circle
//  gridFactory.insertBoundarySegment({1, 2}, std::make_shared<UnitCircleBoundary>(corners0[0], corners0[1]));
//  gridFactory.insertBoundarySegment({2, 3}, std::make_shared<UnitCircleBoundary>(corners0[1], corners0[2]));
//  gridFactory.insertBoundarySegment({3, 4}, std::make_shared<UnitCircleBoundary>(corners0[2], corners0[3]));
//  gridFactory.insertBoundarySegment({4, 5}, std::make_shared<UnitCircleBoundary>(corners0[3], corners0[4]));
//  gridFactory.insertBoundarySegment({5, 6}, std::make_shared<UnitCircleBoundary>(corners0[4], corners0[5]));
//  gridFactory.insertBoundarySegment({6, 1}, std::make_shared<UnitCircleBoundary>(corners0[5], corners0[0]));

  auto grid     = gridFactory.createGrid();


  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basisEmbedded = makeBasis(gridView, power<directorDim>(lagrange<magnetizationOrder>(), BlockedInterleaved()));
  auto basisRie = makeBasis(gridView, power<directorCorrectionDim>(lagrange<magnetizationOrder>(), FlatInterleaved()));
  std::cout << "This gridview cotains: " << std::endl;
  std::cout << gridView.size(2) << " vertices" << std::endl;
  std::cout << gridView.size(1) << " edges" << std::endl;
  std::cout << gridView.size(0) << " elements" << std::endl;
  std::cout << basisRie.size() << " Dofs" << std::endl;

  draw(gridView);

  Ikarus::FiniteElements::MagneticMaterial mat({.K=2e4,.ms =1.432e6});
  std::vector<Ikarus::FiniteElements::MicroMagneticsWithVectorPotential<decltype(basisEmbedded),decltype(basisRie)>> fes;
  auto volumeLoad = [](auto& globalCoord, auto& lamb)
  {Eigen::Vector<double,directorDim> fext;
    fext.setZero();
    fext[1] = 2 * lamb;
    fext[0] = lamb;
    return fext;
  };

  for (auto& element : elements(gridView))
    fes.emplace_back(basisEmbedded,basisRie, element, mat,volumeLoad);

  std::vector<bool> dirichletFlags(basisRie.size(), false);



  Dune::Functions::forEachBoundaryDOF(basisRie, [&](auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[1]) < 1e-8) {
      dirichletFlags[localView.index(localIndex)[0]] = true;
    }
  });

  auto denseAssembler = DenseFlatAssembler(basisRie, fes, dirichletFlags);
  auto sparseAssembler = SparseFlatAssembler(basisRie, fes, dirichletFlags);

  using DirectorVector     = Dune::BlockVector<Ikarus::UnitVector<double, directorDim>>;
  DirectorVector mBlocked(basisEmbedded.size());
  for (auto& msingle : mBlocked) {
    msingle.setValue(Eigen::Vector<double,directorDim>::Random());
  }
  auto mEigen     = Ikarus::LinearAlgebra::viewAsFlatEigenVector(mBlocked);

  double lambda = 0.0;

  auto residualFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements<DirectorVector> req;
    req.sols.emplace_back(disp);
    req.parameter.insert({Ikarus::FEParameter::loadfactor, lambdaLocal});
    req.matrixAffordances = Ikarus::MatrixAffordances::stiffness;
    return denseAssembler.getVector(req);
  };

  auto hessianFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements<DirectorVector> req;
    req.sols.emplace_back(disp);
    req.parameter.insert({Ikarus::FEParameter::loadfactor, lambdaLocal});
    req.matrixAffordances = Ikarus::MatrixAffordances::stiffness;
    return sparseAssembler.getMatrix(req);
  };

  auto energyFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto {
    Ikarus::FErequirements<DirectorVector> req;
    req.sols.emplace_back(disp);
    req.parameter.insert({Ikarus::FEParameter::loadfactor, lambdaLocal});
    req.matrixAffordances = Ikarus::MatrixAffordances::stiffness;
    return denseAssembler.getScalar(req);
  };
  std::cout<<mBlocked;
  auto& h = hessianFunction(mBlocked,lambda);
  std::cout<<h<<std::endl;

  auto& g = residualFunction(mBlocked,lambda);
  std::cout<<g<<std::endl;

  auto e= energyFunction(mBlocked,lambda);
  std::cout<<e<<std::endl;

  assert(g.size()==gridView.size(2)*directorCorrectionDim && "The returned gradient has incorrect size");

  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(energyFunction, residualFunction, hessianFunction),
                                            parameter(mBlocked, lambda));

  checkGradient(nonLinOp, std::function([&](DirectorVector & x, const Eigen::VectorXd& d) {
                  auto dispEigen     = Ikarus::LinearAlgebra::viewAsFlatEigenVector(x);
                  for (auto i = 0U; i < x.size(); ++i) {
                    size_t indexStartI = i * x[0].correctionSize;
                    x[i] += d.segment<directorCorrectionDim>(indexStartI);
                  }
                }));

  assert(nonLinOp.value()==e );
  assert(nonLinOp.derivative().isApprox(g) );
  assert(nonLinOp.secondDerivative().isApprox(h) );
//
  auto nr                      = Ikarus::makeTrustRegion(nonLinOp, std::function([&](DirectorVector & x, const Eigen::VectorXd& d) {
                                      auto dispEigen     = Ikarus::LinearAlgebra::viewAsFlatEigenVector(x);
                                      for (auto i = 0U; i < x.size(); ++i) {
                                        size_t indexStartI = i * x[0].correctionSize;
                                        x[i] += d.segment<directorCorrectionDim>(indexStartI);
                                      }
                                    }));
  nr->setup({.verbosity = 1, .maxiter = 100000, .grad_tol = 1e-8, .corr_tol = 1e-8,.useRand=false,.rho_reg=1e6, .Delta0 = 1});

  //  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();


  //  nr.subscribeAll(nonLinearSolverObserver);

  auto lc = Ikarus::LoadControl(nr, 1, {0, 1});

//  lc.subscribeAll(vtkWriter);
  std::cout << "Energy before: " << nonLinOp.value() << std::endl;
  std::cout<<mBlocked;
  lc.run();
  nonLinOp.update<0>();
  std::cout << "Energy after: " << nonLinOp.value() << std::endl;
  auto wGlobalFunc = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double,directorDim>>(basisEmbedded, mBlocked);
  Dune::SubsamplingVTKWriter vtkWriter(gridView, Dune::refinementLevels(1));
  vtkWriter.addVertexData(wGlobalFunc, Dune::VTK::FieldInfo("m", Dune::VTK::FieldInfo::Type::scalar, directorDim));
  vtkWriter.write("Magnet");

  std::cout<<mBlocked;
}