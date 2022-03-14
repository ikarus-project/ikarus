//
// Created by Alex on 21.07.2021.
//

#include "../../config.h"

#include <dune/alugrid/grid.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
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
#include "ikarus/utils/Observer/genericControlObserver.h"
#include <ikarus/Assembler/SimpleAssemblers.h>
#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/LinearAlgebra/NonLinearOperator.h>
#include <ikarus/utils/functionSanityChecks.h>
#include <ikarus/utils/utils/algorithms.h>

constexpr int magnetizationOrder    = 1;
constexpr int vectorPotOrder        = 1;
constexpr int gridDim               = 2;
constexpr int directorDim           = 3;
constexpr int vectorPotDim          = gridDim == 2 ? 1 : 3;
constexpr int directorCorrectionDim = directorDim - 1;

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  using namespace Ikarus;
  Ikarus::FiniteElements::MagneticMaterial mat({.A = 1.0e-11, .K = 2e4, .ms = 1.432e6});
  const double lx              = sqrt(2 * mat.A / (mat.mu0 * mat.ms * mat.ms));
  const double lengthUnit      = 1e-9;
  const double sizedom1InMeter = 60 * lengthUnit;
  const double sizedom1        = sizedom1InMeter / lx;
  const double sizedom2        = sizedom1;

  using Grid        = Dune::YaspGrid<gridDim>;
  const double Lx   = sizedom1;
  const double Ly   = sizedom2;
  const size_t elex = 1;
  const size_t eley = 1;

  Dune::FieldVector<double, 2> bbox = {Lx, Ly};
  std::array<int, 2> eles           = {elex, eley};
  auto grid                         = std::make_shared<Grid>(bbox, eles);

  grid->globalRefine(5);
  auto gridView = grid->leafGridView();

  spdlog::info("The exchange length is {}.", lx);
  spdlog::info("The domain has a length of {}.", sizedom1);

  using namespace Dune::Functions::BasisFactory;
//  auto basisEmbedded = makeBasis(gridView, power<directorDim>(lagrange<magnetizationOrder>(), BlockedInterleaved()));
  auto basisEmbeddedC
      = makeBasis(gridView, composite(power<directorDim>(lagrange<magnetizationOrder>(), BlockedInterleaved()),
                                      power<vectorPotDim>(lagrange<vectorPotOrder>(), FlatInterleaved()),FlatLexicographic{}));
  auto localViewC = basisEmbeddedC.localView();
  for (auto& element : elements(gridView))
  {
    using namespace Dune::Indices;
    localViewC.bind(element);
    std::cout<<"_0,0,0: "<<localViewC.index(localViewC.tree().child(_0,0).localIndex(0))<<std::endl;
    std::cout<<"_0,1,0: "<<localViewC.index(localViewC.tree().child(_0,1).localIndex(0))<<std::endl;
    std::cout<<"_1,0,0: "<<localViewC.index(localViewC.tree().child(_1,0).localIndex(0))<<std::endl;
    std::cout<<"_1,0,0: "<<localViewC.index(localViewC.tree().child(_1,1).localIndex(0))<<std::endl;
  }
//  auto basisRie = makeBasis(gridView, power<directorCorrectionDim>(lagrange<magnetizationOrder>(), FlatInterleaved()));
  auto basisRieC
      = makeBasis(gridView, composite(power<directorDim>(lagrange<magnetizationOrder>(), FlatInterleaved()),
                                      power<vectorPotDim>(lagrange<vectorPotOrder>(), FlatInterleaved()),FlatLexicographic{}));
  std::cout << "This gridview cotains: " << std::endl;
  std::cout << gridView.size(2) << " vertices" << std::endl;
  std::cout << gridView.size(1) << " edges" << std::endl;
  std::cout << gridView.size(0) << " elements" << std::endl;
  std::cout << basisRieC.size() << " Dofs" << std::endl;

  //  draw(gridView);

  std::vector<Ikarus::FiniteElements::MicroMagneticsWithVectorPotential<decltype(basisEmbeddedC), decltype(basisRieC)>>
      fes;
  auto volumeLoad = [](auto& globalCoord, auto& lamb) {
    Eigen::Vector<double, directorDim> fext;
    fext.setZero();
    fext[0] = lamb;
    fext[1] = lamb;
    return fext;
  };

  for (auto& element : elements(gridView))
    fes.emplace_back(basisEmbeddedC, basisRieC, element, mat, volumeLoad);

  using DirectorVector = Dune::BlockVector<Ikarus::UnitVector<double, directorDim>>;
  DirectorVector mBlocked(basisEmbedded.size());
  for (auto& msingle : mBlocked) {
    msingle.setValue(Eigen::Vector<double, directorDim>::UnitX());
  }

  std::vector<bool> dirichletFlagsEmbedded(basisEmbeddedC.size(), false);
  std::vector<bool> dirichletFlags(basisRieC.size(), false);

  Dune::Functions::forEachBoundaryDOF(subSpaceBasis(basisEmbedded,Dune::Indices::_1), [&](auto&& globalIndex) {
    dirichletFlagsEmbedded[globalIndex[0]] = true;
    if (intersection.geometry().center()[1] < 1e-8)
      mBlocked[localView.index(localIndex)[0]].setValue(Eigen::Vector<double, directorDim>::UnitX());
    else if (intersection.geometry().center()[1] > Ly - 1e-8)
      mBlocked[localView.index(localIndex)[0]].setValue(Eigen::Vector<double, directorDim>::UnitX());
    else if (intersection.geometry().center()[0] > Lx - 1e-8)
      mBlocked[localView.index(localIndex)[0]].setValue(Eigen::Vector<double, directorDim>::UnitZ());
    else if (intersection.geometry().center()[0] < 1e-8)
      mBlocked[localView.index(localIndex)[0]].setValue(-Eigen::Vector<double, directorDim>::UnitZ());
  });

  Dune::Functions::forEachBoundaryDOF(basisRie, [&](auto&& localIndex, auto&& localView, auto&& intersection) {
    dirichletFlags[localView.index(localIndex)[0]] = true;
  });

  auto denseAssembler  = DenseFlatAssembler(basisRie, fes, dirichletFlags);
  auto sparseAssembler = SparseFlatAssembler(basisRie, fes, dirichletFlags);

  double lambda = 0.0;

  auto residualFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements<DirectorVector> req;
    req.sols.emplace_back(disp);
    req.parameter.insert({Ikarus::FEParameter::loadfactor, lambdaLocal});
    req.matrixAffordances = Ikarus::MatrixAffordances::stiffness;
    return denseAssembler.getReducedVector(req);
  };

  auto hessianFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements<DirectorVector> req;
    req.sols.emplace_back(disp);
    req.parameter.insert({Ikarus::FEParameter::loadfactor, lambdaLocal});
    req.matrixAffordances = Ikarus::MatrixAffordances::stiffness;
    return sparseAssembler.getReducedMatrix(req);
  };

  auto energyFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto {
    Ikarus::FErequirements<DirectorVector> req;
    req.sols.emplace_back(disp);
    req.parameter.insert({Ikarus::FEParameter::loadfactor, lambdaLocal});
    req.matrixAffordances = Ikarus::MatrixAffordances::stiffness;
    return denseAssembler.getScalar(req);
  };

  auto& h = hessianFunction(mBlocked, lambda);
  //  std::cout << h << std::endl;

  auto& g = residualFunction(mBlocked, lambda);
  //  std::cout << g << std::endl;

  auto e = energyFunction(mBlocked, lambda);
  //  std::cout << e << std::endl;

  assert(g.size() == gridView.size(2) * directorCorrectionDim - std::ranges::count(dirichletFlags, true)
         && "The returned gradient has incorrect size");

  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(energyFunction, residualFunction, hessianFunction),
                                            parameter(mBlocked, lambda));

  auto updateFunction = std::function([&](DirectorVector& x, const Eigen::VectorXd& d) {
    auto dFull = denseAssembler.createFullVector(d);
    x += dFull;
  });

    checkGradient(nonLinOp, true, updateFunction);
    checkHessian(nonLinOp, true, updateFunction);

  if (not Dune::FloatCmp::eq(nonLinOp.value(), e)) throw std::logic_error("Dune::FloatCmp::eq(nonLinOp.value(), e)");
  if (not nonLinOp.derivative().isApprox(g)) throw std::logic_error("nonLinOp.derivative().isApprox(g)");
  if (not nonLinOp.secondDerivative().isApprox(h)) throw std::logic_error("nonLinOp.secondDerivative().isApprox(h)");

  auto nr = Ikarus::makeTrustRegion(nonLinOp, updateFunction);
  nr->setup({.verbosity = 1,
             .maxiter   = 100000,
             .grad_tol  = 1e-8,
             .corr_tol  = 1e-8,
             .useRand   = false,
             .rho_reg   = 1e6,
             .Delta0    = 1});

  //  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();

  //  nr.subscribeAll(nonLinearSolverObserver);

  auto lc = Ikarus::LoadControl(nr, 100, {0, 100000});

  //  lc.subscribeAll(vtkWriter);
  std::cout << "Energy before: " << nonLinOp.value() << std::endl;
  std::cout << "Energy after: " << nonLinOp.value() << std::endl;


  auto writerObserver = std::make_shared<Ikarus::GenericControlObserver>(ControlMessages::STEP_ENDED,[&](auto i){
    auto wGlobalFunc = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, directorDim>>(
        basisEmbedded, mBlocked);
    Dune::VTKWriter vtkWriter(gridView);
    vtkWriter.addVertexData(wGlobalFunc, Dune::VTK::FieldInfo("m", Dune::VTK::FieldInfo::Type::vector, directorDim));
    vtkWriter.write(std::string("Magnet") + std::to_string(i));
  });
  lc.subscribeAll(writerObserver);
  lc.run();
  nonLinOp.update<0>();


  for (auto& mS : mBlocked)
    if (not Dune::FloatCmp::eq(mS.getValue().norm(), 1.0))
      std::cout << "wrong director found " << mS.getValue().transpose() << std::endl;
  //  std::cout << mBlocked;
}