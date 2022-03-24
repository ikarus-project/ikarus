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
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/iga/nurbsgrid.hh>

#include "spdlog/spdlog.h"

#include <Eigen/Core>

#include "ikarus/Controlroutines/LoadControl.h"
#include "ikarus/FiniteElements/Micromagnetics/MicroMangeticsWithVectorPotential.h"
#include "ikarus/Solver/NonLinearSolver/TrustRegion.hpp"
#include "ikarus/utils/Observer/controlVTKWriter.h"
#include "ikarus/utils/Observer/genericControlObserver.h"
#include <ikarus/Assembler/SimpleAssemblers.h>
#include "ikarus/utils/drawing/griddrawer.h"
#include <ikarus/LinearAlgebra/NonLinearOperator.h>
#include <ikarus/utils/functionSanityChecks.h>
#include <ikarus/utils/utils/algorithms.h>

constexpr int magnetizationOrder    = 1;
constexpr int vectorPotOrder        = 1;
constexpr int gridDim               = 2;
constexpr int directorDim           = 3;
constexpr int vectorPotDim          = gridDim == 2 ? 1 : 3;
constexpr int directorCorrectionDim = directorDim - 1;

using DirectorVector  = Dune::BlockVector<Ikarus::UnitVector<double, directorDim>>;
using VectorPotVector = Dune::BlockVector<Ikarus::RealTuple<double, vectorPotDim>>;
using MultiTypeVector = Dune::MultiTypeBlockVector<DirectorVector, VectorPotVector>;

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  using namespace Ikarus;
  Ikarus::FiniteElements::MagneticMaterial mat({.A = 1.0e-11, .K = 2e4, .ms = 1.432e6});
  const double lx              = sqrt(2 * mat.A / (mat.mu0 * mat.ms * mat.ms));
  const double lengthUnit      = 1e-9;
  const double sizedom1InMeter = 60 * lengthUnit;
  const double sizedom1        = sizedom1InMeter / lx;
  const double sizedom2        = sizedom1 / 2;

    using Grid        = Dune::YaspGrid<gridDim>;
    const double Lx   = sizedom1;
    const double Ly   = sizedom2;
    const double Lz   = sizedom2;
    const size_t elex = 5;
    const size_t eley = elex/2;
    const size_t elez = 1;

    Dune::FieldVector<double, gridDim> bbox;
    std::array<int, gridDim> eles;
    if constexpr (gridDim == 2) {
      bbox = {Lx, Ly};
      eles = {elex, eley};
    } else if constexpr (gridDim == 3) {
    }

    auto grid = std::make_shared<Grid>(bbox, eles);

//  using Grid = Dune::ALUGrid<gridDim, 2, Dune::simplex, Dune::conforming>;
//  auto grid  = Dune::GmshReader<Grid>::read("../../examples/src/testFiles/circle.msh", false);

    grid->globalRefine(0);
  auto gridView = grid->leafGridView();

//  draw(gridView);
  spdlog::info("The exchange length is {}.", lx);
  spdlog::info("The domain has a length of {}.", sizedom1);

  using namespace Dune::Functions::BasisFactory;
  //  auto basisEmbedded = makeBasis(gridView, power<directorDim>(lagrange<magnetizationOrder>(),
  //  BlockedInterleaved()));
  auto basisEmbeddedC
      = makeBasis(gridView, composite(power<directorDim>(lagrange<magnetizationOrder>(), BlockedInterleaved()),
                                      power<vectorPotDim>(lagrange<vectorPotOrder>(), BlockedInterleaved()),
                                      BlockedLexicographic{}));

  //  auto basisRie = makeBasis(gridView, power<directorCorrectionDim>(lagrange<magnetizationOrder>(),
  //  FlatInterleaved()));
  auto basisRieC = makeBasis(
      gridView, composite(power<directorCorrectionDim>(lagrange<magnetizationOrder>(), FlatInterleaved()),
                          power<vectorPotDim>(lagrange<vectorPotOrder>(), FlatInterleaved()), FlatLexicographic{}));
  std::cout << "This gridview contains: " << std::endl;
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

  DirectorVector mBlocked(basisEmbeddedC.size({Dune::Indices::_0}));
  for (auto& msingle : mBlocked) {
    msingle.setValue(0.1*Eigen::Vector<double, directorDim>::Random()+Eigen::Vector<double, directorDim>::UnitZ());
  }

  VectorPotVector aBlocked(basisEmbeddedC.size({Dune::Indices::_1}));
  for (auto& asingle : aBlocked) {
    asingle.setValue(Eigen::Vector<double, vectorPotDim>::Zero());
  }

  MultiTypeVector mAndABlocked(mBlocked, aBlocked);

  std::vector<bool> dirichletFlags(basisRieC.size(), false);
  std::cout << "dirichletFlags.size()" << dirichletFlags.size() << std::endl;
  // Fix vector potential on the whole boundary
  Dune::Functions::forEachBoundaryDOF(Dune::Functions::subspaceBasis(basisRieC, Dune::Indices::_1),
                                      [&](auto&& globalIndex) { dirichletFlags[globalIndex[0]] = true; });

  auto denseAssembler  = DenseFlatAssembler(basisRieC, fes, dirichletFlags);
  auto sparseAssembler = SparseFlatAssembler(basisRieC, fes, dirichletFlags);
  double lambda = 1.0;


  auto residualFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements<MultiTypeVector> req;
    req.sols.emplace_back(disp);
    req.parameter.insert({Ikarus::FEParameter::loadfactor, lambdaLocal});
    return denseAssembler.getReducedVector(req);
  };

  auto hessianFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements<MultiTypeVector> req;
    req.sols.emplace_back(disp);
    req.parameter.insert({Ikarus::FEParameter::loadfactor, lambdaLocal});
    return sparseAssembler.getReducedMatrix(req);
  };

  auto energyFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto {
    Ikarus::FErequirements<MultiTypeVector> req;
    req.sols.emplace_back(disp);
    req.parameter.insert({Ikarus::FEParameter::loadfactor, lambdaLocal});
    return denseAssembler.getScalar(req);
  };

  //  auto& h = hessianFunction(mAndABlocked, lambda);
  //    std::cout <<"hbig"<< h << std::endl;
  //
  //  auto& g = residualFunction(mAndABlocked, lambda);
  //    std::cout <<"g"<< g << std::endl;
  //
  //  auto e = energyFunction(mAndABlocked, lambda);
  //    std::cout <<"e"<< e << std::endl;

  //  assert(g.size() == gridView.size(2) * directorCorrectionDim +gridView.size(2) * vectorPotDim -
  //  std::ranges::count(dirichletFlags, true)
  //         && "The returned gradient has incorrect size");

  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(energyFunction, residualFunction, hessianFunction),
                                            parameter(mAndABlocked, lambda));
  std::cout << "CP4" << std::endl;
  auto updateFunction = std::function([&](MultiTypeVector& multiTypeVector, const Eigen::VectorXd& d) {
    auto dFull = denseAssembler.createFullVector(d);
    multiTypeVector += dFull;
  });

  checkGradient(nonLinOp, true, updateFunction);
  checkHessian(nonLinOp, true, updateFunction);


  auto nr = Ikarus::makeTrustRegion(nonLinOp, updateFunction);
//  auto nr = Ikarus::makeTrustRegion< decltype(nonLinOp),PreConditioner::DiagonalPreconditioner>(nonLinOp, updateFunction);
  nr->setup({.verbosity = 1,
             .maxiter   = 100000,
             .grad_tol  = 1e-8,
             .corr_tol  = 1e-8,
             .useRand   = false,
             .rho_reg   = 1e6,
             .Delta0    = 1000});

  auto lc = Ikarus::LoadControl(nr, 1, {0, 100000});

  auto scalarMagnBasis = makeBasis(gridView, lagrange<magnetizationOrder>());
  auto localViewScalarMagnBasis = scalarMagnBasis.localView();

  std::vector<double> gradMNodalRes(scalarMagnBasis.size());

  auto writerObserver = std::make_shared<Ikarus::GenericControlObserver>(ControlMessages::STEP_ENDED, [&](auto i) {
    auto mGlobalFunc = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, directorDim>>(
        Dune::Functions::subspaceBasis(basisEmbeddedC, Dune::Indices::_0), mAndABlocked);
    auto AGlobalFunc = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, directorDim>>(
        Dune::Functions::subspaceBasis(basisEmbeddedC, Dune::Indices::_1), mAndABlocked);

    Dune::VTKWriter vtkWriter(gridView, Dune::VTK::nonconforming);
    vtkWriter.addVertexData(mGlobalFunc, Dune::VTK::FieldInfo("m", Dune::VTK::FieldInfo::Type::vector, directorDim));
    vtkWriter.addVertexData(AGlobalFunc, Dune::VTK::FieldInfo("A", Dune::VTK::FieldInfo::Type::vector, directorDim));

    Ikarus::ResultRequirements<Ikarus::FErequirements<MultiTypeVector>> resultRequirements;
    resultRequirements.req.sols.emplace_back(mAndABlocked);
    resultRequirements.req.parameter.insert({Ikarus::FEParameter::loadfactor, lambda});
    resultRequirements.resType = ResultType::gradientNormOfMagnetization;
    auto localmFunction  = localFunction(mGlobalFunc);

    std::vector<double> gradMEvaluatedAtNodes(basisRieC.size());
    auto ele = gridView.template begin<0>();
    Eigen::VectorXd result;
    for (auto & fe : fes) {
      localViewScalarMagnBasis.bind(*ele);
      for (int c = 0; c < ele->geometry().corners(); ++c) {
        auto coord = toEigenVector(ele->geometry().corner(c));

         fe.calculateAt(resultRequirements,coord,result);
         gradMNodalRes[localViewScalarMagnBasis.index(localViewScalarMagnBasis.tree().localIndex(c))[0]] =result[0];
      }
      ++ele;
    }

    auto gradmGlobalFunc = Dune::Functions::makeDiscreteGlobalBasisFunction<double>(
        scalarMagnBasis, gradMNodalRes);
    vtkWriter.addVertexData(gradmGlobalFunc, Dune::VTK::FieldInfo("gradMNorm", Dune::VTK::FieldInfo::Type::scalar, 1));
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