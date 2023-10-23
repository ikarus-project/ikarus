//
// Created by Alex on 21.07.2021.
//

#include "config.h"

//#include "src/include/ikarus/utils/algorithms.hh"

#include <dune/alugrid/grid.hh>
#include <dune/common/parametertreeparser.hh>
//#include <dune/fufem/dunepython.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <spdlog/spdlog.h>

#include <Eigen/Core>

#include "ikarus/assembler/simpleAssemblers.hh"
#include "ikarus/controlRoutines/loadControl.hh"
#include <ikarus/finiteElements/mechanics/3dbeam.hh>
#include <ikarus/finiteElements/mechanics/GFE.hh>

#include "ikarus/linearAlgebra/nonLinearOperator.hh"
#include "ikarus/solver/nonLinearSolver/trustRegion.hh"
#include "ikarus/utils/drawing/griddrawer.hh"
#include "ikarus/utils/functionSanityChecks.hh"
#include "ikarus/utils/observer/controlVTKWriter.hh"
#include "ikarus/utils/observer/genericControlObserver.hh"
#include "ikarus/utils/basis.hh"

constexpr int centerLineOrder    = 1;
constexpr int quaternionOrder        = 1;
constexpr int gridDim               = 1;
constexpr int worldDim               = 3;
constexpr int quaternionDim           = 4;
constexpr int quaternionCorrectionDim = quaternionDim - 1;

using CenterLinePositionsVector = Dune::BlockVector<Dune::RealTuple<double, worldDim>>;
using UnitQuaternionVector      = Dune::BlockVector<Dune::UnitVector<double, quaternionDim>>;
using MultiTypeVector = Dune::MultiTypeBlockVector<CenterLinePositionsVector,UnitQuaternionVector >;

int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);
//  Dune::ParameterTree parameterSet;
//  Dune::ParameterTreeParser::readINITree(argv[1], parameterSet);
//
//  const Dune::ParameterTree &gridParameters     = parameterSet.sub("GridParameters");
//  const Dune::ParameterTree &controlParameters  = parameterSet.sub("ControlParameters");
//  const Dune::ParameterTree &materialParameters = parameterSet.sub("MaterialParameters");
//
//  const auto refinement      = gridParameters.get<int>("refinement");
//  const auto innerRadius     = gridParameters.get<double>("innerRadius");
//  const auto mshfilepath     = gridParameters.get<std::string>("mshfilepath");
//  const auto loadSteps       = controlParameters.get<int>("loadSteps");
  const auto loadSteps       = 20;
//  const auto loadFactorRange = controlParameters.get<std::array<double, 2>>("loadFactorRange");
  const auto loadFactorRange = std::array<double, 2>({0,20});
//
//  const auto A  = materialParameters.get<double>("A");
//  const auto K  = materialParameters.get<double>("K");
//  const auto ms = materialParameters.get<double>("ms");

//  Python::start();
//  Python::Reference main = Python::import("__main__");
//  Python::run("import math");

  using namespace Ikarus;
  Ikarus::BeamMaterial mat({.E = 1e8, .nu = 0.5, .A = 1, .I1 = 1, .I2 = 1, .J = 2});


  //  auto isInsidePredicate = [&](auto&& coord) {
  //    if (Dune::power(coord[0], 2) + Dune::power(coord[1], 2) - 1e-4 > Dune::power(innerRadius, 2))
  //      return false;
  //    else
  //      return true;
  //  };

  //  using Grid        = Dune::YaspGrid<gridDim>;
  //  const size_t elex = 60;
  //  const size_t eley = elex / 2;
  //  const size_t elez = 1;
  //  const double Lx   = freeSpaceX;
  //  const double Ly   = freeSpaceY;
  //  const double Lz   = freeSpaceY;
  //
  //  Dune::FieldVector<double, gridDim> bbox;
  //  std::array<int, gridDim> eles{};
  //  if constexpr (gridDim == 2) {
  //    bbox = {Lx, Ly};
  //    eles = {elex, eley};
  //  } else if constexpr (gridDim == 3) {
  //  }
  //
  //  auto grid = std::make_shared<Grid>(bbox, eles);

  Dune::GridFactory<Dune::FoamGrid<1, 3, double>> gridFactory;
  const double h = 1.0;
  const double L = 12.0;
  gridFactory.insertVertex({0, 0,0});
  gridFactory.insertVertex({L,0,0});
  gridFactory.insertElement(Dune::GeometryTypes::line, {0, 1});
  auto grid     = gridFactory.createGrid();
  grid->globalRefine(3);

//  grid->globalRefine(refinement);
  auto gridView = grid->leafGridView();

  //  draw(gridView);
  //  spdlog::info("The domain has a length of {}.", sizedom1);

  using namespace Dune::Functions::BasisFactory;
  //  auto basisEmbedded = makeBasis(gridView, power<directorDim>(lagrange<magnetizationOrder>(),
  //  BlockedInterleaved()));
//  auto basisEmbeddedC
//      = makeBasis(gridView, composite(power<worldDim>(lagrange<centerLineOrder>(), BlockedInterleaved()),
//                                      power<quaternionDim>(lagrange<quaternionOrder>(), BlockedInterleaved()),
//                                      BlockedLexicographic{}));

  //  auto basisRie = makeBasis(gridView, power<quaternionCorrectionDim>(lagrange<magnetizationOrder>(),
  //  FlatInterleaved()));
  auto basis = Ikarus::makeBasis(
      gridView, composite(power<worldDim>(lagrange<centerLineOrder>(), BlockedInterleaved()),
                          power<quaternionDim>(lagrange<quaternionOrder>(), BlockedInterleaved()), BlockedLexicographic{}));
  std::cout << "This gridview contains: " << std::endl;
  std::cout << gridView.size(1) << " edges" << std::endl;
  std::cout << gridView.size(0) << " vertices" << std::endl;
  std::cout << basisRieC.size() << " Dofs" << std::endl;

  //  draw(gridView);

  std::vector<Ikarus::SimoReissnerBeam<decltype(basis)>> fes;
  auto volumeLoad = [](auto &globalCoord, auto &lamb) {
    Eigen::Vector<double, worldDim> fext;
    fext.setZero();
    fext[0] = 0;
    fext[1] = 0;
    fext[2] = -lamb;
    return fext;
  };

  int insideCounter = 0;
  for (auto &element : elements(gridView)) {
    auto geoCoord       = element.geometry().center();
    fes.emplace_back(basisEmbeddedC, basisRieC, element, mat, volumeLoad);
  }

  CenterLinePositionsVector centerLineBlocked(basisEmbeddedC.size({Dune::Indices::_0}));
  for (auto &asingle : centerLineBlocked) {
    asingle.setValue(Eigen::Vector<double, worldDim>::Zero());
  }

  UnitQuaternionVector quatsBlocked(basisEmbeddedC.size({Dune::Indices::_1}));
  for (auto &msingle : quatsBlocked) {
      msingle.setValue(Eigen::Vector<double, quaternionDim>::UnitX());
  }

  MultiTypeVector mAndABlocked( centerLineBlocked,quatsBlocked);

  Ikarus::DirichletValues dirichletValues(basisP->flat());
  // Fix displacement on the whole boundary
  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
//    if (std::abs(intersection.geometry().center()[0]) < 1e-8)
      dirichletFlags[localView.index(localIndex)] = true;
  });

  // Fix displacement on the whole boundary to zero
  Dune::Functions::forEachBoundaryDOF(
      Dune::Functions::subspaceBasis(basisEmbeddedC, Dune::Indices::_0), [&](auto &&globalIndex) {
        mAndABlocked[Dune::Indices::_0][globalIndex[1]].setValue(Eigen::Vector<double, worldDim>::Zero());
      });

//  auto magnetBasis         = Dune::Functions::subspaceBasis(basisRieC, Dune::Indices::_0);
//  auto localView           = magnetBasis.localView();
//  auto seDOFs              = subEntityDOFs(magnetBasis);
//  const auto &gridViewMagn = magnetBasis.gridView();
//  for (auto &&element : elements(gridViewMagn)) {
//    localView.bind(element);
//    for (const auto &intersection : intersections(gridViewMagn, element)) {
//      bool isIntersectionInside = false;
//
//      for (int i = 0; i < intersection.geometry().corners(); ++i) {
//        if (isInsidePredicate(intersection.geometry().corner(i))) {
//          isIntersectionInside = true;
//          break;
//        } else
//          isIntersectionInside = false;
//      }
//
//      if (!isIntersectionInside) {
//        for (auto localIndex : seDOFs.bind(localView, intersection))
//          dirichletFlags[localView.index(localIndex)[0]] = true;
//      }
//    }
//  }

//  auto magnetBasisBlocked   = Dune::Functions::subspaceBasis(basisEmbeddedC, Dune::Indices::_0);
//  auto localView2           = magnetBasisBlocked.localView();
//  auto seDOFs2              = subEntityDOFs(magnetBasisBlocked);
//  const auto &gridViewMagn2 = magnetBasisBlocked.gridView();
//  for (auto &&element : elements(gridViewMagn2)) {
//    localView2.bind(element);
//    for (const auto &intersection : intersections(gridViewMagn2, element))
//      for (auto localIndex : seDOFs2.bind(localView2, intersection)) {
//        bool isIntersectionInside = false;
//
//        for (int i = 0; i < intersection.geometry().corners(); ++i) {
//          if (isInsidePredicate(intersection.geometry().corner(i))) {
//            isIntersectionInside = true;
//            break;
//          } else
//            isIntersectionInside = false;
//        }
//
//        if (!isIntersectionInside) {
//          auto b = mAndABlocked[Dune::Indices::_0][localView2.index(localIndex)[1]].begin();
//          auto e = mAndABlocked[Dune::Indices::_0][localView2.index(localIndex)[1]].end();
//          std::fill(b, e, 0.0);
//        }
//      }
//  }

  auto sparseAssembler = SparseFlatAssembler(fes, dirichletValues);

  double lambda        = 0.0;

  auto residualFunction = [&](auto &&disp, auto &&lambdaLocal) -> auto & {
    Ikarus::FErequirements req = FErequirements<MultiTypeVector>()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacementAndQuaternions, disp)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::VectorAffordances::forces);

    return sparseAssembler.getReducedVector(req);
  };

  auto hessianFunction = [&](auto &&disp, auto &&lambdaLocal) -> auto & {
    Ikarus::FErequirements req = FErequirements<MultiTypeVector>()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacementAndQuaternions, disp)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::MatrixAffordances::stiffness);
    return sparseAssembler.getReducedMatrix(req);
  };

  auto energyFunction = [&](auto &&disp, auto &&lambdaLocal) -> auto {
    Ikarus::FErequirements req = FErequirements<MultiTypeVector>()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacementAndQuaternions, disp)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::ScalarAffordances::mechanicalPotentialEnergy);

    return sparseAssembler.getScalar(req);
  };

  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(energyFunction, residualFunction, hessianFunction),
                                            parameter(mAndABlocked, lambda));
  auto updateFunction = std::function([&](MultiTypeVector &multiTypeVector, const Eigen::VectorXd &d) {
    auto dFull = denseAssembler.createFullVector(d);
    multiTypeVector += dFull;
  });

  checkGradient(nonLinOp, {.draw = true, .writeSlopeStatement = true}, updateFunction);
  checkHessian(nonLinOp, {.draw = true, .writeSlopeStatement = true}, updateFunction);

  auto nr = Ikarus::makeTrustRegion(nonLinOp, updateFunction);
  //  auto nr = Ikarus::makeTrustRegion< decltype(nonLinOp),PreConditioner::DiagonalPreconditioner>(nonLinOp,
  //  updateFunction);
  nr->setup({.verbosity = 1,
             .maxiter   = 100000,
             .grad_tol  = 1e-8,
             .corr_tol  = 1e-16,
             .useRand   = false,
             .rho_reg   = 1e6,
             .Delta0    = 1000});

  auto lc = Ikarus::LoadControl(nr, loadSteps, loadFactorRange);

  auto blockedmidSurfaceBasis2 = Dune::Functions::subspaceBasis(basisEmbeddedC, Dune::Indices::_0);

  auto writerObserver = std::make_shared<Ikarus::ControlSubsamplingVertexVTKWriter,MultiTypeVector>(basis,basis[&](auto i) {
    auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, worldDim>>(
        Dune::Functions::subspaceBasis(basisEmbeddedC, Dune::Indices::_0), mAndABlocked);
    auto quatFunc = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, quaternionDim>>(
        Dune::Functions::subspaceBasis(basisEmbeddedC, Dune::Indices::_1), mAndABlocked);

    Dune::VTKWriter vtkWriter(gridView, Dune::VTK::nonconforming);
    vtkWriter.addVertexData(disp, Dune::VTK::FieldInfo("m", Dune::VTK::FieldInfo::Type::vector, worldDim));
//    vtkWriter.addVertexData(quatFunc, Dune::VTK::FieldInfo("A", Dune::VTK::FieldInfo::Type::vector, vectorPotDim));



    const std::string name3d = "PureBending3D"+std::to_string(gridView.size(0))+"_Order_"+std::to_string(localView.tree().child(Dune::Indices::_0,0).finiteElement().localBasis().order())+"_transverseShearStrainFlag_"+std::to_string(transverseShearStrainFlag)+directorFunction + std::to_string(step);
    std::cout << name3d << std::endl;
    writer2.write(name3d);
  });
  lc.subscribeAll(writerObserver);
  lc.run();
  nonLinOp.update<0>();
  spdlog::info("Energy: {}", nonLinOp.value());
}
