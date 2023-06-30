// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testHelpers.hh"

#include <vector>

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/grid/yaspgrid.hh>

#include "ikarus/finiteElements/mechanics/nonlinearrmshell.hh"
#include "ikarus/finiteElements/physicsHelper.hh"
#include "ikarus/controlRoutines/loadControl.hh"

using Dune::TestSuite;
#include <dune/functions/functionspacebases/interpolate.hh>

#include <Eigen/Core>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/mechanics/nonlinearrmshell.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>
#include <dune/fufem/boundarypatch.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/utils/observer/nonLinearSolverLogger.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <dune/common/tuplevector.hh>
#include <dune/istl/foreach.hh>
#include <dune/iga/nurbsgrid.hh>


constexpr int directorOrder = 1;
constexpr int displacementOrder = 1;

using MidSurfaceVector = Dune::BlockVector<Dune::RealTuple<double, 3>>;
using DirectorVector  = Dune::BlockVector<Dune::UnitVector<double, 3>>;
using MultiTypeVector = Dune::TupleVector<MidSurfaceVector, DirectorVector>;
using MultiTypeVectorRaw = Dune::TupleVector< Dune::BlockVector<Dune::FieldVector<double, 3>>, Dune::BlockVector<Dune::FieldVector<double, 3>>>;

auto squareTest() {
  TestSuite t("SimpleAssemblersTest");
//  using Grid = Dune::YaspGrid<2>;
  using namespace Dune::Indices;
  using namespace Dune;

  constexpr auto dimworld              = 3;
  constexpr auto gridDim              = 2;
  const std::array<int, gridDim> order = {1, 1};

  const std::array<std::vector<double>, gridDim> knotSpans = {{{0, 0, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  const double sizeX = 12;
  const double sizeY = 12;
  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0, 0}, .w = 1}, {.p = {sizeX, 0,0}, .w = 1}},
         {{.p = {0, sizeY}, .w = 1}, {.p = {sizeX, sizeY}, .w = 1}}};

  std::array<int, gridDim> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<gridDim, dimworld>;

  Dune::IGA::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  auto grid               = std::make_shared<Grid>(patchData);

    auto gridView = grid->leafGridView();

    using namespace Dune::Functions::BasisFactory;
    auto basis        = Ikarus::makeBasis(gridView, composite(power<3>(lagrange<displacementOrder >()),power<2>(lagrange<directorOrder>())));
    auto blockedmidSurfaceBasis = subspaceBasis(basis.untouched(),_0);

    MidSurfaceVector mBlocked(basis.untouched().size({Dune::Indices::_0}));
    auto refCoords = [](auto v){ return v;};
//    Functions::interpolate(blockedmidSurfaceBasis, mBlocked, refCoords);
    auto deformationPowerBasis = makeBasis(gridView,power<3>(lagrange<displacementOrder>()));

    std::vector<Eigen::Vector<double, 3> > v;
  for (auto &msingle : mBlocked) {
    msingle.setValue(Eigen::Vector<double, 3>::Zero());
  }


    DirectorVector dBlocked(basis.untouched().size({Dune::Indices::_1}));
    for (auto &dsingle : dBlocked) {
      dsingle.setValue(Eigen::Vector<double, 3>::UnitZ());
    }

    const MultiTypeVector x0(mBlocked, dBlocked);
    MultiTypeVector x =x0;


    auto matParameter = Ikarus::YoungsModulusAndPoissonsRatio({.emodul = 1000, .nu = 0.0});
    const double thickness = 0.1;
    std::vector<Ikarus::NonLinearRMshell<decltype(basis)>> fes;

  auto volumeLoad = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    std::array<Eigen::Vector<double,3>,2> vLoad;
    auto& fext = vLoad[0];
    auto& mext = vLoad[1];
    fext.setZero();
//    fext[0] = lamb;
//    fext[1] = 0.01 * lamb;
    fext[2] = 0 * lamb;
    mext.setZero();
    return vLoad;
  };

  auto neumannBoundaryLoad = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    std::array<Eigen::Vector<double,3>,2> vLoad;
    auto& fext = vLoad[0];
    auto& mext = vLoad[1];
    fext.setZero();
    fext[2] = 0.01 * lamb;
    mext.setZero();
    return vLoad;
  };

  const auto &indexSet = gridView.indexSet();
  // We artificially apply a Neumann load on the complete boundary
  Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
  for (auto &&vertex : vertices(gridView)) {
    bool isNeumann                          = vertex.geometry().corner(0)[0]<1e-8;
    neumannVertices[indexSet.index(vertex)] = isNeumann;
  }

  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);
    for (auto&& ge : elements(gridView))
      fes.emplace_back(basis, ge, matParameter,thickness, x0,volumeLoad, &neumannBoundary, neumannBoundaryLoad);

//    auto basisP = std::make_shared<const decltype(basis)>(basis);
    Ikarus::DirichletValues dirichletValues(basis.flat());
    dirichletValues.fixDOFs([](auto& basis_, auto& dirichletFlags) {
      Dune::Functions::forEachBoundaryDOF(basis_, [&](auto&& localIndex,auto&& localView,auto&& intersection) {
        if(intersection.geometry().center()[0]<1e-8)
          dirichletFlags[localView.index(localIndex)] = true;
      });
    });

    Ikarus::SparseFlatAssembler sparseFlatAssembler(fes, dirichletValues);

    double lambda = 0;
     auto residualFunction = [&](auto &&disp, auto &&lambdaLocal) -> auto & {
      auto req = Ikarus::FErequirements<MultiTypeVector>()
                                       .insertGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector, x)
                                       .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                       .addAffordance(Ikarus::AffordanceCollections::elastoStatics);
      return sparseFlatAssembler.getVector(req);
    };

    auto hessianFunction = [&](auto &&disp, auto &&lambdaLocal) -> auto & {
      auto req = Ikarus::FErequirements<MultiTypeVector>()
                                       .insertGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector, disp)
                                       .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                       .addAffordance(Ikarus::AffordanceCollections::elastoStatics);
      return sparseFlatAssembler.getMatrix(req);
    };

    auto energyFunction = [&](auto &&disp, auto &&lambdaLocal) -> auto {
      auto req = Ikarus::FErequirements<MultiTypeVector>()
                                       .insertGlobalSolution(Ikarus::FESolutions::midSurfaceAndDirector, disp)
                                       .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                       .addAffordance(Ikarus::AffordanceCollections::elastoStatics);


      return sparseFlatAssembler.getScalar(req);
    };

    auto nonLinOp = Ikarus::NonLinearOperator(Ikarus::functions(energyFunction, residualFunction, hessianFunction),
                                              Ikarus::parameter(x, lambda));

    std::cout<<nonLinOp.value()<<std::endl;
    std::cout<<nonLinOp.derivative()<<std::endl;
    std::cout<<nonLinOp.secondDerivative()<<std::endl;

  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::sd_UmfPackLU);
#if solverType == 0
  auto nr = Ikarus::makeNewtonRaphson(nonLinOp.subOperator<1, 2>(), std::move(linSolver));
#endif
#if solverType == 1
  auto nr = Ikarus::makeTrustRegion(nonLinOp);
  nr->setup({.verbosity = 1,
             .maxiter   = 30,
             .grad_tol  = 1e-8,
             .corr_tol  = 1e-8,
             .useRand   = false,
             .rho_reg   = 1e6,
             .Delta0    = 1});
#endif

  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
//  Eigen::VectorXd x0Flat;
//  // define a lambda that copies the entries of a blocked vector into a flat one
//  auto copyToFlat = [&](auto&& blocked, auto&& flat)
//  {
//    Dune::flatVectorForEach(blocked, [&](auto&& entry, auto&& offset)
//    {
//      flat[offset] = entry;
//    });
//  };
//  copyToFlat(x0,x0Flat);


  MultiTypeVectorRaw dispVec;
  auto lambdaForWriting = [&](auto&&vtkWriterL, auto&& basis, auto&& sol,auto && prefixString, int step)
  {
//    MultiTypeVectorRaw dispVec;
    dispVec[_0].resize(x0[_0].size());
    dispVec[_1].resize(x0[_1].size());
    for (int i = 0; i < x0[_0].size(); ++i) {
      for (int j = 0; j < sol[_0][i].getValue().size(); ++j) {

      dispVec[_0][i][j] = sol[_0][i].getValue()[j] -x0[_0][i].getValue()[j];
      }
    }
    for (int i = 0; i < x0[_1].size(); ++i) {
      for (int j = 0; j < sol[_1][i].getValue().size(); ++j) {
        dispVec[_1][i][j] = sol[_1][i].getValue()[j] - x0[_1][i].getValue()[j];
      }
    }

    auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(subspaceBasis(basis,_0),
                                                                                               dispVec);
    Dune::VTK::FieldInfo fieldInfoM{"MidSurfaceDisplacements", Dune::VTK::FieldInfo::Type::vector, 3};
    vtkWriterL.addVertexData(disp, fieldInfoM);

    auto director = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(subspaceBasis(basis,_1),
                                                                                                   sol);
    Dune::VTK::FieldInfo fieldInfoD{"Director", Dune::VTK::FieldInfo::Type::vector, 3};
    vtkWriterL.addVertexData(director, fieldInfoD);
    vtkWriterL.write(prefixString + std::to_string(step++));
  };
  using VtkWriter = ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>,MultiTypeVector>;
  auto vtkWriter = std::make_shared<VtkWriter>(basis.flat(), x,lambdaForWriting, 2);
  vtkWriter->setFileNamePrefix("reisserMindlinshellEX");
//  vtkWriter->setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 3);
  nr->subscribeAll(nonLinearSolverObserver);

  auto lc = Ikarus::LoadControl(nr, 20, {0, 1});

  lc.subscribeAll(vtkWriter);
  std::cout << "Energy before: " << nonLinOp.value() << std::endl;
  lc.run();
  nonLinOp.update<0>();
  std::cout << "Energy after: " << nonLinOp.value() << std::endl;

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest( squareTest());
  return t.exit();
}
