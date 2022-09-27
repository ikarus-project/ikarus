/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#include <config.h>

#include <dune/common/parametertreeparser.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/fufem/dunepython.hh>

#include "spdlog/spdlog.h"

#include <Eigen/Core>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlRoutines/loadControl.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/mechanics/nonLinearElasticityFE.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/manifolds/realTuple.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/solver/nonLinearSolver/trustRegion.hh>
#include <ikarus/utils/algorithms.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <ikarus/utils/observer/nonLinearSolverLogger.hh>

using namespace Ikarus;
using namespace Dune::Indices;

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  constexpr int gridDim = 2;
  const double L2  = 10.0;

  /// read in parameters
  Dune::ParameterTree parameterSet;
  Dune::ParameterTreeParser::readINITree(argv[1], parameterSet);

  const Dune::ParameterTree& gridParameters     = parameterSet.sub("GridParameters");
  const Dune::ParameterTree& materialParameters = parameterSet.sub("MaterialParameters");
  const Dune::ParameterTree& elementParameters  = parameterSet.sub("ElementParameters");

  const auto E                      = materialParameters.get<double>("E");
  const auto nu                     = materialParameters.get<double>("nu");
  const auto refinement_level       = gridParameters.get<int>("refinement");
  const auto r                      = gridParameters.get<double>("r");      // aspect ratio of the rubber block
  const auto ele_x                  = elementParameters.get<int>("ele_x");
  const auto ele_y                  = elementParameters.get<int>("ele_y");
  const auto numberOfEASParameters  = elementParameters.get<int>("numberOfEASParameters");
  const auto L1                     = r * L2;

  if ((ele_x<=0) or (ele_y<=0))
      DUNE_THROW(Dune::MathError,"Number of elements should be greater than zero");

  /// YaspGrid Example
  using Grid        = Dune::YaspGrid<gridDim>;

  Dune::FieldVector<double, 2> bbox = {L1, L2};
  std::array<int, 2> elems           = {ele_x, ele_y};
  auto grid                         = std::make_shared<Grid>(bbox, elems);
  grid->globalRefine(refinement_level);

  auto gridView = grid->leafGridView();

  const auto& indexSet = gridView.indexSet();

  Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);

  std::string lambdaNeumannVertices = std::string("lambda x: ( x[1] > 9.9999 )");
  Python::start();
  Python::Reference main = Python::import("__main__");
  Python::run("import math");

  Python::runStream() << std::endl << "import sys" << std::endl << "import os" << std::endl;

  auto pythonNeumannVertices = Python::make_function<bool>(Python::evaluate(lambdaNeumannVertices));

  for (auto&& vertex : vertices(gridView)) {
      bool isNeumann                          = pythonNeumannVertices(vertex.geometry().corner(0));
      neumannVertices[indexSet.index(vertex)] = isNeumann;
  }

  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);

  auto neumannBoundaryLoad = [](auto& globalCoord, auto& lamb) {
      Eigen::Vector2d fext;
      fext.setZero();
      fext[1] = -lamb; // compression
      fext[0] = 0;
      return fext;
  };

  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, power<gridDim>(lagrange<1>(), FlatInterleaved()));
  std::cout << "This gridview contains: " << std::endl;
  std::cout << gridView.size(2) << " vertices" << std::endl;
  std::cout << gridView.size(1) << " edges" << std::endl;
  std::cout << gridView.size(0) << " elements" << std::endl;
  std::cout << basis.size() << " Dofs" << std::endl;

//  draw(gridView);

  auto localView = basis.localView();
  std::vector<Ikarus::NonLinearElasticityFE<decltype(basis)>> fes;
  auto volumeLoad = [](auto& globalCoord, auto& lamb) {
    Eigen::Vector2d fext;
    fext.setZero();
    return fext;
  };

  for (auto& element : elements(gridView))
    fes.emplace_back(basis, element, E, nu, &neumannBoundary, neumannBoundaryLoad, volumeLoad);

  std::vector<bool> dirichletFlags(basis.size(), false);

  Dune::Functions::forEachBoundaryDOF(subspaceBasis(basis, _1), [&](auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[1]) < 1e-8) dirichletFlags[localView.index(localIndex)[0]] = true;
  });

  dirichletFlags[ele_x] = true; //single node fix in x-direction

  auto sparseAssembler = SparseFlatAssembler(basis, fes, dirichletFlags);

  Eigen::VectorXd d;
  d.setZero(basis.size());
  double lambda = 0.0;

  auto residualFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::VectorAffordances::forces)
                                     .build();
    return sparseAssembler.getVector(req);
  };

  auto KFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::MatrixAffordances::stiffness)
                                     .build();
    return sparseAssembler.getMatrix(req);
  };

  auto energyFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::ScalarAffordances::mechanicalPotentialEnergy)
                                     .build();
    return sparseAssembler.getScalar(req);
  };

  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(energyFunction, residualFunction, KFunction),
                                            parameter(d, lambda));

  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::sd_UmfPackLU);

  auto nr = Ikarus::makeNewtonRaphson(nonLinOp.subOperator<1, 2>(), std::move(linSolver));
//  auto nr = Ikarus::makeTrustRegion(nonLinOp);
//  nr->setup({.verbosity = 1,
//             .maxiter   = 30,
//             .grad_tol  = 1e-9,
//             .corr_tol  = 1e-9,
//             .useRand   = false,
//             .rho_reg   = 1e6,
//             .Delta0    = 1});

  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();

  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<decltype(basis)>>(basis, d, 0);
  vtkWriter->setFileNamePrefix("BifurcationRubberBlock");
  vtkWriter->setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 2);
  nr->subscribeAll(nonLinearSolverObserver);

  auto lc = Ikarus::LoadControl(nr, 200, {0, 2000});

  lc.subscribeAll(vtkWriter);
  std::cout << "Energy before: " << nonLinOp.value() << std::endl;
  lc.run();
  nonLinOp.update<0>();
  std::cout << "Energy after: " << nonLinOp.value() << std::endl;
}