//
// Created by ac136645 on 6/15/2022.
//
#include <config.h>

#include <chrono>
#include <vector>

#include <dune/common/parametertreeparser.hh>
#include <dune/fufem/boundarypatch.hh>
#include <dune/fufem/dunepython.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

#include <Eigen/Eigenvalues>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/finiteElements/mechanics/enhancedAssumedStrains.hh>
#include <ikarus/finiteElements/mechanics/linearElastic.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/solver/linearSolver/linearSolver.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>

using namespace Ikarus;
using namespace Dune::Indices;

int main(int argc, char** argv) {
  auto start = std::chrono::high_resolution_clock::now();

  Dune::MPIHelper::instance(argc, argv);
  constexpr int gridDim     = 3;
  double lambdaLoad         = 1;
  constexpr int basis_order = 1;

  /// read in parameters
  Dune::ParameterTree parameterSet;
  Dune::ParameterTreeParser::readINITree(argv[1], parameterSet);

  const Dune::ParameterTree& gridParameters     = parameterSet.sub("GridParameters");
  const Dune::ParameterTree& controlParameters  = parameterSet.sub("ControlParameters");
  const Dune::ParameterTree& materialParameters = parameterSet.sub("MaterialParameters");
  const Dune::ParameterTree& elementParameters  = parameterSet.sub("ElementParameters");

  const double E                   = materialParameters.get<double>("E");
  const double nu                  = materialParameters.get<double>("nu");
  const auto numberOfEASParameters = elementParameters.get<int>("numberOfEASParameters");
  const int refinement_level       = gridParameters.get<int>("refinement");

  using Grid = Dune::UGGrid<gridDim>;
  auto grid  = Dune::GmshReader<Grid>::read("../../tests/src/testFiles/cook_3d.msh", false);
  grid->globalRefine(refinement_level);
  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = makeBasis(gridView, power<gridDim>(lagrange<basis_order>(), FlatInterleaved()));

  std::cout << "This gridview contains: " << std::endl;
  std::cout << gridView.size(2) << " vertices" << std::endl;
  std::cout << gridView.size(1) << " edges" << std::endl;
  std::cout << gridView.size(0) << " elements" << std::endl;
  std::cout << basis.size() << " Dofs" << std::endl;

  //  draw(gridView);

  /// fix left-hand side x-dir
  std::vector<bool> dirichletFlags(basis.size(), false);
  forEachBoundaryDOF(subspaceBasis(basis, 0), [&](auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[0]) < 1e-8) dirichletFlags[localView.index(localIndex)[0]] = true;
  });

  /// fix left-hand side y-dir
  forEachBoundaryDOF(subspaceBasis(basis, 1), [&](auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[0]) < 1e-8) dirichletFlags[localView.index(localIndex)[0]] = true;
  });

  /// fix z-dir line at back z=0
  forEachBoundaryDOF(subspaceBasis(basis, 2), [&](auto&& localIndex, auto&& localView, auto&& intersection) {
    const auto intersectionCenter = intersection.geometry().center();
    if (std::abs(intersectionCenter[2]) < 1e-8) dirichletFlags[localView.index(localIndex)[0]] = true;
  });

  //  std::vector<Ikarus::LinearElastic<decltype(basis)>> fes;
  std::vector<Ikarus::EnhancedAssumedStrains<Ikarus::LinearElastic<decltype(basis)>>> fes;
  //
  /// function for volume load- here: returns zero
  auto volumeLoad = [](auto& globalCoord, auto& lamb) {
    Eigen::Vector3d fext;
    fext.setZero();
    return fext;
  };
  //
  /// neumann boundary load in vertical direction
  auto neumannBoundaryLoad = [&](auto& globalCoord, auto& lamb) {
    Eigen::Vector3d F = Eigen::Vector3d::Zero();
    F[1]              = lamb / 16.0;
    return F;
  };

  /// Python function which could be used to obtain the vertices at the right edge
  std::string lambdaNeumannVertices = std::string("lambda x: ( x[0]>47.9999 )");
  Python::start();
  Python::Reference main = Python::import("__main__");
  Python::run("import math");

  Python::runStream() << std::endl << "import sys" << std::endl << "import os" << std::endl;

  const auto& indexSet = gridView.indexSet();

  /// Flagging the vertices on which neumann load is applied as true
  Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
  auto pythonNeumannVertices = Python::make_function<bool>(Python::evaluate(lambdaNeumannVertices));

  for (auto&& vertex : vertices(gridView)) {
    bool isNeumann                          = pythonNeumannVertices(vertex.geometry().corner(0));
    neumannVertices[indexSet.index(vertex)] = isNeumann;
  }
  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);

  for (auto& element : elements(gridView)) {
    auto localView = basis.localView();
    fes.emplace_back(basis, element, E, nu, &volumeLoad, &neumannBoundary, &neumannBoundaryLoad);
    fes.back().setEASType(numberOfEASParameters);
  }

  auto sparseAssembler = SparseFlatAssembler(basis, fes, dirichletFlags);

  auto KFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::MatrixAffordances::stiffness)
                                     .build();
    return sparseAssembler.getMatrix(req);
  };

  auto residualFunction = [&](auto&& disp, auto&& lambdaLocal) -> auto& {
    Ikarus::FErequirements req = FErequirementsBuilder()
                                     .insertGlobalSolution(Ikarus::FESolutions::displacement, disp)
                                     .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                     .addAffordance(Ikarus::VectorAffordances::forces)
                                     .build();
    return sparseAssembler.getVector(req);
  };

  Eigen::VectorXd D_Glob = Eigen::VectorXd::Zero(basis.size());

  auto startAssembly = std::chrono::high_resolution_clock::now();
  auto nonLinOp
      = Ikarus::NonLinearOperator(linearAlgebraFunctions(residualFunction, KFunction), parameter(D_Glob, lambdaLoad));
  auto stopAssembly     = std::chrono::high_resolution_clock::now();
  auto durationAssembly = duration_cast<std::chrono::milliseconds>(stopAssembly - startAssembly);
  spdlog::info("The assembly took {} milliseconds", durationAssembly.count());
  const auto& K    = nonLinOp.derivative();
  const auto& Fext = nonLinOp.value();

  auto startSolver = std::chrono::high_resolution_clock::now();
  /// solve the linear system
  auto linSolver = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::sd_CholmodSupernodalLLT);
  linSolver.compute(K);
  linSolver.solve(D_Glob, -Fext);
  auto stopSolver     = std::chrono::high_resolution_clock::now();
  auto durationSolver = duration_cast<std::chrono::milliseconds>(stopSolver - startSolver);
  spdlog::info("The solver took {} milliseconds", durationSolver.count());

  /// Postprocess
  auto disp = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 3>>(basis, D_Glob);
  Dune::VTKWriter vtkWriter(gridView, Dune::VTK::conforming);
  vtkWriter.addVertexData(disp, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 3));
  vtkWriter.write("Cook_Membrane_3D");
  auto stop     = std::chrono::high_resolution_clock::now();
  auto duration = duration_cast<std::chrono::milliseconds>(stop - start);
  spdlog::info("The total execution took {} milliseconds", duration.count());
}