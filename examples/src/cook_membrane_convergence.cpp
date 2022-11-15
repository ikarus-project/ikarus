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

int main(int argc, char **argv) {
  //  auto start = std::chrono::high_resolution_clock::now();
  Dune::MPIHelper::instance(argc, argv);
  constexpr int gridDim     = 2;
  double lambdaLoad         = 1;
  constexpr int basis_order = 1;

  const double E  = 1.0;
  const double nu = 1.0 / 3.0;

  using Grid = Dune::UGGrid<gridDim>;

  Eigen::Vector<int, 4> easSet;
  easSet << 0, 4, 5, 7;

  std::vector<double> dofsVec;
  std::vector<int> timeVec;
  std::vector<double> dispVec;
  std::vector<std::string> legends;
  /// Draw convergence plots
  using namespace matplot;
  auto f  = figure(true);
  auto ax = gca();
  ax->y_axis().label("Displacement at the top-right tip");
  ax->x_axis().label("Dofs");

  auto f2             = figure(true);
  auto axesSecondPlot = gca();
  axesSecondPlot->y_axis().label("Displacement at the top-right tip");
  axesSecondPlot->x_axis().label("Assembly time in ms");

  for (size_t nep = 0; nep < easSet.size(); ++nep) {
    dofsVec.clear();
    dispVec.clear();
    timeVec.clear();
    auto grid = Dune::GmshReader<Grid>::read("../../tests/src/testFiles/cook.msh", false);
    for (size_t ref = 0; ref < 9; ++ref) {
      auto start                 = std::chrono::high_resolution_clock::now();
      auto gridView              = grid->leafGridView();
      auto numberOfEASParameters = easSet(nep);

      using namespace Dune::Functions::BasisFactory;
      auto basis = makeBasis(gridView, power<gridDim>(lagrange<basis_order>(), FlatInterleaved()));

      /// clamp left-hand side
      std::vector<bool> dirichletFlags(basis.size(), false);
      forEachBoundaryDOF(basis, [&](auto &&localIndex, auto &&localView, auto &&intersection) {
        if (std::abs(intersection.geometry().center()[0]) < 1e-8) dirichletFlags[localView.index(localIndex)[0]] = true;
      });

      std::vector<Ikarus::EnhancedAssumedStrains<Ikarus::LinearElastic<decltype(basis)>>> fes;

      /// function for volume load- here: returns zero
      auto volumeLoad = [](auto &globalCoord, auto &lamb) {
        Eigen::Vector2d fext;
        fext.setZero();
        fext[1] = 2 * lamb * 0;
        fext[0] = lamb * 0;
        return fext;
      };

      /// neumann boundary load in vertical direction
      auto neumannBoundaryLoad = [&](auto &globalCoord, auto &lamb) {
        Eigen::Vector2d F = Eigen::Vector2d::Zero();
        F[1]              = lamb / 16.0;
        return F;
      };

      /// Python function which could be used to obtain the vertices at the right edge
      std::string lambdaNeumannVertices = std::string("lambda x: ( x[0]>47.9999 )");
      Python::start();
      Python::Reference main = Python::import("__main__");
      Python::run("import math");

      Python::runStream() << std::endl << "import sys" << std::endl << "import os" << std::endl;

      const auto &indexSet = gridView.indexSet();

      /// Flagging the vertices on which neumann load is applied as true
      Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
      auto pythonNeumannVertices = Python::make_function<bool>(Python::evaluate(lambdaNeumannVertices));

      for (auto &&vertex : vertices(gridView)) {
        bool isNeumann                          = pythonNeumannVertices(vertex.geometry().corner(0));
        neumannVertices[indexSet.index(vertex)] = isNeumann;
      }

      BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);

      for (auto &element : elements(gridView)) {
        auto localView = basis.localView();
        fes.emplace_back(basis, element, E, nu, &volumeLoad, &neumannBoundary, &neumannBoundaryLoad);
        fes.back().setEASType(numberOfEASParameters);
      }

      auto sparseAssembler = SparseFlatAssembler(basis, fes, dirichletFlags);

      auto KFunction = [&](auto &&disp, auto &&lambdaLocal) -> auto & {
        Ikarus::FErequirements req = FErequirementsBuilder()
                                         .insertGlobalSolution(Ikarus::FESolutions::displacement, disp)
                                         .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                         .addAffordance(Ikarus::MatrixAffordances::stiffness)
                                         .build();
        return sparseAssembler.getMatrix(req);
      };

      auto residualFunction = [&](auto &&disp, auto &&lambdaLocal) -> auto & {
        Ikarus::FErequirements req = FErequirementsBuilder()
                                         .insertGlobalSolution(Ikarus::FESolutions::displacement, disp)
                                         .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal)
                                         .addAffordance(Ikarus::VectorAffordances::forces)
                                         .build();
        return sparseAssembler.getVector(req);
      };

      Eigen::VectorXd D_Glob = Eigen::VectorXd::Zero(basis.size());

      auto startAssembly    = std::chrono::high_resolution_clock::now();
      auto nonLinOp         = Ikarus::NonLinearOperator(linearAlgebraFunctions(residualFunction, KFunction),
                                                        parameter(D_Glob, lambdaLoad));
      auto stopAssembly     = std::chrono::high_resolution_clock::now();
      auto durationAssembly = duration_cast<std::chrono::milliseconds>(stopAssembly - startAssembly);
      spdlog::info("The assembly took {:>6d} milliseconds with {} EAS parameters and {:>7d} dofs",
                   durationAssembly.count(), numberOfEASParameters, basis.size());
      ;
      timeVec.push_back(durationAssembly.count());
      const auto &K    = nonLinOp.derivative();
      const auto &Fext = nonLinOp.value();

      /// solve the linear system
      auto linSolver   = Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::sd_CholmodSupernodalLLT);
      auto startSolver = std::chrono::high_resolution_clock::now();

      linSolver.compute(K);
      linSolver.solve(D_Glob, -Fext);
      auto stopSolver     = std::chrono::high_resolution_clock::now();
      auto durationSolver = duration_cast<std::chrono::milliseconds>(stopSolver - startSolver);
      //          spdlog::info("The solver took {} milliseconds with {} EAS parameters and {} refinement level",
      //          durationSolver.count(),numberOfEASParameters,ref);

      /// Postprocess
      auto dispGlobalFunc
          = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis, D_Glob);
      Dune::VTKWriter vtkWriter(gridView, Dune::VTK::conforming);
      vtkWriter.addVertexData(dispGlobalFunc,
                              Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
      vtkWriter.write("Cook_MembraneConvergence" + std::to_string(ref));
      auto localView = basis.localView();
      auto localw    = localFunction(dispGlobalFunc);
      double uy_fe   = 0.0;
      Eigen::Vector2d req_pos;
      req_pos[0] = 48.0;
      req_pos[1] = 60.0;
      for (auto &ele : elements(gridView)) {
        localView.bind(ele);
        localw.bind(ele);
        const auto geo = localView.element().geometry();
        for (size_t i = 0; i < 4; ++i) {
          if ((geo.corner(i)[0] == req_pos[0]) and (geo.corner(i)[1] == req_pos[1])) {
            const auto local_pos = geo.local(toFieldVector(req_pos));
            uy_fe                = toEigenVector(localw(local_pos)).eval()[1];
          }
        }
      }
      dofsVec.push_back(basis.size());
      dispVec.push_back(uy_fe);

      auto stop     = std::chrono::high_resolution_clock::now();
      auto duration = duration_cast<std::chrono::milliseconds>(stop - start);
      spdlog::info("The total execution took {:>6d} milliseconds with {} EAS parameters and {:>7d} dofs",
                   duration.count(), numberOfEASParameters, basis.size());
      grid->globalRefine(1);
    }

    legends.push_back("Q1E" + std::to_string(easSet[nep]));
    auto p = ax->semilogx(dofsVec, dispVec);

    p->line_width(2);
    switch (easSet(nep)) {
      case 0:
        p->marker(line_spec::marker_style::asterisk);
        break;
      case 4:
        p->marker(line_spec::marker_style::circle);
        break;
      case 5:
        p->marker(line_spec::marker_style::cross);
        break;
      case 7:
        p->marker(line_spec::marker_style::diamond);
        break;
    }

    ax->hold(true);

    auto p2 = axesSecondPlot->semilogx(timeVec, dispVec);
    p2->line_width(2);
    p2->marker(line_spec::marker_style::asterisk);
    axesSecondPlot->hold(true);
  }
  ax->legend(legends);
  axesSecondPlot->legend(legends);
  auto legend  = ax->legend();
  auto legend2 = axesSecondPlot->legend();
  legend->location(legend::general_alignment::bottomright);
  legend2->location(legend::general_alignment::bottomright);
  f->show();
  f2->show();
}
