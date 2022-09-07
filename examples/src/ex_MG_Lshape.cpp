//
// Created by Alex on 21.07.2021.
//

#include <config.h>

#include <dune/alugrid/grid.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
//#include <dune/grid/yaspgrid.hh>
//#include <dune/iga/nurbsgrid.hh>
#include <functional>
#include <algorithm>
#include <chrono>

#include <dune/grid/uggrid.hh>

#include "spdlog/spdlog.h"

#include <Eigen/Core>

#include <ikarus/solver/linearSolver/geometricMultigrid/geometricMultiGridSolver.hh>
#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlRoutines/loadControl.hh>
#include <ikarus/finiteElements/feRequirements.hh>
#include <ikarus/finiteElements/mechanics/nonLinearElasticityFE.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/manifolds/realTuple.hh>
#include <ikarus/solver/linearSolver/geometricMultigrid/gridTransfer.h>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/solver/nonLinearSolver/trustRegion.hh>
#include <ikarus/utils/algorithms.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <ikarus/utils/observer/nonLinearSolverLogger.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/common/parametertreeparser.hh>

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);

  Dune::ParameterTree parameterSet;
  Dune::ParameterTreeParser::readINITree(argv[1], parameterSet);

  const Dune::ParameterTree &gridParameters     = parameterSet.sub("GridParameters");
  const Dune::ParameterTree &basisParameters     = parameterSet.sub("BasisParameters");
  const auto refinement      = gridParameters.get<int>("refinement");
  const auto meshType      = gridParameters.get<int>("meshType");
  const auto basisOrder      = basisParameters.get<int>("basisOrder");


  using namespace Ikarus;
  constexpr int gridDim = 2;
  //  //  /// ALUGrid Example
    using Grid = Dune::UGGrid<gridDim>;
    Dune::GridFactory<Grid> gridFactory;

    const double L = 1.0;

    gridFactory.insertVertex({0, 0});
    gridFactory.insertVertex({L/2, 0});
    gridFactory.insertVertex({0, L/2});
    gridFactory.insertVertex({L/2, L/2});
    gridFactory.insertVertex({0, L});
    gridFactory.insertVertex({L/2, L});
    gridFactory.insertVertex({L, L});
    gridFactory.insertVertex({L, L/2});


    gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {0, 1, 2, 3});
    gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {2, 3,4,5});
    if(meshType == 0)
    gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {3,7,5,6});
    else if(meshType == 1)
    {
      gridFactory.insertElement(Dune::GeometryTypes::triangle, {3,7,6});
      gridFactory.insertElement(Dune::GeometryTypes::triangle, {3,6,5});
    }
    auto grid     = gridFactory.createGrid();
    grid->globalRefine(refinement);
    auto leafGridView = grid->leafGridView();

//    draw(leafGridView);
//    draw(grid->levelGridView(1));
//    draw(grid->levelGridView(0));

//
  using namespace Dune::Functions::BasisFactory;

    auto preBasisFactory = power<gridDim>(lagrange(basisOrder), FlatInterleaved());



    auto coarseBasis = makeBasis(grid->levelGridView(0),preBasisFactory);
    auto fineBasis = makeBasis(grid->leafGridView(),preBasisFactory);
    spdlog::info("Dofs on finest level: {}",fineBasis.size());
      double lambdaLoad = 1;

        std::vector<Ikarus::NonLinearElasticityFE<decltype(coarseBasis)>> feVectorCoarse;
        std::function<Eigen::Vector<double, 2>(const Eigen::Vector<double, 2>&,
                                                              const double&)> volumeLoad_ =[](auto& globalCoord, auto& lamb) {
          Eigen::Vector2d fext;
          fext.setZero();
          fext[1] = 2 * lamb;
          fext[0] = 0;
          return fext;
        };
        typename Ikarus::NonLinearElasticityFE<decltype(coarseBasis)>::Settings settings({.emod_=1000, .nu_=0.3, .volumeLoad=volumeLoad_});
        for (auto& element : elements(coarseBasis.gridView()))
          feVectorCoarse.emplace_back(coarseBasis, element, settings);
        auto startSolver = std::chrono::high_resolution_clock::now();
      Ikarus::GeometricMultiGridSolver solver(grid.get(),preBasisFactory,feVectorCoarse);



      Eigen::VectorXd d,FextFine;
      solver.solve(d,(-FextFine).eval());
      auto stopSolver = std::chrono::high_resolution_clock::now();
      auto durationSolver = duration_cast<std::chrono::milliseconds>(stopSolver - startSolver);
      spdlog::info("The solver took {} milliseconds",durationSolver.count());
      Eigen::VectorXd dFull;
      solver.transformToFineFull(d,dFull);



      /// Postprocess
      auto disp
          = Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(fineBasis, dFull);
      Dune::VTKWriter vtkWriter(fineBasis.gridView(), Dune::VTK::conforming);
      vtkWriter.addVertexData(disp, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));

      vtkWriter.write("LShapeMultigrid");

      spdlog::info( "dofs: {}, maximum displacement: {}",fineBasis.size(),std::ranges::max(dFull));

}