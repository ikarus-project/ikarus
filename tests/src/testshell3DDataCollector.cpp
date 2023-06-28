// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testCommon.hh"
#include "testHelpers.hh"

#include <dune/common/parametertreeparser.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/iga/nurbsbasis.hh>

#include "spdlog/spdlog.h"

#include <Eigen/Core>

#include <ikarus/assembler/simpleAssemblers.hh>
#include <ikarus/controlRoutines/loadControl.hh>
#include <ikarus/controlRoutines/pathFollowingTechnique.hh>
#include <ikarus/finiteElements/mechanics/fesettings.hh>
#include <ikarus/finiteElements/mechanics/kirchhoffloveshell.hh>
#include <ikarus/io/resultFunction.hh>
#include <ikarus/linearAlgebra/dirichletValues.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/solver/nonLinearSolver/trustRegion.hh>
#include <ikarus/utils/algorithms.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controlVTKWriter.hh>
#include <ikarus/utils/observer/nonLinearSolverLogger.hh>
#include <ikarus/io/shell3DDataCollector.hh>

using Dune::TestSuite;
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <dune/geometry/virtualrefinement.hh>
auto test3DDataCollector() {
  TestSuite t("test3DDataCollector ");
  constexpr auto dimworld        = 3;
  const std::array<int, 2> order = {1, 1};

  const std::array<std::vector<double>, 2> knotSpans = {{{0, 0, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0, 0}, .w = 1}, {.p = {0, 2, 0}, .w = 1}}, {{.p = {10, 0, 0}, .w = 1}, {.p = {10, 2, 0}, .w = 1}}};

  std::array<int, 2> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<2, dimworld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<2, dimworld>;

  Dune::IGA::NURBSPatchData<2, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  for (int i = 0; i < 2; ++i)
    patchData = degreeElevate(patchData, i, 1);

  Dune::ParameterTree parameterSet;
  Dune::ParameterTreeParser::readINITree("/tmp/Ikarus/tests/src/shell.parset",
                                         parameterSet);

  const auto E              = parameterSet.get<double>("E");
  const auto nu             = parameterSet.get<double>("nu");
  const auto thickness      = parameterSet.get<double>("thickness");
  const auto loadFactor     = parameterSet.get<double>("loadFactor");
  const auto simulationFlag = parameterSet.get<int>("simulationFlag");
  const auto refine         = parameterSet.get<int>("refine");
  const auto plotInPlaneRefine         = parameterSet.get<int>("plotInPlaneRefine");
  auto grid = std::make_shared<Grid>(patchData);

  grid->globalRefineInDirection(0,refine);
  auto gridView = grid->leafGridView();



  using GridView = decltype(gridView);

  Dune::Vtk::Shell3DDataCollector dataCollector1(gridView,1,Dune::RefinementIntervals(plotInPlaneRefine));

  Dune::VtkUnstructuredGridWriter writer2(dataCollector1, Dune::Vtk::FormatTypes::ASCII);
  writer2.write("KLSHELL3DTEST");

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);


  test3DDataCollector();

  typedef Dune::VirtualRefinement<2, double> Refinement;
  Refinement& ref = Dune::buildRefinement<2, double>(Dune::GeometryTypes::quadrilateral, Dune::GeometryTypes::quadrilateral);
  std::cout << "Index\tEcke0\tEcke1\tEcke2\tEcke3" << std::endl;
  Refinement::ElementIterator eEnd = ref.eEnd(Dune::RefinementIntervals(2));
  {
    for (Refinement::ElementIterator i = ref.eBegin(Dune::RefinementIntervals(2)); i != eEnd; ++i) {
      std::cout << i.index() << "\t" << i.vertexIndices()[0] << "\t" << i.vertexIndices()[1] << "\t" << i.vertexIndices()[2]
           << "\t" << i.vertexIndices()[3] <<  std::endl;
      std::cout << std::endl;
    }
  }

  std::cout << "Index\tx\ty" << std::endl;
     Refinement::VertexIterator vEnd = ref.vEnd(Dune::RefinementIntervals(2));
    for(Refinement::VertexIterator i = ref.vBegin(Dune::RefinementIntervals(2)); i != vEnd; ++i) {

    std::cout << i.index() << "\t" << i.coords()[0] << "\t" << i.coords()[1] << std::endl;
     }
     std::cout << std::endl;

}
