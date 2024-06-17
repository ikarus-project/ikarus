// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcommon.hh"

#include <memory>

#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/vtk/datacollectors/discontinuousdatacollector.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/vtkreader.hh>

#include <Eigen/Core>

#include <ikarus/assembler/simpleassemblers.hh>
#include "ikarus/finiteelements/feresulttypes.hh"
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mechanics/loads/volume.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/io/resultfunction.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/init.hh>
using Dune::TestSuite;

auto expectedStressValues(Dune::Vtk::DataTypes precision) -> std::vector<double> {
  if (precision == Dune::Vtk::DataTypes::FLOAT32)
    return {-0.612563, -3.06282,   -0.178688,  -0.605305,  -3.02653,   1.40567e-16, -0.165843,   -2.97347,
            -0.164172, -0.158585,  -2.93718,   0.0145161,  -0.605305,  -3.02653,    1.40567e-16, -0.612563,
            -3.06282,  0.178688,   -0.158585,  -2.93718,   -0.0145161, -0.165843,   -2.97347,    0.164172,
            0.250011,  -0.894201,  0.122203,   0.218461,   -1.05195,   0.0145161,   -0.0192071,  -0.948045,
            0.0591021, -0.0507578, -1.1058,    -0.0485852, 0.218461,   -1.05195,    -0.0145161,  0.250011,
            -0.894201, -0.122203,  -0.0507578, -1.1058,    0.0485852,  -0.0192071,  -0.948045,   -0.0591021};
  // Float64
  return {
      -0.612563441928246,    -3.06281720964123,    -0.1786881422641234,  -0.6053053722981661,  -3.02652686149083,
      1.405665430359793e-16, 0.2500112482565204,   -0.8942014658868807,  0.1222034745688248,   0.2184605681421391,
      -1.051954866458787,    -0.01451613926015979, -0.6053053722981661,  -3.02652686149083,    1.405665430359793e-16,
      -0.6125634419282462,   -3.06281720964123,    0.1786881422641237,   0.2184605681421391,   -1.051954866458787,
      -0.01451613926015979,  0.2500112482565203,   -0.894201465886881,   -0.1222034745688246,  0.2500112482565204,
      -0.8942014658868807,   0.1222034745688248,   0.2184605681421391,   -1.051954866458787,   -0.01451613926015979,
      -0.01920709001514156,  -0.948045133541213,   0.05910211434006233,  -0.05075777012952284, -1.105798534113119,
      0.04858522096860256,   0.2184605681421391,   -1.051954866458787,   -0.01451613926015979, 0.2500112482565203,
      -0.894201465886881,    -0.1222034745688246,  -0.05075777012952284, -1.105798534113119,   0.04858522096860256,
      -0.01920709001514165,  -0.9480451335412133,  -0.05910211434006224};
}

auto testFile(const std::string& fileName, Dune::Vtk::DataTypes precision) {
  TestSuite t("Test read in grid");
  auto epsilon = precision == Dune::Vtk::DataTypes::FLOAT32 ? 1e-6 : 1e-8;

  using Grid = Dune::UGGrid<2>;
  Dune::GridFactory<Grid> factory;

  Dune::Vtk::VtkReader reader{factory};
  reader.read(fileName);
  auto grid     = reader.createGrid();
  auto gridView = grid->leafGridView();

  t.check(gridView.size(0) == 4) << "GridView should have 4 elements, but has " << gridView.size(0) << " elements";

  // get point Data
  auto stressData = reader.getPointData("linearStress");

  t.check(stressData.numComponents() == 3) << "Num components should be 3, but is " << stressData.numComponents();
  t.check(stressData.dataType() == precision) << "Precision is wrong";
  std::cout << "New\n";
  auto localStressFunction = localFunction(stressData);
  auto results             = expectedStressValues(precision);
  for (int i = 0; const auto& ele : Dune::elements(gridView)) {
    localStressFunction.bind(ele);
    auto refEle = Dune::referenceElement(ele);
    for (const auto c : Dune::range(ele.geometry().corners())) {
      auto corner       = refEle.geometry<2>(c).center();
      auto stressInFile = localStressFunction(corner);
      for (const auto j : Dune::range(3)) {
        t.check(Dune::FloatCmp::eq(stressInFile[j], results[i], epsilon))
            << "Result read in is " << stressInFile[j] << "but should be " << results[i];
        ++i;
      }
    }
  }

  return t;
}

auto runTestCase() {
  TestSuite t("Test ResultFunction");
  std::string fileName = "ResultFunctionTest";

  using Grid = Dune::UGGrid<2>;
  using GridView = Grid::LeafGridView;

  constexpr double Lx                     = 4.0;
  constexpr double Ly                     = 4.0;
  const Dune::FieldVector<double, 2> bbox = {Lx, Ly};
  const std::array<unsigned int, 2> elementsPerDirection   = {2, 2};

  auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid({0, 0}, bbox, elementsPerDirection);
  auto gridView                           = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>()));

  Ikarus::DirichletValues dirichletValues(basis.flat());
  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[1]) < 1e-8)
      dirichletFlags[localView.index(localIndex)] = true;
  });

  auto vL      = []([[maybe_unused]] auto& globalCoord, auto& lamb) { return Eigen::Vector2d{0, -1}; };
  auto skills_ = Ikarus::skills(Ikarus::linearElastic({.emodul = 100, .nu = 0.2}), Ikarus::volumeLoad<2>(vL));

  using LinearElastic = decltype(Ikarus::makeFE(basis, skills_));
  std::vector<LinearElastic> fes;

  for (auto&& element : elements(gridView)) {
    fes.emplace_back(Ikarus::makeFE(basis, skills_));
    fes.back().bind(element);
  }

  /// Create a sparse assembler
  auto sparseAssembler = Ikarus::SparseFlatAssembler(fes, dirichletValues);

  Eigen::VectorXd D_Glob = Eigen::VectorXd::Zero(basis.flat().size());

  auto req = Ikarus::FERequirements().addAffordance(Ikarus::AffordanceCollections::elastoStatics);

  auto residualFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return sparseAssembler.getVector(req);
  };

  auto KFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return sparseAssembler.getMatrix(req);
  };

  double lambdaLoad = 1.0;
  auto nonLinOp =
      Ikarus::NonLinearOperator(Ikarus::functions(residualFunction, KFunction), Ikarus::parameter(D_Glob, lambdaLoad));
  const auto& K    = nonLinOp.derivative();
  const auto& Fext = nonLinOp.value();

  auto linSolver = Ikarus::LinearSolver(Ikarus::SolverTypeTag::sd_CholmodSupernodalLLT);
  linSolver.compute(K);
  linSolver.solve(D_Glob, -Fext);
  req.insertGlobalSolution(Ikarus::FESolutions::displacement, D_Glob);

  auto stressFunction    = Ikarus::makeResultFunction<Ikarus::ResultTypes::linearStress>(&fes, req);
  auto stressVtkFunction = Ikarus::makeResultVtkFunction<Ikarus::ResultTypes::linearStress>(&fes, req);

  // Test
  // Dune::VTKWriter will always default to FLOAT32
  Dune::VTKWriter vtkWriter(gridView, Dune::VTK::nonconforming);
  vtkWriter.addVertexData(stressFunction);
  vtkWriter.write(fileName + "_Common_VTK_32");
  t.subTest(testFile(fileName + "_Common_VTK_32" + ".vtu", Dune::Vtk::DataTypes::FLOAT32));

  // Passing float64 doens't change it
  Dune::VTKWriter vtkWriter4(gridView, Dune::VTK::nonconforming, Dune::VTK::Precision::float64);
  vtkWriter4.addVertexData(stressFunction);
  vtkWriter4.write(fileName + "_Common_VTK_64");
  t.subTest(testFile(fileName + "_Common_VTK_64" + ".vtu", Dune::Vtk::DataTypes::FLOAT32));

  // Vtk::VtkWriter can be told to use a special Datatype
  Dune::Vtk::UnstructuredGridWriter<GridView, Dune::Vtk::DiscontinuousDataCollector<GridView>>
      vtkWriter2(gridView, Dune::Vtk::FormatTypes::ASCII, Dune::Vtk::DataTypes::FLOAT32);
  vtkWriter2.addPointData(stressVtkFunction);
  auto vtkfileName = vtkWriter2.write(fileName + "_Vtk_32");
  t.subTest(testFile(vtkfileName, Dune::Vtk::DataTypes::FLOAT32));

  Dune::Vtk::UnstructuredGridWriter vtkwriter5(gridView);
  vtkwriter5.write(fileName + "_5");
  //
  // Dune::Vtk::UnstructuredGridWriter vtkWriter3(gridView, Dune::Vtk::FormatTypes::ASCII, Dune::Vtk::DataTypes::FLOAT64);
  // vtkWriter3.addPointData(stressVtkFunction);
  // auto vtkfileName3 = vtkWriter3.write(fileName + "_Vtk_64");
  // t.subTest(testFile(vtkfileName3, Dune::Vtk::DataTypes::FLOAT64));

  return t;
}

int main(const int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(runTestCase());

  return t.exit();
}