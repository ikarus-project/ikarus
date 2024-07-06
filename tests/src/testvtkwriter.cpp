// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "dummyproblem.hh"
#include "testhelpers.hh"

#include <memory>

#include "dune/common/indices.hh"
#include "dune/functions/functionspacebases/subspacebasis.hh"
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/common/entity.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/vtk/datacollectors/discontinuousdatacollector.hh>
#include <dune/vtk/types.hh>
#include <dune/vtk/vtkreader.hh>

#include <Eigen/Core>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/feresulttypes.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mechanics/loads/volume.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/io/resultfunction.hh>
#include <ikarus/io/vtkwriter.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

template <typename GridView, typename Assembler>
auto testInstantiationAndTemplateArgumentDeduction(const GridView& gridView, std::shared_ptr<Assembler> assembler) {
  static_assert(Ikarus::Concepts::FlatAssembler<typename decltype(assembler)::element_type>);
  // Create Vtk::Writer
  static_assert(std::is_class_v<Ikarus::Vtk::Writer<Assembler, false>>);
  static_assert(std::is_same_v<typename Ikarus::Vtk::Writer<Assembler, false>::DataCollector,
                               Dune::Vtk::ContinuousDataCollector<GridView>>);

  // Unfortionatly we can instantiate the structured Writer and DataCollector even for a unstructured Grid (UG)
  static_assert(std::is_class_v<Ikarus::Vtk::Writer<Assembler, true>>);
  static_assert(std::is_same_v<typename Ikarus::Vtk::Writer<Assembler, true>::DataCollector,
                               Dune::Vtk::YaspDataCollector<GridView>>);

  auto writer = Ikarus::Vtk::Writer<Assembler>(assembler);
  auto writerArgs =
      Ikarus::Vtk::Writer<Assembler>(assembler, Dune::Vtk::FormatTypes::BINARY, Dune::Vtk::DataTypes::FLOAT32);

  auto writerDC     = Ikarus::Vtk::Writer<Assembler, false, Dune::Vtk::DiscontinuousDataCollector<GridView>>(assembler);
  auto writerDCArgs = Ikarus::Vtk::Writer<Assembler, false, Dune::Vtk::DiscontinuousDataCollector<GridView>>(
      assembler, Dune::Vtk::FormatTypes::BINARY, Dune::Vtk::DataTypes::FLOAT32);

  Dune::Vtk::DiscontinuousDataCollector<GridView> dc{gridView};
  static_assert(Ikarus::Concepts::DataCollector<decltype(dc)>);

  auto writerDCAsArg =
      Ikarus::Vtk::Writer<Assembler, false, Dune::Vtk::DiscontinuousDataCollector<GridView>>(assembler, dc);
  auto writerDCAsArgM =
      Ikarus::Vtk::Writer<Assembler, false, Dune::Vtk::DiscontinuousDataCollector<GridView>>(assembler, std::move(dc));
  auto writerDCAsArgArgs = Ikarus::Vtk::Writer<Assembler, false, Dune::Vtk::DiscontinuousDataCollector<GridView>>(
      assembler, dc, Dune::Vtk::FormatTypes::BINARY, Dune::Vtk::DataTypes::FLOAT32);

  // using CTAD
  auto writerCTAD     = Ikarus::Vtk::Writer(assembler);
  auto writerArgsCTAD = Ikarus::Vtk::Writer(assembler, Dune::Vtk::FormatTypes::BINARY, Dune::Vtk::DataTypes::FLOAT32);
  auto writerDCAsArgCTAD  = Ikarus::Vtk::Writer(assembler, dc);
  auto writerDCAsArgMCTAD = Ikarus::Vtk::Writer(assembler, std::move(dc));
  auto writerDCAsArgArgsCTAD =
      Ikarus::Vtk::Writer(assembler, dc, Dune::Vtk::FormatTypes::BINARY, Dune::Vtk::DataTypes::FLOAT32);
}

auto runTest() {
  TestSuite t("Test ResultFunction");
  std::string fileName = "ResultFunctionTest";

  using Grid     = Dune::UGGrid<2>;
  using GridView = Grid::LeafGridView;

  DummyProblem<Grid> testCase{};

  auto& gridView       = testCase.gridView();
  auto sparseAssembler = testCase.sparseAssembler();
  auto& req            = testCase.requirement();
  auto& basis          = testCase.basis();
  auto& D_Glob         = req.globalSolution();

  // Tests
  testInstantiationAndTemplateArgumentDeduction(gridView, sparseAssembler);

  Dune::Vtk::DiscontinuousDataCollector dc{gridView};

  auto writer = Ikarus::Vtk::Writer(sparseAssembler);

  using ::asCellData;
  using ::asPointData;

  writer.setDatatype(Dune::Vtk::DataTypes::FLOAT64);
  writer.setFormat(Dune::Vtk::FormatTypes::ASCII);

  writer.addResult<Ikarus::ResultTypes::linearStress>(asPointData);
  writer.addResultFunction(Ikarus::makeResultFunction<Ikarus::ResultTypes::linearStress>(sparseAssembler), asCellData);

  writer.addInterpolation<2>(D_Glob, basis.flat(), "displacement", asPointData);
  writer.addInterpolation<1>(D_Glob, Dune::Functions::subspaceBasis(basis.flat(), Dune::index_constant<0>()),
                             "displacement_u", asPointData);

  writer.addPointData(
      Dune::Functions::makeDiscreteGlobalBasisFunction<Dune::FieldVector<double, 2>>(basis.flat(), D_Glob),
      Dune::Vtk::FieldInfo("displacements_gf", 2, Dune::Vtk::RangeTypes::VECTOR));

  auto vtkFileName = writer.write(fileName);

  Dune::GridFactory<Grid> factory;
  Dune::Vtk::VtkReader reader{factory};

  reader.read(vtkFileName);
  auto gridRead     = reader.createGrid();
  auto gridViewRead = gridRead->leafGridView();

  t.check(gridView.size(0) == gridViewRead.size(0));

  auto stressData = reader.getPointData("linearStress");
  t.check(stressData.numComponents() == 3)
      << testLocation() << "Num components should be 3, but is " << stressData.numComponents();
  t.check(stressData.dataType() == Dune::Vtk::DataTypes::FLOAT64)
      << testLocation() << std::source_location::current().line() << "Precision should be float64";

  auto stressDataCell = reader.getCellData("linearStress");
  t.check(stressDataCell.numComponents() == 3)
      << testLocation() << "Num components should be 3, but is " << stressDataCell.numComponents();
  t.check(stressDataCell.dataType() == Dune::Vtk::DataTypes::FLOAT64)
      << testLocation() << "Precision is should be float64";

  // Writing out as a vector does always seem to yield 3 (or potentially more) entries
  auto displacementData = reader.getPointData("displacement");
  t.check(displacementData.numComponents() == 3)
      << testLocation() << "Num components should be 3, but is " << displacementData.numComponents();

  auto displacementDataGF = reader.getPointData("displacements_gf");
  t.check(displacementDataGF.numComponents() == 3)
      << testLocation() << "Num components should be 3, but is " << displacementDataGF.numComponents();

  auto writer2 = Ikarus::Vtk::Writer(sparseAssembler, Dune::Vtk::DiscontinuousDataCollector<GridView>{gridView});
  writer2.addAllResults(asPointData);
  writer2.addAllResults(asCellData);
  auto vtkFileName2 = writer2.write(fileName + "_2");

  reader.read(vtkFileName2);
  t.checkNoThrow([&]() { reader.getPointData("linearStress"); });
  t.checkNoThrow([&]() { reader.getCellData("linearStress"); });

  return t;
}

int main(const int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(runTest());

  return t.exit();
}