//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>

#include "testHelpers.h"

//#include <fstream>
#include <vector>

#include <dune/typetree/leafnode.hh>
//#include <Eigen/Core>
#include <dune/grid/yaspgrid.hh>
//#include <ikarus/FEManager/DefaultFEManager.h>
//#include <ikarus/FiniteElements/ElasticityFE.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
//#include <ikarus/Geometries/GeometryType.h>
//#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/typetree/powernode.hh>

GTEST_TEST(Basis, Basistest) {
  using Grid        = Dune::YaspGrid<2>;
  const double L    = 1;
  const double h    = 1;
  const size_t elex = 2;
  const size_t eley = 1;

  Dune::FieldVector<double, 2> bbox = {L, h};
  std::array<int, 2> eles           = {elex, eley};
  auto grid                         = std::make_shared<Grid>(bbox, eles);
  auto gridView                     = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  using namespace Dune::Indices;
  constexpr int p = 1;
  auto basis      = makeBasis(gridView,
                              composite(power<2>(lagrange<p>(), FlatInterleaved()), lagrange<p - 1>(), FlatLexicographic()));
  EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, decltype(basis)::PreBasis::SubPreBases>::IndexMergingStrategy,
                              Dune::Functions::BasisFactory::FlatInterleaved>));
  auto dispBasis               = subspaceBasis(basis, _0);
  auto pressureBasis           = subspaceBasis(basis, _1);
  auto localView               = basis.localView();
  auto localViewOfDisplacement = dispBasis.localView();
  auto localViewOfPressure     = pressureBasis.localView();

  for (auto& e : elements(gridView)) {
    localView.bind(e);
    localViewOfDisplacement.bind(e);
    localViewOfPressure.bind(e);

    EXPECT_EQ(localView.size(), 9);                         // Total Ansatzfunctions (Dofs)
    EXPECT_EQ(localViewOfDisplacement.size(), 9);           // Total Ansatzfunctions (Dofs)
    EXPECT_EQ(localViewOfPressure.size(), 9);               // Total Ansatzfunctions (Dofs)
    EXPECT_EQ(localView.tree().size(), 9);                  // Total Ansatzfunctions (Dofs)
    EXPECT_EQ(localViewOfDisplacement.tree().size(), 8);    // Displacement Ansatzfunctions (DispDofs)
    EXPECT_EQ(localViewOfPressure.tree().size(), 1);        // Pressure Ansatzfunctions (PressureDofs)
    EXPECT_EQ(localViewOfDisplacement.tree().degree(), 2);  // How many displacement childs are there?
    EXPECT_EQ(localViewOfPressure.tree().degree(), 0);      // How many pressure childs are there?
    EXPECT_EQ(localViewOfDisplacement.tree().child(0).degree(), 0);
    EXPECT_EQ(localViewOfDisplacement.tree().child(1).degree(), 0);
    EXPECT_EQ(localViewOfDisplacement.globalBasis().dimension(), 14);  // How many dofs do we have in total?
    EXPECT_EQ(basis.dimension(), basis.size());                        // size and dimension cooincide for Flat basis
  }
}
