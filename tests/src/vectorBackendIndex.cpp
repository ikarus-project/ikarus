// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later


#include "config.h"

#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/grid/yaspgrid.hh>

using Dune::TestSuite;

auto multiIndexTest() {
  TestSuite t("SimpleAssemblersTest");
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox       = {4, 2};
  std::array<int, 2> elementsPerDirection = {2, 1};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);

    auto gridView = grid->leafGridView();

    using namespace Dune::Functions::BasisFactory;
    auto basis = makeBasis(gridView, power<2>(lagrange<1>(), BlockedInterleaved()));


  std::vector<double> values;
  values.resize(basis.dimension());
  auto valBackend =Dune::Functions::istlVectorBackend(values);
  auto localView = basis.localView();
  for (int eleIndex=0; auto& ele :elements(gridView)) {
    localView.bind(ele);
    std::cout<<"Ele: "<<eleIndex<<std::endl;
    const auto& fe = localView.tree().child(0).finiteElement();
    for (size_t i = 0; i < fe.size(); ++i) {
      for (int j = 0; j < 2; ++j) {
        auto globalIndex=localView.index(localView.tree().child(j).localIndex(i));
        std::cout<<"i: "<<i<<" j: "<<j<< " gI: "<<globalIndex<<std::endl;
        auto valIterator = &valBackend[globalIndex];
        std::cout<<valIterator-values.data()<<std::endl;
      }
    }
++eleIndex;
  }
  return t;
}

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t;

  t.subTest(multiIndexTest());
  return t.exit();
}
