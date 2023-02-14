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
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/typetree/treecontainer.hh>

using Dune::TestSuite;

template<typename Basis, int size>
size_t calcFlatindex(Basis& basis,Dune::ReservedVector<long unsigned int, size> globalIndex)
{
//  using Node = std::remove_cvref_t<decltype(basis.preBasis().subPreBasis())>::Node;
//  Dune::ReservedVector<size_t,globalIndex.size>
//  if constexpr  (not Node::isLeaf)
    auto offset = basis.size({globalIndex})*globalIndex.back();
    auto res = offset;
    for (int i = 0; i < globalIndex.back(); ++i) {
      res
    }
    globalIndex.pop_back();
    if (not globalIndex.empty())
    return offset+calcFlatindex(basis,globalIndex);
    else
     return offset;
//  auto result = globalIndex[0];
//  if constexpr (Node::isComposite) {
//    Dune::Hybrid::forEach(std::make_index_sequence<Node::children>(), [&](auto i) {
//      result *= calcFlatindex(basis.preBasis().subPreBasis(i)) * globalIndex[i];
//      });
//  } else if (Node::isPower) {
//    for (int i = 0; i < Node::children; ++i) {
//      result *=  basis.preBasis().subPreBasis(i).dimension * globalIndex[i];
//    }
//  }else if (Node::isLeaf)
//  {
//    result+=globalIndex[index];
//  }
//  return result;
}


auto multiIndexTest() {
  TestSuite t("SimpleAssemblersTest");
  using Grid = Dune::YaspGrid<2>;

  Dune::FieldVector<double, 2> bbox       = {4, 2};
  std::array<int, 2> elementsPerDirection = {2, 1};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);

    auto gridView = grid->leafGridView();

    using namespace Dune::Functions::BasisFactory;
    auto basis = makeBasis(gridView, composite(power<2>(lagrange<1>(), BlockedLexicographic()),lagrange<2>(), BlockedLexicographic()));


  std::vector<double> values;
  values.resize(basis.dimension());
  auto& indexSet = gridView.indexSet();

  for (auto& vert : vertices(gridView)) {
    values[indexSet.index(vert)]= indexSet.index(vert);
  }
  auto valBackend =Dune::Functions::istlVectorBackend(values);
//  valBackend.resize(basis);

  std::vector<int> vec;
  auto localView= basis.localView();
//  auto ele0=elements(gridView).begin();
//  localView.bind(*ele0);
  auto vTest =Dune::TypeTree::makeTreeContainer<int>(localView.tree());
//  Dune::TypeTree::forEachNode()

  for (int eleIndex=0; auto& ele :elements(gridView)) {
    localView.bind(ele);
    std::cout<<"Ele: "<<eleIndex<<std::endl;
    using namespace Dune::Indices;
    const auto& fe = localView.tree().child(_0,0).finiteElement();
    for (size_t i = 0; i < fe.size(); ++i) {
      for (int j = 0; j < 2; ++j) {
        auto globalIndex=localView.index(localView.tree().child(_0,j).localIndex(i));
        std::cout<<"i: "<<i<<" j: "<<j<< " gI: "<<globalIndex<<" "<<calcFlatindex(basis,globalIndex)<<std::endl;
//        int res=0;
//        static_assert(globalIndex.capacity()==3);
//        auto vTest =Dune::TypeTree::makeTreeContainer<int>(localView.tree());
//        auto globalFlatIndex= 0;
      }
    }
    const auto& fe2 = localView.tree().child(_1).finiteElement();
    for (size_t i = 0; i < fe2.size(); ++i) {
        auto globalIndex=localView.index(localView.tree().child(_1).localIndex(i));
        std::cout<<"i: "<<i<< " gI: "<<globalIndex<<" "<<calcFlatindex(basis,globalIndex)<<std::endl;
//        int res=0;
//        static_assert(globalIndex.capacity()==3);
//        auto vTest =Dune::TypeTree::makeTreeContainer<int>(localView.tree());
//        auto globalFlatIndex= 0;

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
