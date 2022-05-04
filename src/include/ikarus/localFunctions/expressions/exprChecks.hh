//
// Created by Alex on 04.05.2022.
//

#pragma once
#include <ikarus/localFunctions/leafNodeCollection.hh>
namespace Ikarus {

template <typename LocalFunctionImpl>
bool checkIfAllLeafNodeHaveTheSameBasisState(const LocalFunctionImpl& lf)  {
  using namespace Dune::Indices;
  auto leafNodeCollection = collectLeafNodeLocalFunctions(lf);
  bool isValid            = true;
  if constexpr (leafNodeCollection.size() > 0) {
//    std::cout<<"Test"<<std::endl;
    const bool isBound = leafNodeCollection.node(_0).basis().isBound();
//    std::cout<<isBound<<std::endl;
//    std::cout<<"Size: "<<leafNodeCollection.size()<<std::endl;
    unsigned int integrationRuleSize = isBound ? leafNodeCollection.node(_0).basis().integrationPointSize() : 0;
//    std::cout<<"integrationRuleSize: "<<integrationRuleSize<<std::endl;
    Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<leafNodeCollection.size()>{}),
                          [&]<typename I>(I &&i) {
                            if constexpr (I::value==0) { //Skip first value
                            } else {
                              auto nodeBasis = leafNodeCollection.node(i).basis();
//                              std::cout<<"I: "<<I::value<<std::endl;
                              if (nodeBasis.isBound()!=isBound)
                                isValid = false;
                              else {
                                if (nodeBasis.integrationPointSize()!=integrationRuleSize)
                                  isValid = false;
                              }
                            }

                          });
  }
  return isValid;
}

} // Ikarus
