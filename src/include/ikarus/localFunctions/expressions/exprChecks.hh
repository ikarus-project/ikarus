// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once
#include <ikarus/localFunctions/leafNodeCollection.hh>
namespace Ikarus {

  template <typename LocalFunctionImpl>
  bool checkIfAllLeafNodeHaveTheSameBasisState(const LocalFunctionImpl& lf) {
    using namespace Dune::Indices;
    auto leafNodeCollection = collectLeafNodeLocalFunctions(lf);
    bool isValid            = true;
    if constexpr (leafNodeCollection.size() > 0) {
      const bool isBound               = leafNodeCollection.node(_0).basis().isBound();
      unsigned int integrationRuleSize = isBound ? leafNodeCollection.node(_0).basis().integrationPointSize() : 0;
      Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<leafNodeCollection.size()>{}),
                            [&]<typename I>(I&& i) {
                              if constexpr (I::value == 0) {  // Skip first value
                              } else {
                                auto nodeBasis = leafNodeCollection.node(i).basis();
                                if (nodeBasis.isBound() != isBound)
                                  isValid = false;
                                else {
                                  if (nodeBasis.integrationPointSize() != integrationRuleSize) isValid = false;
                                }
                              }
                            });
    }
    return isValid;
  }

}  // namespace Ikarus
