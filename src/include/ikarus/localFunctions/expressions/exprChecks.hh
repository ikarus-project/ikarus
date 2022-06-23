/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */



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
