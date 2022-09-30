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

#include "rebind.hh"

#include <ikarus/localFunctions/localFunctionInterface.hh>

namespace Ikarus {

  template <template <typename, typename...> class Op, typename E1, typename... Args>
  struct UnaryLocalFunctionExpression : public Ikarus::LocalFunctionInterface<Op<E1, Args...>> {
    std::tuple<E1> expr;

    using E1Raw = std::remove_cvref_t<E1>;

    template <size_t ID_ = 0>
    static constexpr int orderID = Op<E1, Args...>::template orderID<ID_>;

    const E1& m() const { return std::get<0>(expr); }

    E1& m() { return std::get<0>(expr); }

    using Ids = decltype(Std::makeNestedTupleFlat(std::make_tuple(std::declval<typename E1Raw::Ids>())));

    auto clone() const { return Op<decltype(m().clone()), Args...>(m().clone()); }

    /** Rebind the value type of the underlying local function with the id ID */
    template <typename OtherType, size_t ID = 0>
    auto rebindClone(OtherType&&, [[maybe_unused]] Dune::index_constant<ID>&& id = Dune::index_constant<0>()) const {
      return rebind<Op, E1, OtherType, ID, Args...>(m(), Dune::index_constant<ID>());
    }

    constexpr explicit UnaryLocalFunctionExpression(E1&& u) requires IsLocalFunction<E1> : expr(std::forward<E1>(u)) {}

    static constexpr bool isLeaf  = false;
    static constexpr int children = 1;
  };

}  // namespace Ikarus