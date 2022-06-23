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

  template <template <typename, typename> class Op, typename E1, typename E2>
  struct BinaryLocalFunctionExpression : public Ikarus::LocalFunctionInterface<Op<E1, E2>> {
    using E1Raw = std::remove_cvref_t<E1>;
    using E2Raw = std::remove_cvref_t<E2>;

    constexpr const E1& l() const { return std::get<0>(expr); }
    constexpr const E2& r() const { return std::get<1>(expr); }

    constexpr E1& l() { return std::get<0>(expr); }

    constexpr E2& r() { return std::get<1>(expr); }

    auto clone() const { return Op<decltype(l().clone()), decltype(r().clone())>(l().clone(), r().clone()); }

    /** Rebind the value type of the underlying local function with the id ID */
    template <typename OtherType, size_t ID = 0>
    auto rebindClone(OtherType&&, Dune::index_constant<ID>&& id = Dune::index_constant<0UL>()) const {
      return rebind<Op, E1, E2, OtherType>(l(), r(), Dune::index_constant<ID>());
    }

    /* Creates a tuple of all subtype ids, size l or r is not a tuple, tuple_cat may not work.
     * Thus we artifically wrap them inside a tuple  */
    using Ids = decltype(Std::makeNestedTupleFlat(
        std::make_tuple(std::declval<typename E1Raw::Ids>(), std::declval<typename E2Raw::Ids>())));

    /** The function order wrt. the coefficients */
    template <size_t ID_ = 0>
    static constexpr int orderID = Op<E1, E2>::template orderID<ID_>;

    constexpr BinaryLocalFunctionExpression(E1&& u, E2&& v) requires IsLocalFunction<E1, E2>
        : expr(std::forward<E1>(u), std::forward<E2>(v)) {}

    static constexpr bool isLeaf  = false;
    static constexpr int children = 2;

  private:
    std::tuple<E1, E2> expr;
  };

}  // namespace Ikarus