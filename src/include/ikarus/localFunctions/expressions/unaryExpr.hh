// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

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
