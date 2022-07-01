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
#include <cstddef>
#include <tuple>

#include <ikarus/utils/traits.hh>
namespace Ikarus {

  /*  Rebinds the underlying local function coeff coordinate type to a new type, if the ids have a match
   * For binaryExpr
   * */
  template <template <typename, typename> class Op, typename E1, typename E2, typename OtherType, size_t ID>
  auto rebind(const std::remove_cvref_t<E1>& u, const std::remove_cvref_t<E2>& v,
              Dune::index_constant<ID>&& id = Dune::index_constant<0>()) {
    return Op<decltype(u.rebindClone(OtherType(), std::forward<Dune::index_constant<ID>>(id))),
              decltype(v.rebindClone(OtherType(), std::forward<Dune::index_constant<ID>>(id)))>(
        u.rebindClone(OtherType(), std::forward<Dune::index_constant<ID>>(id)),
        v.rebindClone(OtherType(), std::forward<Dune::index_constant<ID>>(id)));
  }

  /*  Rebinds the underlying local function coeff coordinate type to a new type, if the ids have a match
   * For unaryExpr
   * */
  template <template <typename, typename...> class Op, typename E1, typename OtherType, size_t ID, typename... Args>
  auto rebind(const std::remove_cvref_t<E1>& u, Dune::index_constant<ID>&& id = Dune::index_constant<0>()) {
    return Op<decltype(u.rebindClone(OtherType(), std::forward<Dune::index_constant<ID>>(id))), Args...>(
        u.rebindClone(OtherType(), std::forward<Dune::index_constant<ID>>(id)));
  }

}  // namespace Ikarus
