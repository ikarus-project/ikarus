//
// Created by lex on 4/30/22.
//

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
  template <template <typename> class Op, typename E1, typename OtherType, size_t ID>
  auto rebind(const std::remove_cvref_t<E1>& u, Dune::index_constant<ID>&& id = Dune::index_constant<0>()) {
    return Op<decltype(u.rebindClone(OtherType(), std::forward<Dune::index_constant<ID>>(id)))>(
        u.rebindClone(OtherType(), std::forward<Dune::index_constant<ID>>(id)));
  }

}  // namespace Ikarus
