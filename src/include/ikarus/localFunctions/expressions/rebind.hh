//
// Created by lex on 4/30/22.
//

#pragma once
#include <cstddef>
#include <ikarus/utils/traits.hh>
#include <tuple>
namespace Ikarus
{

  /*  Rebinds the underlying local function coeff coordinate type to a new type, if the ids have a match
   * For binaryExpr
   * */   template <template <typename, typename > class Op, typename E1, typename E2, typename OtherType,
            size_t ID>
     auto rebind(const std::remove_cvref_t<E1>& u, const std::remove_cvref_t<E2>& v,Dune::index_constant<ID>&& id=Dune::index_constant<0>()) {
      using E1Raw = std::remove_cvref_t<E1>;
      using E2Raw = std::remove_cvref_t<E2>;
      constexpr bool e1isLeaf = E1Raw::isLeaf;
      constexpr bool e1ReboundNeeded = Std::hasType<Dune::index_constant<ID>,typename E1Raw::Ids>::value;
      constexpr bool e2isLeaf = E2Raw::isLeaf;
      constexpr bool e2ReboundNeeded = Std::hasType<Dune::index_constant<ID>,typename E2Raw::Ids>::value;

      using E1New                    = typename E1Raw::template Rebind<OtherType>::other;
      using E2New                    = typename E2Raw::template Rebind<OtherType>::other;
      if constexpr (e1isLeaf) {
        if constexpr (e1ReboundNeeded and not e2ReboundNeeded) {
          return Op<E1New, E2Raw>(u.rebindClone(OtherType()), v.clone());
        }else if constexpr (e2ReboundNeeded)
        {
          if constexpr (e2isLeaf)
            return Op<E1New, E2New>(u.rebindClone(OtherType()), v.rebindClone(OtherType()));
          else
            return Op<E1New, E2New>(u.rebindClone(OtherType()), v.rebindClone(id,OtherType()));
        }
      } else if constexpr (not e1isLeaf) {
        if constexpr (e1ReboundNeeded and not e2ReboundNeeded) {
          return Op<E1New, std::remove_cvref_t<E2>>(u.rebindClone(id,OtherType()), v.clone());
        }else if constexpr (e2ReboundNeeded)
        {
          if constexpr (e2isLeaf)
            return Op<E1New, E2New>(u.rebindClone(id,OtherType()), v.rebindClone(OtherType()));
          else
            return Op<E1New, E2New>(u.rebindClone(id,OtherType()), v.rebindClone(id,OtherType()));
        }
      }

    }

    /*  Rebinds the underlying local function coeff coordinate type to a new type, if the ids have a match
   * For unaryExpr
   * */   template <template <typename > class Op, typename E1, typename OtherType,
              size_t ID>
     auto rebind(const std::remove_cvref_t<E1>& u,Dune::index_constant<ID>&& id=Dune::index_constant<0>()) {
      constexpr bool e1isLeaf = std::remove_cvref_t<E1>::isLeaf;
      using E1New                    = typename std::remove_cvref_t<E1>::template Rebind<OtherType>::other;
      if constexpr (e1isLeaf)
          return Op<E1New>(u.rebindClone(OtherType(),std::forward<Dune::index_constant<ID>>(id)));
      else
        return Op<std::remove_cvref_t<E1>>(u.rebindClone(OtherType(),std::forward<Dune::index_constant<ID>>(id)));

    }

}
