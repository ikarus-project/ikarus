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
//      constexpr bool e1ReboundNeeded = Std::hasType<Dune::index_constant<ID>,typename E1Raw::Ids>::value;
      constexpr bool e2isLeaf = E2Raw::isLeaf;
//      constexpr bool e2ReboundNeeded = Std::hasType<Dune::index_constant<ID>,typename E2Raw::Ids>::value;

      if constexpr (e1isLeaf) {
        using E1New                    = typename E1Raw::template Rebind<OtherType>::other;
        if constexpr (e2isLeaf)
        {
          using E2New                    = typename E2Raw::template Rebind<OtherType>::other;
          return Op<E1New, E2New>(u.rebindClone(OtherType(),std::forward<Dune::index_constant<ID>>(id)),
                                  v.rebindClone(OtherType(),std::forward<Dune::index_constant<ID>>(id)));
        }
        else //FIXME simplify all here with decltype
          return Op<E1New, E2Raw>(u.rebindClone(OtherType(),std::forward<Dune::index_constant<ID>>(id)),
                                  v.rebindClone(OtherType(),std::forward<Dune::index_constant<ID>>(id)));

      } else if constexpr (not e1isLeaf) {
        if constexpr (e2isLeaf)
        {
          using E2New                    = typename E2Raw::template Rebind<OtherType>::other;
          return Op<E1Raw, E2New>(u.rebindClone(OtherType(),std::forward<Dune::index_constant<ID>>(id)),
                                  v.rebindClone(OtherType()),std::forward<Dune::index_constant<ID>>(id));
        }
        else
          return Op<E1Raw, E2Raw>(u.rebindClone(OtherType(),std::forward<Dune::index_constant<ID>>(id)),
                                  v.rebindClone(OtherType()),std::forward<Dune::index_constant<ID>>(id));
      }

    }

    /*  Rebinds the underlying local function coeff coordinate type to a new type, if the ids have a match
   * For unaryExpr
   * */   template <template <typename > class Op, typename E1, typename OtherType,
              size_t ID>
     auto rebind(const std::remove_cvref_t<E1>& u,Dune::index_constant<ID>&& id=Dune::index_constant<0>()) {
      constexpr bool e1isLeaf = std::remove_cvref_t<E1>::isLeaf;
      if constexpr (e1isLeaf) {
        using E1New                    = typename std::remove_cvref_t<E1>::template Rebind<OtherType>::other;
        return Op<E1New>(u.rebindClone(OtherType(), std::forward<Dune::index_constant<ID>>(id)));
      } else
        return Op<std::remove_cvref_t<E1>>(u.rebindClone(OtherType(),std::forward<Dune::index_constant<ID>>(id)));

    }

}
