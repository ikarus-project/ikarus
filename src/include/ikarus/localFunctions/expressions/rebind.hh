//
// Created by lex on 4/30/22.
//

#pragma once
#include <cstddef>
#include <ikarus/utils/traits.hh>
#include <tuple>
namespace Ikarus
{

  /*  Rebinds the underlying local function coeff coordinate type to a new type, if the ids have a match   */
  template <template <typename E1, typename E2> class Op,typename E1, typename E2, typename OtherType,typename IDTuple>
  struct Rebind {
    static auto helper()
    {
      constexpr bool E1mustBeRebound = E1::isLeaf and Std::any_of(IDTuple(),E1::Ids);
      constexpr bool E2mustBeRebound = E2::isLeaf and Std::any_of(IDTuple(),E2::Ids);
        using E1New = typename E1::template Rebind<OtherType>::other;
        using E2New = typename E2::template Rebind<OtherType>::other;
        return Op<std::conditional_t<E1mustBeRebound,E1New,E1>,std::conditional_t<E2mustBeRebound,E2New,E2>>();
    }

    using other = decltype(helper());
  };
}
