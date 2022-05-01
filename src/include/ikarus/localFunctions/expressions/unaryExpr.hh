//
// Created by Alex on 21.04.2022.
//

#pragma once

#include <ikarus/localFunctions/localFunctionInterface.hh>
#include "rebind.hh"
namespace Ikarus
{

template <template<typename> class Op, typename E1>
struct UnaryLocalFunctionExpression : public Ikarus::LocalFunctionInterface<Op<E1>> {

  std::tuple<E1> expr;

  using M = E1;

  const E1& m() const{
    return std::get<0>(expr);
  }

   E1& m(){
    return std::get<0>(expr);
  }

  using Ids =decltype(Std::makeNestedTupleFlat(std::make_tuple(std::declval<typename std::remove_cvref_t<E1>::Ids>() )));

  auto clone ()const{
    return Op<decltype(m().clone())>(m().clone());
  }

  template <typename OtherType,size_t ID=0 >
  auto rebindClone( OtherType&& ,Dune::index_constant<ID>&& id =Dune::index_constant<0>()) const
  {
    return rebind<Op,E1,OtherType>(m(),Dune::index_constant<ID>() );
  }

  constexpr explicit UnaryLocalFunctionExpression(E1 && u)
  requires IsLocalFunction<E1>
  : expr(std::forward<E1>(u)) {  }

  static constexpr bool isLeaf = false;
  static constexpr int children = 1;

};

}