//
// Created by Alex on 21.04.2022.
//

#pragma once

#include <ikarus/localFunctions/localFunctionInterface.hh>
namespace Ikarus
{

template <typename Op, typename E1>
struct UnaryLocalFunctionExpression : public Ikarus::LocalFunctionInterface<Op> {

  using E1StorageType = std::conditional_t<IsNonArithmeticLeafNode<E1>, const E1 & ,E1>;
  std::tuple<E1StorageType> expr;

  const auto& m()const{
    return std::get<0>(expr);
  }

  using Ids =decltype(Std::makeNestedTupleFlat(std::make_tuple(std::declval<typename E1::Ids>() )));

  constexpr UnaryLocalFunctionExpression(E1 const& u)
  requires IsLocalFunction<E1>
  : expr(u) {  }

  static constexpr bool isLeaf = false;
  static constexpr int children = 1;

};

}