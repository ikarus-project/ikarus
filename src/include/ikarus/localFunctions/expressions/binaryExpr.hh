//
// Created by Alex on 21.04.2022.
//

#pragma once

#include <ikarus/localFunctions/localFunctionInterface.hh>
#include <ikarus/localFunctions/expressions/unaryExpr.hh>

namespace Ikarus {

  template <typename Op, typename E1, typename E2>
  struct BinaryLocalFunctionExpression : public Ikarus::LocalFunctionInterface<Op> {

    // if the E1 expression is an arithmetic type (int or double) then we store the value otherwise we store a const reference
    using E1StorageType = std::conditional_t<IsNonArithmeticLeafNode<E1>,const E1&,E1>;
    using E2StorageType = std::conditional_t<IsNonArithmeticLeafNode<E2>, const E2&,E2>;

    std::tuple<E1StorageType,E2StorageType> expr;

    using L = E1;
    using R = E2;

    const auto& l()const{
      return std::get<0>(expr);
    }
    const auto& r()const{
      return std::get<1>(expr);
    }

    /* Creates a tuple of all subtype ids, size l or r is not a tuple, tuple_cat may not work.
     * Thus we artifically wrap them inside a tuple  */
    using Ids =decltype(Std::makeNestedTupleFlat(std::make_tuple(std::declval<typename E1::Ids>(),std::declval<typename E2::Ids>() )));

    constexpr BinaryLocalFunctionExpression(E1 const& u,E2 const& v)
      requires  IsLocalFunction<E1,E2>
    : expr(u,v) {}


    static constexpr bool isLeaf = false;
    static constexpr int children = 2;
  };

}  // namespace Ikarus