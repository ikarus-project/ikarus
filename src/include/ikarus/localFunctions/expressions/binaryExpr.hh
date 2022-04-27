//
// Created by Alex on 21.04.2022.
//

#pragma once

#include <ikarus/localFunctions/localFunctionInterface.hh>
namespace Ikarus {

  template <typename Op, typename E1, typename E2>
  class BinaryLocalFunctionExpression : public Ikarus::LocalFunctionInterface<Op> {
  protected:
    //if the E1 expression is an arithmetic type (int or double) then we store the value otherwise we store a const reference
    E1 const& l; //std::conditional_t<std::is_arithmetic_v<E1>,E1,E1 const&>
    E2 const& r;

  public:
    BinaryLocalFunctionExpression(Ikarus::LocalFunctionInterface<E1> const& u,
                                  Ikarus::LocalFunctionInterface<E2> const& v)
      requires(!std::is_arithmetic_v<E1>)
    : l(static_cast<E1 const&>(u)), r(static_cast<E2 const&>(v)) {}

    BinaryLocalFunctionExpression(E1 const& u, Ikarus::LocalFunctionInterface<E2> const& v)
      requires std::is_arithmetic_v<E1>
    : l(u), r(static_cast<E2 const&>(v)) {}

    static constexpr bool isLeaf = false;
  };

}  // namespace Ikarus