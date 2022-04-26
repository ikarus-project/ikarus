//
// Created by Alex on 21.04.2022.
//

#pragma once

#include <ikarus/localFunctions/localFunctionInterface.hh>
namespace Ikarus {

  template <typename Op, typename E1, typename E2>
  class BinaryLocalFunctionExpression : public Ikarus::LocalFunctionInterface<Op> {
  protected:
    E1 const& _u;
    E2 const& _v;

  public:
    BinaryLocalFunctionExpression(Ikarus::LocalFunctionInterface<E1> const& u,
                                  Ikarus::LocalFunctionInterface<E2> const& v)
      requires(!std::is_arithmetic_v<E1>)
    : _u(static_cast<E1 const&>(u)), _v(static_cast<E2 const&>(v)) {}
    BinaryLocalFunctionExpression(E1 const& u, Ikarus::LocalFunctionInterface<E2> const& v)
      requires std::is_arithmetic_v<E1>
    : _u(u), _v(static_cast<E2 const&>(v)) {}

    static constexpr bool isLeaf = false;
  };

}  // namespace Ikarus