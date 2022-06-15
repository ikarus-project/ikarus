//
// Created by lex on 6/14/22.
//

#pragma once

#include <ikarus/localFunctions/expressions/scalarUnaryExpressions/scalarUnaryExpression.hh>

namespace Ikarus {

  struct SqrtFunc {
    template <typename ScalarType>
    static auto value(const ScalarType& v) {
      return sqrt(v);
    }

    template <typename ScalarType>
    static auto derivative(const ScalarType& v) {
      return ScalarType(1.0) / (ScalarType(2.0) * sqrt(v));
    }

    template <typename ScalarType>
    static auto secondDerivative(const ScalarType& v) {
      using std::pow;
      return ScalarType(-1.0) / (ScalarType(4.0) * pow(v, 3.0 / 2.0));
    }

    template <typename ScalarType>
    static auto thirdDerivative(const ScalarType& v) {
      using std::pow;
      return ScalarType(3.0) / (ScalarType(8.0) * pow(v, 5.0 / 2.0));
    }
  };
  template <typename E1>
  requires IsLocalFunction<E1>
  constexpr auto sqrt(E1&& u) {
    static_assert(std::remove_cvref_t<E1>::valueSize == 1,
                  "Sqrt expression only defined for scalar valued local functions.");
    return ScalarUnaryExpression<E1, SqrtFunc>(std::forward<E1>(u));
  }
}  // namespace Ikarus
